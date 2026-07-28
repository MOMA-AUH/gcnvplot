"""Parsing of GATK GermlineCNVCaller segment VCFs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .models import PlotOverlay, TranscriptAnnotation
from .transcripts import contig_aliases
from .utils import open_text


@dataclass(frozen=True)
class CnvSegment:
    """One called deletion or duplication segment."""

    contig: str
    start: int
    end: int
    svtype: str
    copy_number: int | None
    record_id: str | None
    sample: str
    filter_status: str

    @property
    def label(self) -> str:
        """Return a compact user-facing description."""
        copy_number = "" if self.copy_number is None else f" · CN={self.copy_number}"
        filtered = (
            ""
            if self.filter_status in {".", "PASS"}
            else f" · FILTER={self.filter_status}"
        )
        return (
            f"{self.svtype} · {self.contig}:{self.start:,}-{self.end:,}"
            f"{copy_number}{filtered}"
        )

    def as_overlay(self, contig: str | None = None) -> PlotOverlay:
        """Return the native plot overlay for this segment."""
        kind = "deletion" if self.svtype == "DEL" else "duplication"
        return PlotOverlay(
            (contig or self.contig, self.start, self.end),
            kind,
            self.label,
        )


def _parse_info(text: str) -> dict[str, str | None]:
    info: dict[str, str | None] = {}
    if text == ".":
        return info
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = None
    return info


def _integer(value: str | None) -> int | None:
    if value is None or value in {"", "."}:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _called_svtype(
    alt_text: str,
    genotype: str | None,
    copy_number: int | None,
) -> str | None:
    alternatives = alt_text.split(",") if alt_text not in {"", "."} else []
    structural = {
        index: allele[1:-1].upper()
        for index, allele in enumerate(alternatives, start=1)
        if allele.upper() in {"<DEL>", "<DUP>"}
    }
    if not structural:
        return None

    if genotype not in {None, "", ".", "./.", ".|."}:
        alleles = [
            int(value)
            for value in re.split(r"[|/]", genotype)
            if value.isdigit()
        ]
        called = {
            structural[index]
            for index in alleles
            if index > 0 and index in structural
        }
        return called.pop() if len(called) == 1 else None

    # Some GATK segment outputs omit GT for a called record. Only accept this
    # when the record has a copy-number value and one unambiguous CNV type;
    # never infer the type by comparing CN with a diploid baseline.
    types = set(structural.values())
    if copy_number is not None and len(types) == 1:
        return next(iter(types))
    return None


def parse_segments_vcf(
    path: Path,
    *,
    sample_name: str | None = None,
    include_filtered: bool = False,
) -> tuple[str, list[CnvSegment]]:
    """Read called DEL/DUP segments for one sample from a VCF or VCF.GZ."""
    selected_sample: str | None = None
    selected_column: int | None = None
    segments: list[CnvSegment] = []

    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                if not samples:
                    raise ValueError(f"{path}: segments VCF has no sample column")
                if len(samples) != len(set(samples)):
                    raise ValueError(f"{path}: segments VCF has duplicate sample names")
                if sample_name is None:
                    if len(samples) != 1:
                        raise ValueError(
                            f"{path}: segments VCF contains {len(samples)} samples; "
                            "specify --sample-name"
                        )
                    selected_sample = samples[0]
                elif sample_name not in samples:
                    raise ValueError(
                        f"{path}: VCF sample {sample_name!r} not found; "
                        f"available: {', '.join(samples)}"
                    )
                else:
                    selected_sample = sample_name
                selected_column = 9 + samples.index(selected_sample)
                continue
            if line.startswith("#"):
                continue
            if not line.strip():
                continue
            if selected_column is None or selected_sample is None:
                raise ValueError(f"{path}:{line_number}: VCF header is missing #CHROM")

            fields = line.rstrip("\n").split("\t")
            if len(fields) <= selected_column or len(fields) < 10:
                raise ValueError(f"{path}:{line_number}: malformed VCF record")
            contig, position, record_id, _ref, alt, _qual, filters, info_text = fields[:8]
            if not include_filtered and filters not in {".", "PASS"}:
                continue
            try:
                start = int(position)
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: VCF POS is not an integer"
                ) from error
            info = _parse_info(info_text)
            end = _integer(info.get("END"))
            if end is None:
                raise ValueError(
                    f"{path}:{line_number}: segment is missing integer INFO/END"
                )
            if end < start:
                raise ValueError(f"{path}:{line_number}: segment END is before POS")

            format_keys = fields[8].split(":")
            sample_values = fields[selected_column].split(":")
            if len(sample_values) > len(format_keys):
                raise ValueError(
                    f"{path}:{line_number}: sample has more values than FORMAT fields"
                )
            sample_data = dict(zip(format_keys, sample_values))
            copy_number = _integer(sample_data.get("CN"))
            svtype = _called_svtype(alt, sample_data.get("GT"), copy_number)
            if svtype is None:
                continue
            segments.append(
                CnvSegment(
                    contig=contig,
                    start=start,
                    end=end,
                    svtype=svtype,
                    copy_number=copy_number,
                    record_id=None if record_id == "." else record_id,
                    sample=selected_sample,
                    filter_status=filters,
                )
            )

    if selected_sample is None:
        raise ValueError(f"{path}: VCF header is missing #CHROM")
    return selected_sample, segments


def segments_for_transcript(
    segments: list[CnvSegment],
    transcript: TranscriptAnnotation,
) -> list[CnvSegment]:
    """Return segments overlapping a transcript using common contig aliases."""
    transcript_aliases = set(contig_aliases(transcript.contig))
    return [
        segment
        for segment in segments
        if transcript_aliases.intersection(contig_aliases(segment.contig))
        and segment.start <= transcript.end
        and segment.end >= transcript.start
    ]
