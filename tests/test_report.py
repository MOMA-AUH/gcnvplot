"""Tests for transcript metadata, segment parsing, and HTML reports."""

from __future__ import annotations

import gzip
import sqlite3
from pathlib import Path

import pytest

import gcnvplot
import gcnvplot.cli as cli
import gcnvplot.plotting as plotting


GTF_LINES = [
    'NC_000001.11\ttest\tgene\t100\t250\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id ""; db_xref "HGNC:HGNC:1";',
    'NC_000001.11\ttest\ttranscript\t100\t250\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "CUSTOM_TX.1"; db_xref "GenBank:CUSTOM_TX.1"; tag "MANE Select";',
    'NC_000001.11\ttest\texon\t100\t150\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "CUSTOM_TX.1";',
    'NC_000001.11\ttest\texon\t200\t250\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; transcript_id "CUSTOM_TX.1";',
]


def _write_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    read_counts = tmp_path / "sample.tsv"
    read_counts.write_text(
        "CONTIG\tSTART\tEND\tCOUNT\n"
        "chr1\t100\t150\t10\n"
        "chr1\t200\t250\t5\n",
        encoding="utf-8",
    )
    background = tmp_path / "background.tsv"
    background.write_text(
        "# normalization=median-of-ratios\n"
        "# samples=12\n"
        "# lower_percentile=5\n"
        "# upper_percentile=95\n"
        "CONTIG\tSTART\tEND\tBASELINE_MEDIAN\tN\tBG_NORM_MEAN\t"
        "BG_NORM_MEDIAN\tBG_NORM_SD\tBG_NORM_P5\tBG_NORM_P95\n"
        "chr1\t100\t150\t10\t12\t10\t10\t1\t8\t12\n"
        "chr1\t200\t250\t10\t12\t10\t10\t1\t8\t12\n",
        encoding="utf-8",
    )
    gtf = tmp_path / "annotations.gtf"
    gtf.write_text(
        "#!genome-build GRCh38.p14\n" + "\n".join(GTF_LINES) + "\n",
        encoding="utf-8",
    )
    transcript_db = tmp_path / "annotations.sqlite"
    assert gcnvplot.index_gtf(gtf, transcript_db) == (1, 2)
    return read_counts, background, gtf, transcript_db


def _write_genes(
    path: Path,
    *,
    gene_name: str = "GENE1",
    stable_id: str = "HGNC:1",
    transcript: str = "",
) -> None:
    path.write_text(
        "gene_name\tstable_id\ttranscript\n"
        f"{gene_name}\t{stable_id}\t{transcript}\n",
        encoding="utf-8",
    )


def _write_segments(
    path: Path,
    *,
    filters: str = "PASS",
    sample: str = "SAMPLE",
) -> None:
    text = (
        "##fileformat=VCFv4.2\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        f"chr1\t120\tcall1\tN\t<DEL>,<DUP>\t100\t{filters}\t"
        "END=220\tGT:CN\t1:1\n"
    )
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def test_index_stores_generic_transcript_gene_ids_and_mane(tmp_path: Path) -> None:
    _counts, _background, _gtf, transcript_db = _write_inputs(tmp_path)
    with gcnvplot.TranscriptIndex(transcript_db) as index:
        annotation = index.resolve_gene(
            "gene1",
            "HGNC:1",
            available_contigs={"chr1"},
        )
        explicit = index.resolve_gene(
            "GENE1",
            "G1",
            "CUSTOM_TX",
            available_contigs={"chr1"},
        )
    assert annotation.transcript_id == "CUSTOM_TX.1"
    assert annotation.is_mane_select
    assert annotation.stable_ids == ("G1", "HGNC:1")
    assert explicit == annotation


def test_index_has_one_current_schema_without_compatibility_tables(
    tmp_path: Path,
) -> None:
    _counts, _background, _gtf, transcript_db = _write_inputs(tmp_path)
    with sqlite3.connect(transcript_db) as connection:
        schema_version = connection.execute(
            "SELECT value FROM metadata WHERE key = 'schema_version'"
        ).fetchone()[0]
        columns = {
            row[1]
            for row in connection.execute("PRAGMA table_info(transcripts)")
        }
    assert schema_version == "2"
    assert "is_mane_select" in columns


def test_gene_list_accepts_arbitrary_transcript_identifier(tmp_path: Path) -> None:
    genes = tmp_path / "genes.tsv"
    _write_genes(genes, transcript="CUSTOM_TX.1")
    assert gcnvplot.read_gene_requests(genes) == [
        gcnvplot.GeneRequest("GENE1", "HGNC:1", "CUSTOM_TX.1")
    ]


@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_parse_segments_and_native_overlay(tmp_path: Path, suffix: str) -> None:
    segments_path = tmp_path / f"segments{suffix}"
    _write_segments(segments_path)
    sample, segments = gcnvplot.parse_segments_vcf(segments_path)
    assert sample == "SAMPLE"
    assert len(segments) == 1
    assert segments[0].svtype == "DEL"
    assert segments[0].as_overlay("chr1").kind == "deletion"


def test_multisample_segments_require_explicit_selection(tmp_path: Path) -> None:
    path = tmp_path / "multi.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
        "chr1\t120\tcall\tN\t<DEL>\t100\tPASS\tEND=220\tGT:CN\t1:1\t0:2\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="specify --sample-name"):
        gcnvplot.parse_segments_vcf(path)
    sample, segments = gcnvplot.parse_segments_vcf(path, sample_name="S2")
    assert sample == "S2"
    assert segments == []


def test_report_cli_has_one_sample_name_option(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["report", "--help"])
    assert exc_info.value.code == 0
    help_text = capsys.readouterr().out
    assert "--sample-name" in help_text
    assert "--vcf-sample" not in help_text


def test_filtered_segments_are_explicitly_opted_in(tmp_path: Path) -> None:
    segments_path = tmp_path / "segments.vcf"
    _write_segments(segments_path, filters="LowQS")
    assert gcnvplot.parse_segments_vcf(segments_path)[1] == []
    _sample, segments = gcnvplot.parse_segments_vcf(
        segments_path,
        include_filtered=True,
    )
    assert segments[0].filter_status == "LowQS"
    assert "FILTER=LowQS" in segments[0].label


def test_prepared_data_normalizes_only_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    counts, background, _gtf, transcript_db = _write_inputs(tmp_path)
    original = plotting.normalized_counts
    calls = 0

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(plotting, "normalized_counts", counted)
    prepared = gcnvplot.prepare_plot_data(counts, background)
    with gcnvplot.TranscriptIndex(transcript_db) as index:
        annotation = index.resolve_gene(
            "GENE1",
            "HGNC:1",
            available_contigs={"chr1"},
        )
    for _ in range(2):
        plotting.render_prepared_plot_svg(
            prepared,
            ("chr1", 100, 250),
            transcript=annotation,
        )
    assert calls == 1


def test_transcript_plot_adapts_refseq_contig_to_count_contig(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    counts, background, _gtf, transcript_db = _write_inputs(tmp_path)
    output = tmp_path / "plot.svg"

    assert cli.main(
        [
            "plot",
            "--read-counts",
            str(counts),
            "--background",
            str(background),
            "--transcript",
            "CUSTOM_TX.1",
            "--transcript-db",
            str(transcript_db),
            "--output",
            str(output),
        ]
    ) == 0
    assert "chr1:100-250" in output.read_text(encoding="utf-8")
    capsys.readouterr()

    svg = gcnvplot.render_plot_svg(
        counts,
        background,
        transcript_id="CUSTOM_TX.1",
        transcript_index=transcript_db,
    )
    assert "chr1:100-250" in svg


def test_cli_builds_self_contained_report_with_native_cnv_overlay(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    counts, background, _gtf, transcript_db = _write_inputs(tmp_path)
    genes = tmp_path / "genes.tsv"
    _write_genes(genes)
    segments = tmp_path / "segments.vcf"
    report_sample = 'SAMPLE </script> "one"'
    _write_segments(segments, sample=report_sample)
    output = tmp_path / "report.html"

    assert cli.main(
        [
            "report",
            "--read-counts",
            str(counts),
            "--background",
            str(background),
            "--genes",
            str(genes),
            "--transcript-db",
            str(transcript_db),
            "--segments-vcf",
            str(segments),
            "--sample-name",
            report_sample,
            "--output",
            str(output),
        ]
    ) == 0

    page = output.read_text(encoding="utf-8")
    assert "<template" in page
    assert '<ul id="gene-suggestions" role="listbox"' in page
    assert "const geneNames = Object.values(genes)" in page
    assert "gene-options" not in page
    assert 'class="overlay-band overlay-deletion"' in page
    assert ">Highlight</text>" in page
    assert "chr1:120-220" in page
    assert ">Exons</text>" in page
    assert ">1-2</text>" in page
    assert "gcnvplot " in page
    assert "GRCh38.p14" in page
    assert "12 background samples" not in page
    assert "SAMPLE &lt;/script&gt; &quot;one&quot;" in page
    assert "SAMPLE </script>" not in page
    assert "https://" not in page
    assert capsys.readouterr().out == (
        "Report genes: 1\n"
        "Overlapping CNV segments: 1\n"
        f"Wrote: {output}\n"
    )


def test_report_failure_does_not_replace_existing_output(tmp_path: Path) -> None:
    counts, background, _gtf, transcript_db = _write_inputs(tmp_path)
    genes = tmp_path / "genes.tsv"
    _write_genes(genes, gene_name="WRONG")
    output = tmp_path / "report.html"
    output.write_text("keep me", encoding="utf-8")

    with pytest.raises(ValueError, match="Report gene validation failed"):
        gcnvplot.build_report(
            counts,
            background,
            genes,
            transcript_db,
            output,
        )
    assert output.read_text(encoding="utf-8") == "keep me"
