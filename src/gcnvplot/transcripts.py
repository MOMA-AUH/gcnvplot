"""Transcript annotation loading and SQLite indexing."""

from __future__ import annotations

import sqlite3
from pathlib import Path

from .models import TranscriptAnnotation, TranscriptExon
from .utils import open_text


def parse_gtf_attributes(text: str) -> dict[str, str]:
    """Parse a GTF attribute column into a dictionary."""
    attributes: dict[str, str] = {}
    for item in text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        key, value = item.split(" ", 1)
        attributes[key] = value.strip().strip('"')
    return attributes


def transcript_id_matches(value: str | None, target: str) -> bool:
    """Return True when two transcript IDs match with or without version suffixes."""
    if value is None:
        return False
    if value == target:
        return True
    return value.split(".", 1)[0] == target.split(".", 1)[0]


def _resolve_transcript_annotation(
    annotations: dict[str, TranscriptAnnotation],
    transcript_id: str,
    source: Path,
) -> TranscriptAnnotation:
    """Resolve one transcript by exact ID or unversioned ID when unambiguous."""
    exact = annotations.get(transcript_id)
    if exact is not None:
        return exact

    root_id = transcript_id.split(".", 1)[0]
    matches = [
        annotation
        for candidate_id, annotation in annotations.items()
        if candidate_id.split(".", 1)[0] == root_id
    ]
    if not matches:
        raise ValueError(f"{source}: transcript {transcript_id} not found or has no exons")
    if len(matches) > 1:
        versions = ", ".join(sorted(annotation.transcript_id for annotation in matches))
        raise ValueError(f"{source}: transcript {transcript_id} is ambiguous; matches: {versions}")
    return matches[0]


def load_transcript_annotations(gtf_path: Path) -> dict[str, TranscriptAnnotation]:
    """Load transcript annotations from a GTF keyed by exact transcript ID."""
    records: dict[str, dict[str, object]] = {}

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqname, _source, feature, start_text, end_text, _score, feature_strand, _frame, attributes_text = fields
            attributes = parse_gtf_attributes(attributes_text)
            transcript_id = attributes.get("transcript_id")
            if transcript_id is None:
                continue

            record = records.setdefault(
                transcript_id,
                {
                    "contig": None,
                    "strand": None,
                    "gene_name": None,
                    "gene_id": None,
                    "exon_coords": [],
                },
            )

            feature_contig = seqname
            feature_strand = feature_strand
            contig = record["contig"]
            strand = record["strand"]
            if contig is None:
                record["contig"] = feature_contig
            elif contig != feature_contig:
                raise ValueError(f"{gtf_path}: transcript {transcript_id} spans multiple contigs")

            if strand is None:
                record["strand"] = feature_strand
            elif strand != feature_strand:
                raise ValueError(f"{gtf_path}: transcript {transcript_id} has inconsistent strands")

            if record["gene_name"] is None:
                record["gene_name"] = attributes.get("gene_name") or attributes.get("gene_id")
            if record["gene_id"] is None:
                record["gene_id"] = attributes.get("gene_id")

            if feature == "exon":
                exon_coords = record["exon_coords"]
                assert isinstance(exon_coords, list)
                exon_coords.append((int(start_text), int(end_text)))

    annotations: dict[str, TranscriptAnnotation] = {}
    for transcript_id, record in records.items():
        exon_coords = record["exon_coords"]
        assert isinstance(exon_coords, list)
        if not exon_coords:
            continue

        contig = record["contig"]
        strand = record["strand"]
        if contig is None or strand is None:
            raise ValueError(f"{gtf_path}: transcript {transcript_id} is missing contig or strand information")

        ordered_coords = sorted(exon_coords, key=lambda item: item[0], reverse=strand == "-")
        exons = [
            TranscriptExon(number=index + 1, start=start, end=end)
            for index, (start, end) in enumerate(ordered_coords)
        ]
        gene_name = record["gene_name"]
        gene_id = record["gene_id"]
        transcript_gene_name = str(gene_name or gene_id or transcript_id)
        start = min(exon.start for exon in exons)
        end = max(exon.end for exon in exons)
        annotations[transcript_id] = TranscriptAnnotation(
            transcript_id=transcript_id,
            gene_name=transcript_gene_name,
            contig=str(contig),
            strand=str(strand),
            start=start,
            end=end,
            exons=exons,
        )

    return annotations


def load_transcript_annotation(gtf_path: Path, transcript_id: str) -> TranscriptAnnotation:
    """Load one transcript annotation from a GTF."""
    annotations = load_transcript_annotations(gtf_path)
    return _resolve_transcript_annotation(annotations, transcript_id, gtf_path)


def index_gtf(gtf_path: Path, output: Path) -> tuple[int, int]:
    """Index transcript annotations from a GTF into a SQLite database."""
    annotations = load_transcript_annotations(gtf_path)
    transcript_count = len(annotations)
    exon_count = sum(len(annotation.exons) for annotation in annotations.values())

    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        output.unlink()

    with sqlite3.connect(output) as connection:
        connection.executescript(
            """
            CREATE TABLE metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );
            CREATE TABLE transcripts (
                transcript_id TEXT PRIMARY KEY,
                transcript_id_root TEXT NOT NULL,
                gene_name TEXT NOT NULL,
                contig TEXT NOT NULL,
                strand TEXT NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER NOT NULL
            );
            CREATE TABLE exons (
                transcript_id TEXT NOT NULL,
                exon_number INTEGER NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER NOT NULL,
                PRIMARY KEY (transcript_id, exon_number),
                FOREIGN KEY (transcript_id) REFERENCES transcripts(transcript_id) ON DELETE CASCADE
            );
            CREATE INDEX idx_transcripts_root ON transcripts (transcript_id_root);
            CREATE INDEX idx_exons_transcript ON exons (transcript_id);
            """
        )
        connection.executemany(
            "INSERT INTO metadata(key, value) VALUES(?, ?)",
            [
                ("source_gtf", str(gtf_path)),
                ("schema_version", "1"),
            ],
        )
        connection.executemany(
            """
            INSERT INTO transcripts(
                transcript_id,
                transcript_id_root,
                gene_name,
                contig,
                strand,
                start,
                end
            ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    annotation.transcript_id,
                    annotation.transcript_id.split(".", 1)[0],
                    annotation.gene_name,
                    annotation.contig,
                    annotation.strand,
                    annotation.start,
                    annotation.end,
                )
                for annotation in annotations.values()
            ],
        )
        connection.executemany(
            "INSERT INTO exons(transcript_id, exon_number, start, end) VALUES (?, ?, ?, ?)",
            [
                (annotation.transcript_id, exon.number, exon.start, exon.end)
                for annotation in annotations.values()
                for exon in annotation.exons
            ],
        )

    return transcript_count, exon_count


class TranscriptIndex:
    """Reusable SQLite-backed transcript annotation index."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self._connection = sqlite3.connect(self.path)
        self._connection.row_factory = sqlite3.Row

    def close(self) -> None:
        """Close the underlying SQLite connection."""
        self._connection.close()

    def __enter__(self) -> "TranscriptIndex":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    def get(self, transcript_id: str) -> TranscriptAnnotation:
        """Return one transcript annotation by exact ID or unversioned ID when unambiguous."""
        transcript_row = self._connection.execute(
            """
            SELECT transcript_id, gene_name, contig, strand, start, end
            FROM transcripts
            WHERE transcript_id = ?
            """,
            (transcript_id,),
        ).fetchone()
        if transcript_row is None:
            root_id = transcript_id.split(".", 1)[0]
            transcript_rows = self._connection.execute(
                """
                SELECT transcript_id, gene_name, contig, strand, start, end
                FROM transcripts
                WHERE transcript_id_root = ?
                ORDER BY transcript_id
                """,
                (root_id,),
            ).fetchall()
            if not transcript_rows:
                raise ValueError(f"{self.path}: transcript {transcript_id} not found")
            if len(transcript_rows) > 1:
                versions = ", ".join(str(row["transcript_id"]) for row in transcript_rows)
                raise ValueError(f"{self.path}: transcript {transcript_id} is ambiguous; matches: {versions}")
            transcript_row = transcript_rows[0]

        exon_rows = self._connection.execute(
            """
            SELECT exon_number, start, end
            FROM exons
            WHERE transcript_id = ?
            ORDER BY exon_number
            """,
            (transcript_row["transcript_id"],),
        ).fetchall()
        if not exon_rows:
            raise ValueError(f"{self.path}: transcript {transcript_row['transcript_id']} has no indexed exons")

        exons = [
            TranscriptExon(number=int(row["exon_number"]), start=int(row["start"]), end=int(row["end"]))
            for row in exon_rows
        ]
        return TranscriptAnnotation(
            transcript_id=str(transcript_row["transcript_id"]),
            gene_name=str(transcript_row["gene_name"]),
            contig=str(transcript_row["contig"]),
            strand=str(transcript_row["strand"]),
            start=int(transcript_row["start"]),
            end=int(transcript_row["end"]),
            exons=exons,
        )


def load_transcript_annotation_from_db(db_path: Path, transcript_id: str) -> TranscriptAnnotation:
    """Load one transcript annotation from a SQLite transcript database."""
    with TranscriptIndex(db_path) as transcript_index:
        return transcript_index.get(transcript_id)
