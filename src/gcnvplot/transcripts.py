"""Transcript annotation loading and SQLite indexing."""

from __future__ import annotations

import os
import sqlite3
import tempfile
from collections import defaultdict
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
from urllib.parse import quote

from .models import TranscriptAnnotation, TranscriptExon
from .utils import open_text

_REFSEQ_CHROMOSOMES = {
    **{f"NC_{number:06d}": str(number) for number in range(1, 23)},
    "NC_000023": "X",
    "NC_000024": "Y",
    "NC_012920": "MT",
}


def _parse_gtf_attributes_multi(text: str) -> dict[str, list[str]]:
    """Parse GTF attributes while preserving repeated keys."""
    attributes: dict[str, list[str]] = defaultdict(list)
    for item in text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(None, 1)
        if len(parts) != 2:
            continue
        key, value = parts
        attributes[key].append(value.strip().strip('"'))
    return dict(attributes)


def parse_gtf_attributes(text: str) -> dict[str, str]:
    """Parse a GTF attribute column into a dictionary."""
    return {
        key: values[-1]
        for key, values in _parse_gtf_attributes_multi(text).items()
    }


def _first(attributes: dict[str, list[str]], *keys: str) -> str | None:
    for key in keys:
        for value in attributes.get(key, []):
            value = value.strip()
            if value:
                return value
    return None


def normalize_gene_identifier(value: str) -> str:
    """Normalize a gene identifier without imposing one identifier system."""
    identifier = value.strip()
    if identifier.upper().startswith("HGNC:HGNC:"):
        identifier = identifier[5:]
    return identifier


def _gene_identifiers(
    attributes: dict[str, list[str]],
    *,
    include_cross_references: bool,
) -> set[str]:
    identifiers: set[str] = set()
    keys = ["gene_id", "hgnc_id", "HGNC_ID", "hgnc"]
    if include_cross_references:
        keys.extend(("db_xref", "Dbxref"))
    for key in keys:
        for value in attributes.get(key, []):
            identifier = normalize_gene_identifier(value)
            if identifier:
                identifiers.add(identifier)
    return identifiers


def _is_mane_select(attributes: dict[str, list[str]]) -> bool:
    values = attributes.get("tag", []) + attributes.get("MANE_status", [])
    return any(value.casefold().replace("_", " ") == "mane select" for value in values)


def contig_aliases(contig: str) -> tuple[str, ...]:
    """Return common aliases for a chromosome or RefSeq contig."""
    aliases = [contig]
    root = contig.split(".", 1)[0]
    chromosome = _REFSEQ_CHROMOSOMES.get(root)
    if chromosome is not None:
        aliases.extend((f"chr{chromosome}", chromosome))
        if chromosome == "MT":
            aliases.extend(("chrM", "M"))
    elif contig.startswith("chr"):
        aliases.append(contig[3:])
    else:
        aliases.append(f"chr{contig}")
    return tuple(dict.fromkeys(aliases))


def adapt_transcript_contig(
    annotation: TranscriptAnnotation,
    available_contigs: set[str],
    *,
    label: str | None = None,
) -> TranscriptAnnotation:
    """Adapt an annotation contig to an equivalent contig present in plot data."""
    if annotation.contig in available_contigs:
        return annotation
    for alias in contig_aliases(annotation.contig):
        if alias in available_contigs:
            return replace(annotation, contig=alias)
    available = ", ".join(sorted(available_contigs)[:8])
    subject = label or annotation.transcript_id
    raise ValueError(
        f"{subject}: annotation contig {annotation.contig!r} does not match "
        f"the read-count/background contigs ({available})"
    )


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
    gene_metadata: dict[str, dict[str, object]] = {}

    with open_text(gtf_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqname, _source, feature, start_text, end_text, _score, feature_strand, _frame, attributes_text = fields
            attributes = _parse_gtf_attributes_multi(attributes_text)
            gene_key = _first(attributes, "gene_id", "gene", "gene_name")
            gene_name = _first(attributes, "gene_name", "gene")
            identifiers = _gene_identifiers(
                attributes,
                include_cross_references=feature == "gene",
            )
            if gene_key is not None:
                metadata = gene_metadata.setdefault(
                    gene_key,
                    {"gene_name": None, "stable_ids": set()},
                )
                if gene_name is not None:
                    metadata["gene_name"] = metadata["gene_name"] or gene_name
                stable_ids = metadata["stable_ids"]
                assert isinstance(stable_ids, set)
                stable_ids.update(identifiers)

            transcript_id = _first(attributes, "transcript_id")
            if transcript_id is None:
                continue

            record = records.setdefault(
                transcript_id,
                {
                    "contig": None,
                    "strand": None,
                    "gene_name": None,
                    "gene_id": None,
                    "stable_ids": set(),
                    "is_mane_select": False,
                    "exon_coords": set(),
                },
            )

            contig = record["contig"]
            strand = record["strand"]
            if contig is None:
                record["contig"] = seqname
            elif contig != seqname:
                raise ValueError(f"{gtf_path}: transcript {transcript_id} spans multiple contigs")

            if strand is None:
                record["strand"] = feature_strand
            elif strand != feature_strand:
                raise ValueError(f"{gtf_path}: transcript {transcript_id} has inconsistent strands")
            if feature_strand not in {"+", "-"}:
                raise ValueError(
                    f"{gtf_path}:{line_number}: transcript {transcript_id} "
                    f"has unsupported strand {feature_strand!r}"
                )

            if record["gene_name"] is None:
                record["gene_name"] = gene_name or gene_key
            if record["gene_id"] is None:
                record["gene_id"] = gene_key
            stable_ids = record["stable_ids"]
            assert isinstance(stable_ids, set)
            stable_ids.update(identifiers)
            record["is_mane_select"] = bool(record["is_mane_select"]) or _is_mane_select(attributes)

            if feature == "exon":
                exon_coords = record["exon_coords"]
                assert isinstance(exon_coords, set)
                try:
                    start = int(start_text)
                    end = int(end_text)
                except ValueError as error:
                    raise ValueError(
                        f"{gtf_path}:{line_number}: transcript {transcript_id} "
                        "has non-integer exon coordinates"
                    ) from error
                if start > end:
                    raise ValueError(
                        f"{gtf_path}:{line_number}: transcript {transcript_id} "
                        "has exon start after end"
                    )
                exon_coords.add((start, end))

    annotations: dict[str, TranscriptAnnotation] = {}
    for transcript_id, record in records.items():
        exon_coords = record["exon_coords"]
        assert isinstance(exon_coords, set)
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
        metadata = gene_metadata.get(str(gene_id), {})
        transcript_gene_name = str(
            gene_name or metadata.get("gene_name") or gene_id or transcript_id
        )
        stable_ids = record["stable_ids"]
        assert isinstance(stable_ids, set)
        inherited_ids = metadata.get("stable_ids", set())
        assert isinstance(inherited_ids, set)
        stable_ids.update(inherited_ids)
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
            stable_ids=tuple(sorted(str(item) for item in stable_ids)),
            is_mane_select=bool(record["is_mane_select"]),
        )

    return annotations


def load_transcript_annotation(gtf_path: Path, transcript_id: str) -> TranscriptAnnotation:
    """Load one transcript annotation from a GTF."""
    annotations = load_transcript_annotations(gtf_path)
    return _resolve_transcript_annotation(annotations, transcript_id, gtf_path)


def _gtf_header_metadata(gtf_path: Path) -> dict[str, str]:
    metadata: dict[str, str] = {}
    with open_text(gtf_path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if not line.startswith("#!"):
                continue
            item = line[2:].strip()
            parts = item.split(None, 1)
            if len(parts) == 2:
                metadata[parts[0]] = parts[1]
    return metadata


@contextmanager
def _atomic_database_path(output: Path) -> Iterator[Path]:
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.",
        suffix=".tmp",
        dir=output.parent,
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    temporary_path.unlink()
    try:
        yield temporary_path
        os.replace(temporary_path, output)
    finally:
        temporary_path.unlink(missing_ok=True)


def index_gtf(gtf_path: Path, output: Path) -> tuple[int, int]:
    """Index transcript annotations from a GTF into a SQLite database."""
    annotations = load_transcript_annotations(gtf_path)
    transcript_count = len(annotations)
    exon_count = sum(len(annotation.exons) for annotation in annotations.values())

    with _atomic_database_path(output) as database_path, sqlite3.connect(
        database_path
    ) as connection:
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
                is_mane_select INTEGER NOT NULL,
                contig TEXT NOT NULL,
                strand TEXT NOT NULL,
                start INTEGER NOT NULL,
                end INTEGER NOT NULL
            );
            CREATE TABLE transcript_gene_ids (
                transcript_id TEXT NOT NULL,
                stable_id TEXT NOT NULL COLLATE NOCASE,
                PRIMARY KEY (transcript_id, stable_id),
                FOREIGN KEY (transcript_id) REFERENCES transcripts(transcript_id) ON DELETE CASCADE
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
            CREATE INDEX idx_transcripts_gene ON transcripts (gene_name COLLATE NOCASE);
            CREATE INDEX idx_transcript_gene_ids_id ON transcript_gene_ids (stable_id COLLATE NOCASE);
            CREATE INDEX idx_exons_transcript ON exons (transcript_id);
            """
        )
        metadata = {
            "source_gtf": str(gtf_path),
            "schema_version": "2",
            **_gtf_header_metadata(gtf_path),
        }
        connection.executemany(
            "INSERT INTO metadata(key, value) VALUES(?, ?)",
            sorted(metadata.items()),
        )
        connection.executemany(
            """
            INSERT INTO transcripts(
                transcript_id,
                transcript_id_root,
                gene_name,
                is_mane_select,
                contig,
                strand,
                start,
                end
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    annotation.transcript_id,
                    annotation.transcript_id.split(".", 1)[0],
                    annotation.gene_name,
                    int(annotation.is_mane_select),
                    annotation.contig,
                    annotation.strand,
                    annotation.start,
                    annotation.end,
                )
                for annotation in annotations.values()
            ],
        )
        connection.executemany(
            "INSERT INTO transcript_gene_ids(transcript_id, stable_id) VALUES (?, ?)",
            [
                (annotation.transcript_id, stable_id)
                for annotation in annotations.values()
                for stable_id in annotation.stable_ids
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
        if not self.path.exists():
            raise ValueError(f"{self.path}: transcript database not found")
        if not self.path.is_file():
            raise ValueError(f"{self.path}: transcript database is not a file")

        self._connection: sqlite3.Connection | None = None
        try:
            database_uri = f"file:{quote(str(self.path.resolve()), safe='/')}?mode=ro"
            self._connection = sqlite3.connect(database_uri, uri=True)
            self._connection.row_factory = sqlite3.Row
            self._validate_schema()
        except ValueError:
            self.close()
            raise
        except sqlite3.Error as error:
            self.close()
            raise ValueError(f"{self.path}: invalid transcript database: {error}") from error

    def _validate_schema(self) -> None:
        assert self._connection is not None
        required_columns = {
            "metadata": {"key", "value"},
            "transcripts": {
                "transcript_id",
                "transcript_id_root",
                "gene_name",
                "is_mane_select",
                "contig",
                "strand",
                "start",
                "end",
            },
            "transcript_gene_ids": {"transcript_id", "stable_id"},
            "exons": {"transcript_id", "exon_number", "start", "end"},
        }
        rows = self._connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        ).fetchall()
        tables = {str(row["name"]) for row in rows}
        missing = set(required_columns) - tables
        if missing:
            missing_text = ", ".join(sorted(missing))
            raise ValueError(
                f"{self.path}: invalid transcript database; "
                f"missing table(s): {missing_text}"
            )
        for table, columns in required_columns.items():
            column_rows = self._connection.execute(
                f"PRAGMA table_info({table})"
            ).fetchall()
            existing_columns = {str(row["name"]) for row in column_rows}
            missing_columns = columns - existing_columns
            if missing_columns:
                missing_text = ", ".join(sorted(missing_columns))
                raise ValueError(
                    f"{self.path}: invalid transcript database; "
                    f"table {table} missing column(s): {missing_text}"
                )

    def close(self) -> None:
        """Close the underlying SQLite connection."""
        if self._connection is not None:
            self._connection.close()
            self._connection = None

    def metadata(self) -> dict[str, str]:
        """Return metadata stored when the transcript index was built."""
        assert self._connection is not None
        return {
            str(row["key"]): str(row["value"])
            for row in self._connection.execute(
                "SELECT key, value FROM metadata ORDER BY key"
            )
        }

    def __enter__(self) -> "TranscriptIndex":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    def get(self, transcript_id: str) -> TranscriptAnnotation:
        """Return one transcript annotation by exact ID or unversioned ID when unambiguous."""
        assert self._connection is not None
        transcript_row = self._connection.execute(
            """
            SELECT transcript_id, gene_name, is_mane_select, contig, strand, start, end
            FROM transcripts
            WHERE transcript_id = ?
            """,
            (transcript_id,),
        ).fetchone()
        if transcript_row is None:
            root_id = transcript_id.split(".", 1)[0]
            transcript_rows = self._connection.execute(
                """
                SELECT transcript_id, gene_name, is_mane_select, contig, strand, start, end
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
        identifier_rows = self._connection.execute(
            """
            SELECT stable_id
            FROM transcript_gene_ids
            WHERE transcript_id = ?
            ORDER BY stable_id COLLATE NOCASE
            """,
            (transcript_row["transcript_id"],),
        ).fetchall()
        return TranscriptAnnotation(
            transcript_id=str(transcript_row["transcript_id"]),
            gene_name=str(transcript_row["gene_name"]),
            contig=str(transcript_row["contig"]),
            strand=str(transcript_row["strand"]),
            start=int(transcript_row["start"]),
            end=int(transcript_row["end"]),
            exons=exons,
            stable_ids=tuple(str(row["stable_id"]) for row in identifier_rows),
            is_mane_select=bool(transcript_row["is_mane_select"]),
        )

    def resolve_gene(
        self,
        gene_name: str,
        stable_id: str,
        transcript_id: str | None = None,
        *,
        available_contigs: set[str] | None = None,
    ) -> TranscriptAnnotation:
        """Resolve a transcript after cross-checking its gene name and stable ID."""
        assert self._connection is not None
        parameters: list[object] = [gene_name, normalize_gene_identifier(stable_id)]
        where = (
            "t.gene_name = ? COLLATE NOCASE "
            "AND EXISTS ("
            "SELECT 1 FROM transcript_gene_ids AS i "
            "WHERE i.transcript_id = t.transcript_id "
            "AND i.stable_id = ? COLLATE NOCASE)"
        )
        if transcript_id is None:
            where += " AND t.is_mane_select = 1"
        else:
            where += " AND (t.transcript_id = ? OR t.transcript_id_root = ?)"
            parameters.extend((transcript_id, transcript_id.split(".", 1)[0]))
        rows = self._connection.execute(
            f"SELECT t.transcript_id FROM transcripts AS t WHERE {where} "
            "ORDER BY t.transcript_id",
            parameters,
        ).fetchall()
        if not rows:
            requested = (
                f"transcript {transcript_id}"
                if transcript_id is not None
                else "a MANE Select transcript"
            )
            raise ValueError(
                f"{gene_name} ({stable_id}): could not find {requested}"
            )
        if len(rows) > 1:
            matches = ", ".join(str(row["transcript_id"]) for row in rows)
            raise ValueError(
                f"{gene_name} ({stable_id}): ambiguous transcripts: {matches}"
            )
        annotation = self.get(str(rows[0]["transcript_id"]))
        if available_contigs is None:
            return annotation
        return adapt_transcript_contig(
            annotation,
            available_contigs,
            label=gene_name,
        )


def load_transcript_annotation_from_db(db_path: Path, transcript_id: str) -> TranscriptAnnotation:
    """Load one transcript annotation from a SQLite transcript database."""
    with TranscriptIndex(db_path) as transcript_index:
        return transcript_index.get(transcript_id)
