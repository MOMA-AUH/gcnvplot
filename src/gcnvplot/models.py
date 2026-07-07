"""Shared data models for gcnvplot."""

from __future__ import annotations

from dataclasses import dataclass

Interval = tuple[str, int, int]


@dataclass(frozen=True)
class TranscriptExon:
    """One exon in transcript order."""

    number: int
    start: int
    end: int


@dataclass(frozen=True)
class TranscriptAnnotation:
    """Transcript model used for plotting."""

    transcript_id: str
    gene_name: str
    contig: str
    strand: str
    start: int
    end: int
    exons: list[TranscriptExon]


@dataclass(frozen=True)
class BackgroundSummary:
    """Background intervals plus metadata used in plots."""

    rows: dict[Interval, dict[str, str]]
    lower_percentile: int | None = None
    upper_percentile: int | None = None
