"""Shared data models for gcnvplot."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

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
    stable_ids: tuple[str, ...] = ()
    is_mane_select: bool = False


@dataclass(frozen=True)
class BackgroundSummary:
    """Background intervals plus metadata used in plots."""

    rows: dict[Interval, dict[str, str]]
    lower_percentile: int | None = None
    upper_percentile: int | None = None
    metadata: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class PlotOverlay:
    """A typed genomic band rendered on a plot."""

    interval: Interval
    kind: Literal["highlight", "deletion", "duplication"]
    label: str | None = None
