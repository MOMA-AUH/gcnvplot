"""Compatibility facade for gcnvplot helpers.

The implementation is split across focused modules. Imports from
``gcnvplot.core`` are kept working for callers that used earlier releases.
"""

from __future__ import annotations

from .background import (
    BACKGROUND_FIELDS,
    interval_baselines,
    load_background,
    log2_ratio,
    normalized_counts,
    parse_background,
    size_factor_from_baseline,
    write_background,
)
from .cli_handlers import build_background, index_transcripts, plot_sample, rows_for_plot
from .models import BackgroundSummary, Interval, TranscriptAnnotation, TranscriptExon
from .plotting import _coerce_highlight, _coerce_region, plot_rows, render_plot_svg, write_plot
from .read_counts import READ_COUNT_FIELDS, parse_read_counts, read_path_list
from .rendering import _build_plot_svg, nice_ticks, write_svg_plot
from .transcripts import (
    TranscriptIndex,
    _resolve_transcript_annotation,
    index_gtf,
    load_transcript_annotation,
    load_transcript_annotation_from_db,
    load_transcript_annotations,
    parse_gtf_attributes,
    transcript_id_matches,
)
from .utils import fmt, median, open_text, overlaps, parse_region, percentile, stdev

__all__ = [
    "BACKGROUND_FIELDS",
    "BackgroundSummary",
    "Interval",
    "READ_COUNT_FIELDS",
    "TranscriptAnnotation",
    "TranscriptExon",
    "TranscriptIndex",
    "_build_plot_svg",
    "_coerce_highlight",
    "_coerce_region",
    "_resolve_transcript_annotation",
    "build_background",
    "fmt",
    "index_gtf",
    "index_transcripts",
    "interval_baselines",
    "load_background",
    "load_transcript_annotation",
    "load_transcript_annotation_from_db",
    "load_transcript_annotations",
    "log2_ratio",
    "median",
    "nice_ticks",
    "normalized_counts",
    "open_text",
    "overlaps",
    "parse_background",
    "parse_gtf_attributes",
    "parse_read_counts",
    "parse_region",
    "percentile",
    "plot_rows",
    "plot_sample",
    "read_path_list",
    "render_plot_svg",
    "rows_for_plot",
    "size_factor_from_baseline",
    "stdev",
    "transcript_id_matches",
    "write_background",
    "write_plot",
    "write_svg_plot",
]
