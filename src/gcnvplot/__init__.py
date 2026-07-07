"""Top-level package for gcnvplot."""

from .core import (
    BackgroundSummary,
    Interval,
    TranscriptAnnotation,
    TranscriptExon,
    TranscriptIndex,
    index_gtf,
    load_background,
    parse_read_counts,
    parse_region,
    render_plot_svg,
    write_plot,
)

__version__ = "0.5.0"

__all__ = [
    "BackgroundSummary",
    "Interval",
    "TranscriptAnnotation",
    "TranscriptExon",
    "TranscriptIndex",
    "__version__",
    "index_gtf",
    "load_background",
    "parse_read_counts",
    "parse_region",
    "render_plot_svg",
    "write_plot",
]
