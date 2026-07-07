"""Top-level package for gcnvplot."""

from .background import load_background
from .models import BackgroundSummary, Interval, TranscriptAnnotation, TranscriptExon
from .plotting import render_plot_svg, write_plot
from .read_counts import parse_read_counts
from .transcripts import TranscriptIndex, index_gtf
from .utils import parse_region

__version__ = "1.0.0"

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
