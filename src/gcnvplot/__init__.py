"""Top-level package for gcnvplot."""

from .background import create_background, load_background, write_background
from .models import BackgroundSummary, Interval, PlotOverlay, TranscriptAnnotation, TranscriptExon
from .plotting import (
    PreparedPlotData,
    prepare_plot_data,
    render_plot_svg,
    render_prepared_plot_svg,
    write_plot,
)
from .read_counts import parse_read_counts
from .report import GeneRequest, ReportBuildResult, build_report, read_gene_requests
from .segments import CnvSegment, parse_segments_vcf
from .transcripts import TranscriptIndex, index_gtf
from .utils import parse_region

__version__ = "1.1.0"

__all__ = [
    "BackgroundSummary",
    "Interval",
    "PlotOverlay",
    "PreparedPlotData",
    "GeneRequest",
    "ReportBuildResult",
    "CnvSegment",
    "TranscriptAnnotation",
    "TranscriptExon",
    "TranscriptIndex",
    "__version__",
    "create_background",
    "build_report",
    "index_gtf",
    "load_background",
    "parse_read_counts",
    "parse_segments_vcf",
    "prepare_plot_data",
    "parse_region",
    "read_gene_requests",
    "render_plot_svg",
    "render_prepared_plot_svg",
    "write_background",
    "write_plot",
]
