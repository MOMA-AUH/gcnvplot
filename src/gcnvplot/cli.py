"""Command-line interface for gcnvplot."""

from __future__ import annotations

import argparse
from pathlib import Path

from . import __version__
from .core import build_background, parse_region, plot_sample


def build_parser() -> argparse.ArgumentParser:
    """Build and return the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="gcnvplot",
        description="Plot GATK germline CNV read-count signals against a background cohort.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser(
        "create-background",
        help="Create interval-wise normalized background statistics.",
    )
    build.add_argument(
        "--read-counts-list",
        type=Path,
        required=True,
        help="Text file with one read_counts.tsv path per line.",
    )
    build.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output TSV path for the background cohort summary.",
    )
    build.set_defaults(handler=build_background)

    plot = subparsers.add_parser(
        "plot",
        help="Plot one sample against a background cohort.",
    )
    plot.add_argument("--read-counts", type=Path, required=True, help="GATK read_counts.tsv for the sample to plot.")
    plot.add_argument("--background", type=Path, required=True, help="Background TSV from create-background.")
    region_or_transcript = plot.add_mutually_exclusive_group(required=True)
    region_or_transcript.add_argument(
        "--region",
        type=parse_region,
        help="Genomic region to plot, e.g. chr1:100-200.",
    )
    region_or_transcript.add_argument(
        "--transcript",
        help="Transcript ID to plot using exon coordinates from a GTF file.",
    )
    plot.add_argument(
        "--gtf",
        type=Path,
        help="GTF file for --transcript; gzipped files are supported.",
    )
    plot.add_argument(
        "--highlight",
        type=parse_region,
        help="Optional genomic region to highlight as a semi-transparent band.",
    )
    plot.add_argument(
        "--sample-name",
        help="Optional sample label to show in the plot info panel.",
    )
    plot.add_argument("--output", type=Path, required=True, help="SVG output path.")
    plot.set_defaults(handler=plot_sample)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "transcript", None) is not None and getattr(args, "gtf", None) is None:
        parser.error("--gtf is required when --transcript is used")
    handler = getattr(args, "handler", None)
    if handler is None:
        parser.print_help()
        return 1
    return handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
