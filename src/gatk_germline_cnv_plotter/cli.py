"""Command-line interface for gatk-germline-cnv-plotter."""

from __future__ import annotations

import argparse
from pathlib import Path

from . import __version__


def build_parser() -> argparse.ArgumentParser:
    """Build and return the top-level argument parser."""
    parser = argparse.ArgumentParser(
        prog="gatk-germline-cnv-plotter",
        description="Plot GATK germline CNV outputs.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    plot_parser = subparsers.add_parser(
        "plot",
        help="Generate a CNV plot from input files.",
        description="Generate a CNV plot from GATK germline CNV outputs.",
    )
    plot_parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to the primary input file or directory.",
    )
    plot_parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Path to the output plot file.",
    )
    plot_parser.set_defaults(handler=_run_plot_command)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the CLI."""
    parser = build_parser()
    args = parser.parse_args(argv)
    handler = getattr(args, "handler", None)
    if handler is None:
        parser.print_help()
        return 1
    return handler(args)


def _run_plot_command(args: argparse.Namespace) -> int:
    """Run the placeholder plot command."""
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print("Plot generation is not implemented yet.")
    return 0
