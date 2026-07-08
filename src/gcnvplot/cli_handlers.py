"""Argparse command handlers for gcnvplot."""

from __future__ import annotations

import argparse
import sys

from .background import parse_background, write_background
from .models import BackgroundSummary, Interval, TranscriptAnnotation
from .plotting import (
    format_interval_set_mismatch,
    interval_set_mismatches_for_region,
    plot_rows,
)
from .read_counts import parse_read_counts, read_path_list
from .rendering import write_svg_plot
from .transcripts import index_gtf, load_transcript_annotation_from_db


def index_transcripts(args: argparse.Namespace) -> int:
    """CLI handler for indexing transcript annotations into SQLite."""
    transcript_count, exon_count = index_gtf(args.gtf, args.output)
    print(f"Indexed transcripts: {transcript_count}")
    print(f"Indexed exons: {exon_count}")
    print(f"Wrote: {args.output}")
    return 0


def build_background(args: argparse.Namespace) -> int:
    """Build interval-wise background summaries."""
    files = read_path_list(args.read_counts_list)
    if not files:
        raise SystemExit(f"{args.read_counts_list}: no read-count TSV paths found")

    if args.allow_interval_mismatches:
        print(
            "Warning: interval-set validation is disabled; "
            "mismatched background intervals can affect normalization.",
            file=sys.stderr,
        )

    written_intervals = write_background(
        files,
        args.output,
        strict_intervals=not args.allow_interval_mismatches,
    )

    print(f"Background samples: {len(files)}")
    print(f"Background intervals: {written_intervals}")
    print(f"Wrote: {args.output}")
    return 0


def rows_for_plot(
    args: argparse.Namespace,
    region: Interval,
    background_summary: BackgroundSummary,
    transcript: TranscriptAnnotation | None = None,
) -> list[dict[str, object]]:
    """CLI helper to join one sample with the background for plotting."""
    counts = parse_read_counts(args.read_counts)
    strict_intervals = not getattr(args, "allow_interval_mismatches", False)
    sample_only, background_only = interval_set_mismatches_for_region(
        counts,
        region,
        background_summary,
    )
    if sample_only or background_only:
        message = format_interval_set_mismatch(sample_only, background_only)
        if strict_intervals:
            raise ValueError(message)
        print(f"Warning: {message}", file=sys.stderr)

    rows, missing_background = plot_rows(counts, region, background_summary, transcript=transcript)
    if missing_background:
        if strict_intervals:
            interval_text = "interval is" if missing_background == 1 else "intervals are"
            raise ValueError(
                f"{missing_background} selected {interval_text} missing usable background statistics"
            )
        print(
            f"Warning: skipped intervals missing usable background statistics: {missing_background}",
            file=sys.stderr,
        )
    return rows


def plot_sample(args: argparse.Namespace) -> int:
    """Plot one sample against the background."""
    transcript: TranscriptAnnotation | None = None
    region = args.region
    if args.transcript is not None:
        transcript = load_transcript_annotation_from_db(args.transcript_db, args.transcript)
        region = (transcript.contig, transcript.start, transcript.end)

    background_summary = parse_background(args.background)
    rows = rows_for_plot(args, region, background_summary, transcript=transcript)
    if not rows:
        raise SystemExit("No intervals with background statistics overlap the selected region.")

    write_svg_plot(
        rows,
        args.output,
        region,
        sample_name=args.sample_name,
        highlight=args.highlight,
        transcript=transcript,
    )
    print(f"Plotted intervals: {len(rows)}")
    print(f"Wrote: {args.output}")
    return 0
