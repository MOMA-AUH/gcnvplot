"""Plot data preparation and public plotting API."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from .background import load_background, log2_ratio, normalized_counts
from .models import BackgroundSummary, Interval, PlotOverlay, TranscriptAnnotation
from .read_counts import parse_read_counts
from .rendering import _build_plot_svg
from .transcripts import TranscriptIndex, adapt_transcript_contig
from .utils import overlaps, parse_region


@dataclass(frozen=True)
class PreparedPlotData:
    """Parsed and normalized sample data reusable across many plot regions."""

    counts: dict[Interval, int]
    background: BackgroundSummary
    sample_norm: dict[Interval, float]
    counts_by_contig: dict[str, tuple[Interval, ...]]
    background_by_contig: dict[str, frozenset[Interval]]


def prepare_plot_data(
    read_counts: str | Path | dict[Interval, int],
    background: str | Path | BackgroundSummary,
) -> PreparedPlotData:
    """Parse, normalize, and index one sample for repeated plotting."""
    background_summary = (
        background
        if isinstance(background, BackgroundSummary)
        else load_background(Path(background))
    )
    counts = (
        read_counts
        if isinstance(read_counts, dict)
        else parse_read_counts(Path(read_counts))
    )
    baselines = {
        interval: float(row["BASELINE_MEDIAN"])
        for interval, row in background_summary.rows.items()
        if float(row["BASELINE_MEDIAN"]) > 0
    }
    sample_norm = normalized_counts(counts, baselines)
    count_index: dict[str, list[Interval]] = defaultdict(list)
    for interval in counts:
        count_index[interval[0]].append(interval)
    background_index: dict[str, set[Interval]] = defaultdict(set)
    for interval in background_summary.rows:
        background_index[interval[0]].add(interval)
    return PreparedPlotData(
        counts=counts,
        background=background_summary,
        sample_norm=sample_norm,
        counts_by_contig={
            contig: tuple(sorted(intervals, key=lambda interval: (interval[1], interval[2])))
            for contig, intervals in count_index.items()
        },
        background_by_contig={
            contig: frozenset(intervals)
            for contig, intervals in background_index.items()
        },
    )


def _prepared_plot_rows(
    prepared: PreparedPlotData,
    region: Interval,
    transcript: TranscriptAnnotation | None = None,
) -> tuple[list[dict[str, object]], int]:
    """Build plot rows from data normalized once for the complete report."""
    rows: list[dict[str, object]] = []
    missing_background = 0
    for interval in prepared.counts_by_contig.get(region[0], ()):
        if not overlaps(interval[1], interval[2], region[1], region[2]):
            continue
        bg = prepared.background.rows.get(interval)
        sample_value = prepared.sample_norm.get(interval)
        if (
            bg is None
            or sample_value is None
            or float(bg["BASELINE_MEDIAN"]) <= 0
        ):
            missing_background += 1
            continue
        bg_median = float(bg["BG_NORM_MEDIAN"])
        bg_lower = float(bg["BG_NORM_P5"])
        bg_upper = float(bg["BG_NORM_P95"])
        row: dict[str, object] = {
            "contig": interval[0],
            "start": interval[1],
            "end": interval[2],
            "mid": (interval[1] + interval[2]) / 2,
            "count": prepared.counts[interval],
            "sample_norm": sample_value,
            "bg_median_norm": bg_median,
            "bg_lower_norm": bg_lower,
            "bg_upper_norm": bg_upper,
            "background_n": int(bg["N"]),
        }
        row["overlaps_exon"] = transcript is None or any(
            overlaps(interval[1], interval[2], exon.start, exon.end)
            for exon in transcript.exons
        )
        row["signal"] = log2_ratio(sample_value, bg_median, 0.01)
        row["signal_lower"] = log2_ratio(bg_lower, bg_median, 0.01)
        row["signal_upper"] = log2_ratio(bg_upper, bg_median, 0.01)
        rows.append(row)
    return rows, missing_background


def plot_rows(
    counts: dict[Interval, int],
    region: Interval,
    background_summary: BackgroundSummary,
    transcript: TranscriptAnnotation | None = None,
) -> tuple[list[dict[str, object]], int]:
    """Join sample counts with the background for plotting."""
    prepared = prepare_plot_data(counts, background_summary)
    return _prepared_plot_rows(prepared, region, transcript)


def _coerce_region(value: Interval | str) -> Interval:
    if isinstance(value, str):
        return parse_region(value)
    return value


def _coerce_highlight(value: Interval | str | None) -> Interval | None:
    if value is None:
        return None
    return _coerce_region(value)


def _selected_intervals(intervals: Iterable[Interval], region: Interval) -> set[Interval]:
    return {
        interval
        for interval in intervals
        if interval[0] == region[0]
        and overlaps(interval[1], interval[2], region[1], region[2])
    }


def _format_interval(interval: Interval) -> str:
    contig, start, end = interval
    return f"{contig}:{start}-{end}"


def _format_interval_examples(intervals: set[Interval], limit: int = 3) -> str:
    examples = sorted(intervals, key=lambda item: (item[0], item[1], item[2]))[:limit]
    return ", ".join(_format_interval(interval) for interval in examples)


def interval_set_mismatches_for_region(
    counts: dict[Interval, int],
    region: Interval,
    background_summary: BackgroundSummary,
) -> tuple[set[Interval], set[Interval]]:
    """Return selected sample-only and background-only intervals."""
    sample_intervals = _selected_intervals(counts, region)
    background_intervals = _selected_intervals(background_summary.rows, region)
    return sample_intervals - background_intervals, background_intervals - sample_intervals


def format_interval_set_mismatch(
    sample_only: set[Interval],
    background_only: set[Interval],
) -> str:
    """Return a concise message for selected interval-set mismatches."""
    details = []
    if sample_only:
        details.append(
            f"{len(sample_only)} sample interval(s) missing from the background"
            f" ({_format_interval_examples(sample_only)})"
        )
    if background_only:
        details.append(
            f"{len(background_only)} background interval(s) missing from the sample"
            f" ({_format_interval_examples(background_only)})"
        )
    return "Interval sets do not match for the selected region: " + "; ".join(details)


def _prepared_interval_mismatches(
    prepared: PreparedPlotData,
    region: Interval,
) -> tuple[set[Interval], set[Interval]]:
    sample_intervals = {
        interval
        for interval in prepared.counts_by_contig.get(region[0], ())
        if overlaps(interval[1], interval[2], region[1], region[2])
    }
    background_intervals = {
        interval
        for interval in prepared.background_by_contig.get(region[0], frozenset())
        if overlaps(interval[1], interval[2], region[1], region[2])
    }
    return sample_intervals - background_intervals, background_intervals - sample_intervals


def render_prepared_plot_svg(
    prepared: PreparedPlotData,
    region: Interval,
    *,
    transcript: TranscriptAnnotation | None = None,
    sample_name: str | None = None,
    highlight: Interval | str | None = None,
    overlays: list[PlotOverlay] | None = None,
    strict_intervals: bool = True,
) -> str:
    """Render one region from data prepared once for repeated plotting."""
    sample_only, background_only = _prepared_interval_mismatches(prepared, region)
    if strict_intervals and (sample_only or background_only):
        raise ValueError(format_interval_set_mismatch(sample_only, background_only))
    rows, missing_background = _prepared_plot_rows(prepared, region, transcript)
    if missing_background and strict_intervals:
        interval_text = "interval is" if missing_background == 1 else "intervals are"
        raise ValueError(
            f"{missing_background} selected {interval_text} missing usable background statistics"
        )
    if not rows:
        raise ValueError("No intervals with background statistics overlap the selected region.")
    return _build_plot_svg(
        rows,
        region,
        sample_name=sample_name,
        highlight=_coerce_highlight(highlight),
        transcript=transcript,
        overlays=overlays,
    )


def render_plot_svg(
    read_counts: str | Path | dict[Interval, int],
    background: str | Path | BackgroundSummary,
    *,
    region: Interval | str | None = None,
    transcript_id: str | None = None,
    transcript_index: str | Path | TranscriptIndex | None = None,
    sample_name: str | None = None,
    highlight: Interval | str | None = None,
    overlays: list[PlotOverlay] | None = None,
    strict_intervals: bool = True,
) -> str:
    """Render a gcnvplot SVG for use in Python code or reports.

    Provide exactly one of ``region`` or ``transcript_id``. When
    ``strict_intervals`` is true, the selected sample/background interval sets
    must match exactly and selected intervals without usable background
    statistics raise ``ValueError`` instead of being silently skipped.
    """
    if region is not None and transcript_id is not None:
        raise ValueError("region and transcript_id are mutually exclusive")
    if transcript_id is not None and transcript_index is None:
        raise ValueError("transcript_index is required when transcript_id is provided")
    if transcript_id is None and transcript_index is not None:
        raise ValueError("transcript_index can only be used when transcript_id is provided")
    if transcript_id is None and region is None:
        raise ValueError("region is required when transcript_id is not provided")

    prepared = prepare_plot_data(read_counts, background)
    transcript: TranscriptAnnotation | None = None
    if transcript_id is not None:
        if isinstance(transcript_index, TranscriptIndex):
            transcript = transcript_index.get(transcript_id)
        else:
            with TranscriptIndex(Path(transcript_index)) as index:
                transcript = index.get(transcript_id)
        shared_contigs = set(prepared.counts_by_contig) & set(
            prepared.background_by_contig
        )
        transcript = adapt_transcript_contig(transcript, shared_contigs)
        region_value: Interval = (
            transcript.contig,
            transcript.start,
            transcript.end,
        )
    else:
        assert region is not None
        region_value = _coerce_region(region)

    return render_prepared_plot_svg(
        prepared,
        region_value,
        transcript=transcript,
        sample_name=sample_name,
        highlight=highlight,
        overlays=overlays,
        strict_intervals=strict_intervals,
    )


def write_plot(
    read_counts: str | Path | dict[Interval, int],
    background: str | Path | BackgroundSummary,
    output: str | Path,
    *,
    region: Interval | str | None = None,
    transcript_id: str | None = None,
    transcript_index: str | Path | TranscriptIndex | None = None,
    sample_name: str | None = None,
    highlight: Interval | str | None = None,
    overlays: list[PlotOverlay] | None = None,
    strict_intervals: bool = True,
) -> None:
    """Render and write a gcnvplot SVG."""
    svg = render_plot_svg(
        read_counts,
        background,
        region=region,
        transcript_id=transcript_id,
        transcript_index=transcript_index,
        sample_name=sample_name,
        highlight=highlight,
        overlays=overlays,
        strict_intervals=strict_intervals,
    )
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
