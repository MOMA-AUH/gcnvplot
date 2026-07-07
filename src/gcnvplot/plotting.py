"""Plot data preparation and public plotting API."""

from __future__ import annotations

from pathlib import Path

from .background import load_background, log2_ratio, normalized_counts
from .models import BackgroundSummary, Interval, TranscriptAnnotation
from .read_counts import parse_read_counts
from .rendering import _build_plot_svg
from .transcripts import TranscriptIndex
from .utils import overlaps, parse_region


def plot_rows(
    counts: dict[Interval, int],
    region: Interval,
    background_summary: BackgroundSummary,
    transcript: TranscriptAnnotation | None = None,
) -> tuple[list[dict[str, object]], int]:
    """Join sample counts with the background for plotting."""
    background = background_summary.rows
    baselines = {
        interval: float(row["BASELINE_MEDIAN"])
        for interval, row in background.items()
        if float(row["BASELINE_MEDIAN"]) > 0
    }
    sample_norm = normalized_counts(counts, baselines)
    rows: list[dict[str, object]] = []
    missing_background = 0
    for interval in sorted(counts, key=lambda item: (item[0], item[1], item[2])):
        if interval[0] != region[0] or not overlaps(interval[1], interval[2], region[1], region[2]):
            continue
        baseline = baselines.get(interval)
        if baseline is None or baseline <= 0:
            missing_background += 1
            continue
        bg = background.get(interval)
        if bg is None:
            missing_background += 1
            continue
        sample_value = sample_norm[interval]
        bg_median = float(bg["BG_NORM_MEDIAN"])
        bg_lower = float(bg["BG_NORM_P5"])
        bg_upper = float(bg["BG_NORM_P95"])
        row: dict[str, object] = {
            "contig": interval[0],
            "start": interval[1],
            "end": interval[2],
            "mid": (interval[1] + interval[2]) / 2,
            "count": counts[interval],
            "sample_norm": sample_value,
            "bg_median_norm": bg_median,
            "bg_lower_norm": bg_lower,
            "bg_upper_norm": bg_upper,
            "background_n": int(bg["N"]),
        }
        row["overlaps_exon"] = transcript is None or any(
            overlaps(interval[1], interval[2], exon.start, exon.end) for exon in transcript.exons
        )
        row["signal"] = log2_ratio(sample_value, bg_median, 0.01)
        row["signal_lower"] = log2_ratio(bg_lower, bg_median, 0.01)
        row["signal_upper"] = log2_ratio(bg_upper, bg_median, 0.01)
        rows.append(row)

    if missing_background:
        return rows, missing_background
    return rows, 0


def _coerce_region(value: Interval | str) -> Interval:
    if isinstance(value, str):
        return parse_region(value)
    return value


def _coerce_highlight(value: Interval | str | None) -> Interval | None:
    if value is None:
        return None
    return _coerce_region(value)


def render_plot_svg(
    read_counts: str | Path | dict[Interval, int],
    background: str | Path | BackgroundSummary,
    *,
    region: Interval | str | None = None,
    transcript_id: str | None = None,
    transcript_index: str | Path | TranscriptIndex | None = None,
    sample_name: str | None = None,
    highlight: Interval | str | None = None,
    strict_background: bool = True,
) -> str:
    """Render a gcnvplot SVG for use in Python code or reports.

    Provide exactly one of ``region`` or ``transcript_id``. When
    ``strict_background`` is true, selected intervals without usable background
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

    transcript: TranscriptAnnotation | None = None
    if transcript_id is not None:
        if isinstance(transcript_index, TranscriptIndex):
            transcript = transcript_index.get(transcript_id)
        else:
            with TranscriptIndex(Path(transcript_index)) as index:
                transcript = index.get(transcript_id)
        region_value: Interval = (transcript.contig, transcript.start, transcript.end)
    else:
        assert region is not None
        region_value = _coerce_region(region)

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
    rows, missing_background = plot_rows(
        counts,
        region_value,
        background_summary,
        transcript=transcript,
    )
    if missing_background and strict_background:
        interval_text = "interval is" if missing_background == 1 else "intervals are"
        raise ValueError(
            f"{missing_background} selected {interval_text} missing usable background statistics"
        )
    if not rows:
        raise ValueError("No intervals with background statistics overlap the selected region.")

    return _build_plot_svg(
        rows,
        region_value,
        sample_name=sample_name,
        highlight=_coerce_highlight(highlight),
        transcript=transcript,
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
    strict_background: bool = True,
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
        strict_background=strict_background,
    )
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
