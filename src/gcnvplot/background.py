"""Background cohort normalization and summary handling."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

from .models import BackgroundSummary, Interval
from .read_counts import parse_read_counts
from .utils import fmt, median, open_text, percentile, stdev

BACKGROUND_FIELDS = [
    "CONTIG",
    "START",
    "END",
    "BASELINE_MEDIAN",
    "N",
    "BG_NORM_MEAN",
    "BG_NORM_MEDIAN",
    "BG_NORM_SD",
    "BG_NORM_P5",
    "BG_NORM_P95",
]


def parse_background(path: Path) -> BackgroundSummary:
    """Parse a background TSV created by create-background."""
    lower_percentile: int | None = None
    upper_percentile: int | None = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                if line.startswith("# lower_percentile="):
                    lower_percentile = int(line.split("=", 1)[1])
                elif line.startswith("# upper_percentile="):
                    upper_percentile = int(line.split("=", 1)[1])
                continue
            header = line.rstrip("\n").split("\t")
            if header != BACKGROUND_FIELDS:
                raise ValueError(f"{path}: unexpected background header: {header}")
            break
        else:
            raise ValueError(f"{path}: no background header found")

        background: dict[Interval, dict[str, str]] = {}
        reader = csv.DictReader(handle, fieldnames=BACKGROUND_FIELDS, delimiter="\t")
        for row in reader:
            interval = (row["CONTIG"], int(row["START"]), int(row["END"]))
            background[interval] = row
    return BackgroundSummary(
        rows=background,
        lower_percentile=lower_percentile,
        upper_percentile=upper_percentile,
    )


def load_background(path: Path) -> BackgroundSummary:
    """Load a background summary TSV."""
    return parse_background(path)


def interval_baselines(sample_counts: list[dict[Interval, int]]) -> dict[Interval, float]:
    """Estimate a robust baseline for each interval from background samples."""
    interval_values: dict[Interval, list[float]] = defaultdict(list)
    for counts in sample_counts:
        for interval, count in counts.items():
            if count > 0:
                interval_values[interval].append(float(count))

    baselines: dict[Interval, float] = {}
    for interval, values in interval_values.items():
        baselines[interval] = median(values)
    return baselines


def size_factor_from_baseline(
    counts: dict[Interval, int], baselines: dict[Interval, float]
) -> float:
    """Return a DESeq-style median-of-ratios size factor."""
    ratios = [
        count / baseline
        for interval, count in counts.items()
        if count > 0 and (baseline := baselines.get(interval)) is not None and baseline > 0
    ]
    if not ratios:
        raise ValueError("Cannot normalize a file without usable baseline intervals")
    return median(ratios)


def normalized_counts(counts: dict[Interval, int], baselines: dict[Interval, float]) -> dict[Interval, float]:
    """Return interval counts divided by the median-of-ratios size factor."""
    factor = size_factor_from_baseline(counts, baselines)
    return {interval: count / factor for interval, count in counts.items()}


def log2_ratio(sample_value: float, expected_value: float, pseudocount: float) -> float:
    """Return log2(sample / expected) with a small stabilizing pseudocount."""
    if pseudocount <= 0:
        raise ValueError("pseudocount must be > 0")
    return math.log2((sample_value + pseudocount) / (expected_value + pseudocount))


def write_background(read_count_paths: list[Path], output: Path) -> int:
    """Build and write interval-wise background summaries."""
    sample_counts = [parse_read_counts(path) for path in read_count_paths]
    baselines = interval_baselines(sample_counts)

    interval_values: dict[Interval, list[float]] = defaultdict(list)
    for counts in sample_counts:
        for interval, value in normalized_counts(counts, baselines).items():
            interval_values[interval].append(value)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        handle.write("# normalization=median-of-ratios\n")
        handle.write("# baseline=median-positive-count\n")
        handle.write(f"# samples={len(read_count_paths)}\n")
        handle.write("# lower_percentile=5\n")
        handle.write("# upper_percentile=95\n")
        writer = csv.DictWriter(handle, fieldnames=BACKGROUND_FIELDS, delimiter="\t")
        writer.writeheader()
        written_intervals = 0
        for interval in sorted(interval_values, key=lambda item: (item[0], item[1], item[2])):
            baseline = baselines.get(interval)
            if baseline is None or baseline <= 0:
                continue
            values = interval_values[interval]
            contig, start, end = interval
            row = {
                "CONTIG": contig,
                "START": start,
                "END": end,
                "BASELINE_MEDIAN": baseline,
                "N": len(values),
                "BG_NORM_MEAN": sum(values) / len(values),
                "BG_NORM_MEDIAN": median(values),
                "BG_NORM_SD": stdev(values),
                "BG_NORM_P5": percentile(values, 5),
                "BG_NORM_P95": percentile(values, 95),
            }
            writer.writerow({key: fmt(value) for key, value in row.items()})
            written_intervals += 1

    return written_intervals
