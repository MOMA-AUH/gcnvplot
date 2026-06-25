"""Core data handling and plotting logic for gcnvplot."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import math
from collections import defaultdict
from pathlib import Path

Interval = tuple[str, int, int]

READ_COUNT_FIELDS = ["CONTIG", "START", "END", "COUNT"]
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


def median(values: list[float]) -> float:
    """Return the median of a non-empty list."""
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def percentile(values: list[float], percentile_value: float) -> float:
    """Return a percentile using linear interpolation."""
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = (len(ordered) - 1) * percentile_value / 100
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[int(index)]
    fraction = index - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def stdev(values: list[float]) -> float:
    """Return the sample standard deviation."""
    if len(values) < 2:
        return 0.0
    mean_value = sum(values) / len(values)
    variance = sum((value - mean_value) ** 2 for value in values) / (len(values) - 1)
    return math.sqrt(variance)


def fmt(value: float | str | int) -> str:
    """Format values for TSV or plot labels."""
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "nan"
        return f"{value:.8g}"
    return str(value)


def open_text(path: Path):
    """Open a plain-text or gzipped text file."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(mode="rt", encoding="utf-8")


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Return True when two closed intervals overlap."""
    return a_start <= b_end and b_start <= a_end


def parse_read_counts(path: Path) -> dict[Interval, int]:
    """Parse a GATK CollectReadCounts TSV."""
    counts: dict[Interval, int] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            header = line.rstrip("\n").split("\t")
            if header != READ_COUNT_FIELDS:
                raise ValueError(f"{path}: unexpected read-count header: {header}")
            break
        else:
            raise ValueError(f"{path}: no read-count header found")

        reader = csv.DictReader(handle, fieldnames=READ_COUNT_FIELDS, delimiter="\t")
        for row in reader:
            interval = (row["CONTIG"], int(row["START"]), int(row["END"]))
            counts[interval] = int(row["COUNT"])
    return counts


def read_path_list(path: Path) -> list[Path]:
    """Read a list of sample paths, one per line."""
    paths: list[Path] = []
    base = path.parent
    with path.open(mode="rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            item = Path(line)
            paths.append(item if item.is_absolute() else base / item)
    return paths


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


def parse_region(value: str) -> Interval:
    """Parse a genomic region like chr1:100-200."""
    parts = value.split(":", 1)
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("Expected region like chr1:100-200")
    contig, coords = parts
    bounds = coords.split("-", 1)
    if len(bounds) != 2:
        raise argparse.ArgumentTypeError("Expected region like chr1:100-200")
    start_text, end_text = bounds
    try:
        start = int(start_text.replace(",", ""))
        end = int(end_text.replace(",", ""))
    except ValueError as error:
        raise argparse.ArgumentTypeError("Expected region like chr1:100-200") from error
    if start > end:
        raise argparse.ArgumentTypeError("Region start must be <= end")
    return contig, start, end


def parse_background(path: Path) -> dict[Interval, dict[str, str]]:
    """Parse a background TSV created by create-background."""
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
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
    return background


def build_background(args: argparse.Namespace) -> int:
    """Build interval-wise background summaries."""
    files = read_path_list(args.read_counts_list)
    if not files:
        raise SystemExit(f"{args.read_counts_list}: no read-count TSV paths found")

    sample_counts = [parse_read_counts(path) for path in files]
    baselines = interval_baselines(sample_counts)

    interval_values: dict[Interval, list[float]] = defaultdict(list)
    for counts in sample_counts:
        for interval, value in normalized_counts(counts, baselines).items():
            interval_values[interval].append(value)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        handle.write("# normalization=median-of-ratios\n")
        handle.write("# baseline=median-positive-count\n")
        handle.write(f"# samples={len(files)}\n")
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

    print(f"Background samples: {len(files)}")
    print(f"Background intervals: {written_intervals}")
    print(f"Wrote: {args.output}")
    return 0


def rows_for_plot(args: argparse.Namespace) -> list[dict[str, object]]:
    """Join one sample with the background for plotting."""
    counts = parse_read_counts(args.read_counts)
    background = parse_background(args.background)
    baselines = {
        interval: float(row["BASELINE_MEDIAN"])
        for interval, row in background.items()
        if float(row["BASELINE_MEDIAN"]) > 0
    }
    sample_norm = normalized_counts(counts, baselines)
    region = args.region

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
        row["signal"] = log2_ratio(sample_value, bg_median, 0.01)
        row["signal_lower"] = log2_ratio(bg_lower, bg_median, 0.01)
        row["signal_upper"] = log2_ratio(bg_upper, bg_median, 0.01)
        rows.append(row)

    if missing_background:
        print(f"Skipped intervals missing from background: {missing_background}")
    return rows


def nice_ticks(start: int, end: int, max_ticks: int = 8) -> list[int]:
    """Return visually reasonable genomic coordinate ticks."""
    span = max(1, end - start)
    raw_step = span / max_ticks
    magnitude = 10 ** math.floor(math.log10(raw_step))
    candidates = [1, 2, 5, 10]
    step = min(candidates, key=lambda value: abs(raw_step - value * magnitude)) * magnitude
    first = math.ceil(start / step) * step
    ticks = []
    value = first
    while value <= end:
        ticks.append(int(value))
        value += step
    return ticks


def write_svg_plot(rows: list[dict[str, object]], output: Path, title: str, signal: str, region: Interval) -> None:
    """Write a simple SVG plot for the sample and background."""
    width = 1400
    height = 620
    left = 95
    right = 45
    top = 72
    bottom = 110
    plot_width = width - left - right
    plot_height = height - top - bottom
    region_start = region[1]
    region_end = region[2]

    y_values = [float(row["signal"]) for row in rows]
    band_lower = [float(row["signal_lower"]) for row in rows]
    band_upper = [float(row["signal_upper"]) for row in rows]
    all_y = y_values + band_lower + band_upper
    y_min = min(all_y)
    y_max = max(all_y)
    max_abs = max(abs(y_min), abs(y_max), 0.25)
    y_min = -max_abs
    y_max = max_abs

    def x_for(value: float) -> float:
        return left + (value - region_start) * plot_width / max(1, region_end - region_start)

    def y_for(value: float) -> float:
        return top + (y_max - value) * plot_height / max(1e-12, y_max - y_min)

    points = sorted(rows, key=lambda row: float(row["mid"]))
    upper_points = [(x_for(float(row["mid"])), y_for(float(row["signal_upper"]))) for row in points]
    lower_points = [(x_for(float(row["mid"])), y_for(float(row["signal_lower"]))) for row in reversed(points)]
    ribbon_points = " ".join(f"{x:.2f},{y:.2f}" for x, y in upper_points + lower_points)
    sample_points = " ".join(f"{x_for(float(row['mid'])):.2f},{y_for(float(row['signal'])):.2f}" for row in points)
    baseline_value = 0.0

    y_axis_label = "Log2(sample/background)"
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text { font-family: Arial, Helvetica, sans-serif; fill: #1f2933; }",
        ".axis { stroke: #1f2933; stroke-width: 1.3; }",
        ".grid { stroke: #d8dee4; stroke-width: 1; }",
        ".ribbon { fill: #c8d2df; opacity: 0.85; }",
        ".sample { fill: none; stroke: #127c78; stroke-width: 2.2; }",
        ".point { fill: #127c78; stroke: white; stroke-width: 1.2; }",
        ".baseline { stroke: #5b6472; stroke-width: 2; stroke-dasharray: 6 5; }",
        ".legend-label { font-size: 12px; }",
        "</style>",
        f'<rect width="{width}" height="{height}" fill="white"/>',
        f'<text x="{width / 2}" y="30" text-anchor="middle" font-size="21" font-weight="700">{html.escape(title)}</text>',
        f'<text x="{width / 2}" y="53" text-anchor="middle" font-size="13">{html.escape(y_axis_label)}</text>',
    ]

    for i in range(6):
        value = y_min + i * (y_max - y_min) / 5
        y = y_for(value)
        elements.append(f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{width - right}" y2="{y:.2f}"/>')
        elements.append(f'<text x="{left - 12}" y="{y + 4:.2f}" text-anchor="end" font-size="12">{value:.3g}</text>')

    for tick in nice_ticks(region_start, region_end):
        x = x_for(tick)
        elements.append(f'<line class="grid" x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{height - bottom}" opacity="0.45"/>')
        label_y = height - bottom + 20
        elements.append(
            f'<text x="{x:.2f}" y="{label_y}" text-anchor="end" font-size="11" '
            f'transform="rotate(-35 {x:.2f} {label_y})">{tick:,}</text>'
        )

    elements.extend(
        [
            f'<polygon class="ribbon" points="{ribbon_points}"/>',
            f'<line class="baseline" x1="{left}" y1="{y_for(baseline_value):.2f}" x2="{width - right}" y2="{y_for(baseline_value):.2f}"/>',
            f'<polyline class="sample" points="{sample_points}"/>',
            f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height - bottom}"/>',
            f'<line class="axis" x1="{left}" y1="{height - bottom}" x2="{width - right}" y2="{height - bottom}"/>',
        ]
    )

    for row in points:
        x = x_for(float(row["mid"]))
        y = y_for(float(row["signal"]))
        title_text = (
            f"{row['contig']}:{row['start']}-{row['end']}\n"
            f"COUNT={row['count']}\n"
            f"SAMPLE_NORM={fmt(float(row['sample_norm']))}\n"
            f"BG_MEDIAN_NORM={fmt(float(row['bg_median_norm']))}\n"
            f"SIGNAL={fmt(float(row['signal']))}\n"
            f"BACKGROUND_N={row['background_n']}"
        )
        elements.append(
            f'<circle class="point" cx="{x:.2f}" cy="{y:.2f}" r="4.8">'
            f"<title>{html.escape(title_text)}</title></circle>"
        )

    elements.extend(
        [
            f'<text x="{width / 2}" y="{height - 24}" text-anchor="middle" font-size="14">Genomic coordinate</text>',
            f'<text x="24" y="{top + plot_height / 2}" text-anchor="middle" font-size="14" transform="rotate(-90 24 {top + plot_height / 2})">{html.escape(y_axis_label)}</text>',
            f'<rect x="{width - 322}" y="{top - 8}" width="15" height="15" fill="#c8d2df" opacity="0.85"/>',
            f'<text class="legend-label" x="{width - 300}" y="{top + 4}">background band</text>',
            f'<line x1="{width - 322}" y1="{top + 24}" x2="{width - 304}" y2="{top + 24}" stroke="#5b6472" stroke-width="2" stroke-dasharray="6 5"/>',
            f'<text class="legend-label" x="{width - 300}" y="{top + 28}">expected baseline</text>',
            f'<circle class="point" cx="{width - 314}" cy="{top + 49}" r="5"/>',
            f'<text class="legend-label" x="{width - 300}" y="{top + 53}">sample signal</text>',
            "</svg>",
        ]
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(elements) + "\n", encoding="utf-8")


def plot_sample(args: argparse.Namespace) -> int:
    """Plot one sample against the background."""
    rows = rows_for_plot(args)
    if not rows:
        raise SystemExit("No intervals with background statistics overlap the selected region.")

    write_svg_plot(rows, args.output, "GATK Germline CNV read-count signal", "log2-ratio", args.region)
    print(f"Plotted intervals: {len(rows)}")
    print(f"Wrote: {args.output}")
    return 0
