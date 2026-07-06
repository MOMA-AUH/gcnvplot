"""Core data handling and plotting logic for gcnvplot."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import math
from dataclasses import dataclass
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


@dataclass(frozen=True)
class TranscriptExon:
    """One exon in transcript order."""

    number: int
    start: int
    end: int


@dataclass(frozen=True)
class TranscriptAnnotation:
    """Transcript model used for plotting."""

    transcript_id: str
    gene_name: str
    contig: str
    strand: str
    start: int
    end: int
    exons: list[TranscriptExon]


@dataclass(frozen=True)
class BackgroundSummary:
    """Background intervals plus metadata used in plots."""

    rows: dict[Interval, dict[str, str]]
    lower_percentile: int | None = None
    upper_percentile: int | None = None


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


def parse_gtf_attributes(text: str) -> dict[str, str]:
    """Parse a GTF attribute column into a dictionary."""
    attributes: dict[str, str] = {}
    for item in text.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        key, value = item.split(" ", 1)
        attributes[key] = value.strip().strip('"')
    return attributes


def transcript_id_matches(value: str | None, target: str) -> bool:
    """Return True when two transcript IDs match with or without version suffixes."""
    if value is None:
        return False
    if value == target:
        return True
    return value.split(".", 1)[0] == target.split(".", 1)[0]


def load_transcript_annotation(gtf_path: Path, transcript_id: str) -> TranscriptAnnotation:
    """Load exon coordinates and gene name for one transcript from a GTF."""
    contig: str | None = None
    strand: str | None = None
    gene_name: str | None = None
    gene_id: str | None = None
    exon_coords: list[tuple[int, int]] = []

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqname, _source, feature, start_text, end_text, _score, feature_strand, _frame, attributes_text = fields
            attributes = parse_gtf_attributes(attributes_text)
            if not transcript_id_matches(attributes.get("transcript_id"), transcript_id):
                continue

            feature_contig = seqname
            feature_strand = feature_strand
            if contig is None:
                contig = feature_contig
            elif contig != feature_contig:
                raise ValueError(f"{gtf_path}: transcript {transcript_id} spans multiple contigs")

            if strand is None:
                strand = feature_strand
            elif strand != feature_strand:
                raise ValueError(f"{gtf_path}: transcript {transcript_id} has inconsistent strands")

            if gene_name is None:
                gene_name = attributes.get("gene_name") or attributes.get("gene_id")
            if gene_id is None:
                gene_id = attributes.get("gene_id")

            if feature == "exon":
                exon_coords.append((int(start_text), int(end_text)))

    if not exon_coords:
        raise ValueError(f"{gtf_path}: transcript {transcript_id} not found or has no exons")
    if contig is None or strand is None:
        raise ValueError(f"{gtf_path}: transcript {transcript_id} is missing contig or strand information")

    ordered_coords = sorted(exon_coords, key=lambda item: item[0], reverse=strand == "-")
    exons = [TranscriptExon(number=index + 1, start=start, end=end) for index, (start, end) in enumerate(ordered_coords)]
    transcript_gene_name = gene_name or gene_id or transcript_id
    start = min(exon.start for exon in exons)
    end = max(exon.end for exon in exons)
    return TranscriptAnnotation(
        transcript_id=transcript_id,
        gene_name=transcript_gene_name,
        contig=contig,
        strand=strand,
        start=start,
        end=end,
        exons=exons,
    )


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


def rows_for_plot(
    args: argparse.Namespace,
    region: Interval,
    background_summary: BackgroundSummary,
    transcript: TranscriptAnnotation | None = None,
) -> list[dict[str, object]]:
    """Join one sample with the background for plotting."""
    counts = parse_read_counts(args.read_counts)
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


def write_svg_plot(
    rows: list[dict[str, object]],
    output: Path,
    title: str,
    signal: str,
    region: Interval,
    background_summary: BackgroundSummary,
    sample_name: str | None = None,
    highlight: Interval | None = None,
    transcript: TranscriptAnnotation | None = None,
) -> None:
    """Write a simple SVG plot for the sample and background."""
    width = 1400
    height = 620
    left = 95
    top = 72
    bottom = 110 if transcript is None else 180
    panel_width = 304
    panel_gap = 22
    plot_right = width - panel_width - panel_gap
    plot_width = plot_right - left
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

    highlight_start: float | None = None
    highlight_end: float | None = None
    if highlight is not None and highlight[0] == region[0]:
        clipped_start = max(region_start, highlight[1])
        clipped_end = min(region_end, highlight[2])
        if clipped_start < clipped_end:
            highlight_start = x_for(float(clipped_start))
            highlight_end = x_for(float(clipped_end))

    y_axis_label = "log2(sample/background)"
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text { font-family: Arial, Helvetica, sans-serif; fill: #1f2933; }",
        ".axis { stroke: #1f2933; stroke-width: 1.3; }",
        ".grid { stroke: #d8dee4; stroke-width: 1; }",
        ".ribbon { fill: #c8d2df; opacity: 0.85; }",
        ".sample { fill: none; stroke: #127c78; stroke-width: 2.2; }",
        ".point { stroke-width: 1.6; }",
        ".point-filled { fill: #127c78; stroke: white; }",
        ".point-open { fill: white; stroke: #127c78; }",
        ".baseline { stroke: #5b6472; stroke-width: 2; stroke-dasharray: 6 5; }",
        ".highlight-band { fill: #f59e0b; opacity: 0.28; stroke: #b45309; stroke-width: 1.1; }",
        ".panel { fill: #ffffff; stroke: #d8dee4; stroke-width: 1.2; }",
        ".panel-title { font-size: 14px; font-weight: 700; }",
        ".panel-label { font-size: 12px; font-weight: 700; }",
        ".panel-text { font-size: 12px; }",
        ".panel-separator { stroke: #d8dee4; stroke-width: 1.2; }",
        "</style>",
        f'<rect width="{width}" height="{height}" fill="white"/>',
    ]

    if highlight_start is not None and highlight_end is not None:
        band_x = min(highlight_start, highlight_end)
        band_width = abs(highlight_end - highlight_start)
        elements.append(
            f'<rect class="highlight-band" x="{band_x:.2f}" y="{top}" width="{band_width:.2f}" height="{height - bottom - top:.2f}"/>'
        )

    for i in range(6):
        value = y_min + i * (y_max - y_min) / 5
        y = y_for(value)
        elements.append(f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{plot_right}" y2="{y:.2f}"/>')
        elements.append(f'<text x="{left - 12}" y="{y + 4:.2f}" text-anchor="end" font-size="12">{value:.3g}</text>')

    for tick in nice_ticks(region_start, region_end):
        x = x_for(tick)
        elements.append(f'<line class="grid" x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{height - bottom}" opacity="0.45"/>')
        tick_y1 = height - bottom
        tick_y2 = tick_y1 + 7
        label_y = tick_y2 + 13
        elements.extend(
            [
                f'<line class="x-axis-tick" x1="{x:.2f}" y1="{tick_y1:.2f}" x2="{x:.2f}" y2="{tick_y2:.2f}" stroke="#1f2933" stroke-width="1.1"/>',
                f'<text class="x-axis-tick-label" x="{x:.2f}" y="{label_y:.2f}" text-anchor="middle" font-size="11">{tick:,}</text>',
            ]
        )

    elements.extend(
        [
            f'<polygon class="ribbon" points="{ribbon_points}"/>',
            f'<line class="baseline" x1="{left}" y1="{y_for(baseline_value):.2f}" x2="{plot_right}" y2="{y_for(baseline_value):.2f}"/>',
            f'<polyline class="sample" points="{sample_points}"/>',
            f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{height - bottom}"/>',
            f'<line class="axis" x1="{left}" y1="{height - bottom}" x2="{plot_right}" y2="{height - bottom}"/>',
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
        point_class = "point point-filled" if bool(row["overlaps_exon"]) else "point point-open"
        elements.append(
            f'<circle class="{point_class}" cx="{x:.2f}" cy="{y:.2f}" r="4.8">'
            f"<title>{html.escape(title_text)}</title></circle>"
        )

    if transcript is not None:
        track_y = height - 95
        exon_height = 18
        exon_top = track_y - exon_height / 2
        arrow_forward = transcript.strand == "+"
        arrow_exon_gap = 12.0
        arrow_size_x = 6.0
        arrow_size_y = 4.0

        def add_arrow(x_pos: float) -> None:
            if arrow_forward:
                points = (
                    f"{x_pos - arrow_size_x:.2f},{track_y - arrow_size_y:.2f} "
                    f"{x_pos + arrow_size_x:.2f},{track_y:.2f} "
                    f"{x_pos - arrow_size_x:.2f},{track_y + arrow_size_y:.2f}"
                )
            else:
                points = (
                    f"{x_pos + arrow_size_x:.2f},{track_y - arrow_size_y:.2f} "
                    f"{x_pos - arrow_size_x:.2f},{track_y:.2f} "
                    f"{x_pos + arrow_size_x:.2f},{track_y + arrow_size_y:.2f}"
                )
            elements.append(f'<polygon points="{points}" fill="#5b6472"/>')

        genomic_exons = sorted(transcript.exons, key=lambda exon: exon.start)
        spans = [
            (left.end, right.start)
            for left, right in zip(genomic_exons, genomic_exons[1:])
            if right.start > left.end
        ]

        for span_start, span_end in spans:
            x_start = x_for(float(span_start))
            x_end = x_for(float(span_end))
            if x_end < x_start:
                x_start, x_end = x_end, x_start
            elements.append(
                f'<line x1="{x_start:.2f}" y1="{track_y}" x2="{x_end:.2f}" y2="{track_y}" stroke="#5b6472" stroke-width="2"/>'
            )
            arrow_start = x_start + arrow_exon_gap + arrow_size_x
            arrow_end = x_end - arrow_exon_gap - arrow_size_x
            if arrow_end <= arrow_start:
                continue
            span_width = abs(arrow_end - arrow_start)
            arrow_count = max(1, int(span_width / 90))
            for index in range(arrow_count):
                fraction = (index + 1) / (arrow_count + 1)
                add_arrow(arrow_start + fraction * (arrow_end - arrow_start))

        if highlight_start is not None and highlight_end is not None:
            track_highlight_y = exon_top - 18
            track_highlight_h = exon_height + 30
            band_x = min(highlight_start, highlight_end)
            band_width = abs(highlight_end - highlight_start)
            elements.append(
                f'<rect class="highlight-band" x="{band_x:.2f}" y="{track_highlight_y:.2f}" width="{band_width:.2f}" height="{track_highlight_h:.2f}"/>'
            )

        for exon in transcript.exons:
            exon_start = x_for(exon.start)
            exon_end = x_for(exon.end)
            exon_x = min(exon_start, exon_end)
            exon_width = max(1.0, abs(exon_end - exon_start))
            elements.extend(
                [
                    f'<rect x="{exon_x:.2f}" y="{exon_top:.2f}" width="{exon_width:.2f}" height="{exon_height:.2f}" rx="3" ry="3" fill="#1f2933"/>',
                ]
            )

        label_base_y = exon_top - 4.0
        label_padding = 4.0
        label_min_gap = 6.0
        label_specs: list[tuple[int, float, float, float, float, int]] = []
        last_exon_number = len(transcript.exons)
        for exon in transcript.exons:
            exon_start = x_for(exon.start)
            exon_end = x_for(exon.end)
            exon_x = min(exon_start, exon_end)
            exon_width = max(1.0, abs(exon_end - exon_start))
            exon_mid = exon_x + exon_width / 2
            label_text = str(exon.number)
            label_width = max(8.0, 6.5 * len(label_text))
            label_left = exon_mid - label_width / 2 - label_padding
            label_right = exon_mid + label_width / 2 + label_padding
            priority = 3 if exon.number in {1, last_exon_number} else 1
            label_specs.append((exon.number, exon_mid, label_width, label_left, label_right, priority))

        label_specs.sort(key=lambda spec: (spec[4], spec[3], spec[0]))
        compatible_indices: list[int] = []
        for index, (_exon_number, _exon_mid, _label_width, label_left, _label_right, _priority) in enumerate(label_specs):
            compatible_index = -1
            for previous_index in range(index - 1, -1, -1):
                previous_right = label_specs[previous_index][4]
                if previous_right + label_min_gap < label_left:
                    compatible_index = previous_index
                    break
            compatible_indices.append(compatible_index)

        best_scores: list[float] = [0.0] * len(label_specs)
        take_label: list[bool] = [False] * len(label_specs)
        for index, (_exon_number, _exon_mid, _label_width, _label_left, _label_right, priority) in enumerate(label_specs):
            include_score = priority + (best_scores[compatible_indices[index]] if compatible_indices[index] >= 0 else 0.0)
            exclude_score = best_scores[index - 1] if index > 0 else 0.0
            if include_score > exclude_score:
                best_scores[index] = include_score
                take_label[index] = True
            else:
                best_scores[index] = exclude_score

        label_positions: list[tuple[int, float, float]] = []
        index = len(label_specs) - 1
        while index >= 0:
            if not take_label[index]:
                index -= 1
                continue
            exon_number, exon_mid, _label_width, _label_left, _label_right, _priority = label_specs[index]
            label_positions.append((exon_number, exon_mid, label_base_y))
            index = compatible_indices[index]
        label_positions.reverse()

        for exon_number, label_x, label_y in label_positions:
            elements.append(
                f'<text x="{label_x:.2f}" y="{label_y:.2f}" text-anchor="middle" font-size="11" font-weight="700">{exon_number}</text>'
            )

    elements.extend(
        [
            f'<text x="24" y="{top + plot_height / 2}" text-anchor="middle" font-size="14" transform="rotate(-90 24 {top + plot_height / 2})">{html.escape(y_axis_label)}</text>',
        ]
    )

    info_sections: list[list[tuple[str, str]]] = []
    if sample_name is not None:
        info_sections.append([("Sample", sample_name)])

    main_section: list[tuple[str, str]] = []
    if transcript is not None:
        main_section.extend(
            [
                ("Gene", transcript.gene_name),
                ("Transcript", transcript.transcript_id),
            ]
        )
    main_section.append(("Region", f"{region[0]}:{region[1]:,}-{region[2]:,}"))
    info_sections.append(main_section)

    detail_section: list[tuple[str, str]] = []
    if highlight_start is not None and highlight_end is not None:
        detail_section.append(("Highlight", f"{highlight[0]}:{highlight[1]:,}-{highlight[2]:,}"))
    if transcript is not None and highlight is not None and highlight[0] == region[0]:
        visible_highlight_start = max(region[1], highlight[1])
        visible_highlight_end = min(region[2], highlight[2])
        if visible_highlight_start < visible_highlight_end:
            covered_exons = [
                exon.number
                for exon in transcript.exons
                if overlaps(exon.start, exon.end, visible_highlight_start, visible_highlight_end)
            ]
            if covered_exons:
                runs: list[str] = []
                run_start = covered_exons[0]
                run_prev = covered_exons[0]
                for exon_number in covered_exons[1:]:
                    if exon_number == run_prev + 1:
                        run_prev = exon_number
                        continue
                    runs.append(str(run_start) if run_start == run_prev else f"{run_start}-{run_prev}")
                    run_start = run_prev = exon_number
                runs.append(str(run_start) if run_start == run_prev else f"{run_start}-{run_prev}")
                detail_section.append(("Exons", ",".join(runs)))
    if detail_section:
        info_sections.append(detail_section)

    panel_x = plot_right + panel_gap
    panel_y = top - 8
    info_row_count = sum(len(section) for section in info_sections)
    separator_count = max(0, len(info_sections) - 1)
    panel_h = max(132, 56 + info_row_count * 20 + separator_count * 14 + 72)
    elements[3:3] = [
        f'<rect class="panel" x="{panel_x}" y="{panel_y}" width="{panel_width}" height="{panel_h}"/>',
        f'<text class="panel-title" x="{panel_x + 16}" y="{panel_y + 24}">Plot Info</text>',
    ]

    info_y = panel_y + 49
    for section_index, section in enumerate(info_sections):
        for label, value in section:
            elements.append(f'<text class="panel-label" x="{panel_x + 16}" y="{info_y}">{html.escape(label)}</text>')
            elements.append(f'<text class="panel-text" x="{panel_x + 88}" y="{info_y}">{html.escape(value)}</text>')
            info_y += 20
        if section_index < len(info_sections) - 1:
            separator_y = info_y - 6
            elements.append(
                f'<line class="panel-separator" x1="{panel_x + 16}" y1="{separator_y:.2f}" x2="{panel_x + panel_width - 16}" y2="{separator_y:.2f}"/>'
            )
            info_y += 14

    legend_y = max(info_y + 18, panel_y + panel_h - 58)
    panel_h = legend_y + 54 - panel_y
    elements[3] = f'<rect class="panel" x="{panel_x}" y="{panel_y}" width="{panel_width}" height="{panel_h}"/>'
    elements.extend(
        [
            f'<circle class="point point-filled" cx="{panel_x + 24}" cy="{legend_y}" r="5"/>',
            f'<text class="panel-text" x="{panel_x + 40}" y="{legend_y + 4}">Overlaps exon</text>',
            f'<circle class="point point-open" cx="{panel_x + 24}" cy="{legend_y + 24}" r="5"/>',
            f'<text class="panel-text" x="{panel_x + 40}" y="{legend_y + 28}">Outside exon</text>',
            "</svg>",
        ]
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(elements) + "\n", encoding="utf-8")


def plot_sample(args: argparse.Namespace) -> int:
    """Plot one sample against the background."""
    transcript: TranscriptAnnotation | None = None
    region = args.region
    if args.transcript is not None:
        transcript = load_transcript_annotation(args.gtf, args.transcript)
        region = (transcript.contig, transcript.start, transcript.end)

    background_summary = parse_background(args.background)
    rows = rows_for_plot(args, region, background_summary, transcript=transcript)
    if not rows:
        raise SystemExit("No intervals with background statistics overlap the selected region.")

    write_svg_plot(
        rows,
        args.output,
        "Normalized GATK GermlineCNVCaller Read Counts",
        "log2-ratio",
        region,
        background_summary=background_summary,
        sample_name=args.sample_name,
        highlight=args.highlight,
        transcript=transcript,
    )
    print(f"Plotted intervals: {len(rows)}")
    print(f"Wrote: {args.output}")
    return 0
