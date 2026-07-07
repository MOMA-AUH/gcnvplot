"""SVG rendering for gcnvplot."""

from __future__ import annotations

import html
import math
from pathlib import Path

from .models import Interval, TranscriptAnnotation
from .utils import fmt, overlaps


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


def _build_plot_svg(
    rows: list[dict[str, object]],
    region: Interval,
    sample_name: str | None = None,
    highlight: Interval | None = None,
    transcript: TranscriptAnnotation | None = None,
) -> str:
    """Build the SVG markup for a plot from prepared rows."""
    base_width = 1400
    height = 620
    left = 95
    top = 72
    bottom = 110 if transcript is None else 180
    panel_width = 304
    panel_gap = 22
    if sample_name is not None:
        estimated_sample_width = math.ceil(len(sample_name) * 7.2)
        panel_width = max(panel_width, 88 + estimated_sample_width + 16)
    plot_right = base_width - 304 - panel_gap
    width = max(base_width, plot_right + panel_gap + panel_width)
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
        ".exon-covered { fill: #1f2933; }",
        ".exon-uncovered { fill: #1f2933; }",
        ".exon-uncovered-marker { fill: #b91c1c; }",
        ".baseline { stroke: #5b6472; stroke-width: 2; stroke-dasharray: 6 5; }",
        ".highlight-band { fill: #f59e0b; opacity: 0.28; stroke: #b45309; stroke-width: 1.1; }",
        ".panel { fill: #ffffff; stroke: #d8dee4; stroke-width: 1.2; }",
        ".panel-title { font-size: 14px; font-weight: 700; }",
        ".panel-label { font-size: 12px; font-weight: 700; }",
        ".panel-text { font-size: 12px; }",
        ".panel-label-emphasis { font-size: 13px; font-weight: 800; }",
        ".panel-text-emphasis { font-size: 13px; font-weight: 700; }",
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
        uncovered_marker_top = exon_top + exon_height + 4.0
        arrow_forward = transcript.strand == "+"
        arrow_exon_gap = 12.0
        arrow_size_x = 6.0
        arrow_size_y = 4.0
        uncovered_marker_size = 6.5
        covered_exon_numbers = {
            exon.number
            for exon in transcript.exons
            if any(overlaps(exon.start, exon.end, int(row["start"]), int(row["end"])) for row in rows)
        }

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

        if highlight_start is not None and highlight_end is not None:
            track_highlight_y = exon_top - 18
            track_highlight_h = exon_height + 30
            band_x = min(highlight_start, highlight_end)
            band_width = abs(highlight_end - highlight_start)
            elements.append(
                f'<rect class="highlight-band" x="{band_x:.2f}" y="{track_highlight_y:.2f}" width="{band_width:.2f}" height="{track_highlight_h:.2f}"/>'
            )

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

        for exon in transcript.exons:
            exon_start = x_for(exon.start)
            exon_end = x_for(exon.end)
            exon_x = min(exon_start, exon_end)
            exon_width = max(1.0, abs(exon_end - exon_start))
            exon_mid = exon_x + exon_width / 2
            exon_class = "exon-covered" if exon.number in covered_exon_numbers else "exon-uncovered"
            elements.extend(
                [
                    f'<rect class="{exon_class}" x="{exon_x:.2f}" y="{exon_top:.2f}" width="{exon_width:.2f}" height="{exon_height:.2f}" rx="3" ry="3"/>',
                ]
            )
            if exon.number not in covered_exon_numbers:
                marker_points = (
                    f"{exon_mid - uncovered_marker_size:.2f},{uncovered_marker_top + uncovered_marker_size:.2f} "
                    f"{exon_mid + uncovered_marker_size:.2f},{uncovered_marker_top + uncovered_marker_size:.2f} "
                    f"{exon_mid:.2f},{uncovered_marker_top:.2f}"
                )
                elements.append(f'<polygon class="exon-uncovered-marker" points="{marker_points}"/>')

        label_base_y = exon_top - 4.0
        label_padding = 3.0
        label_min_gap = 4.0
        label_specs: list[tuple[int, float, float, float, float, int]] = []
        last_exon_number = len(transcript.exons)
        for exon in transcript.exons:
            exon_start = x_for(exon.start)
            exon_end = x_for(exon.end)
            exon_x = min(exon_start, exon_end)
            exon_width = max(1.0, abs(exon_end - exon_start))
            exon_mid = exon_x + exon_width / 2
            label_text = str(exon.number)
            # Use a narrower width estimate so closely spaced exon numbers are
            # only pruned when they would truly collide in the SVG.
            label_width = max(7.0, 5.8 * len(label_text))
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
            label_class = "panel-label-emphasis" if label == "Sample" else "panel-label"
            value_class = "panel-text-emphasis" if label == "Sample" else "panel-text"
            elements.append(
                f'<text class="{label_class}" x="{panel_x + 16}" y="{info_y}">{html.escape(label)}</text>'
            )
            elements.append(f'<text class="{value_class}" x="{panel_x + 88}" y="{info_y}">{html.escape(value)}</text>')
            info_y += 20
        if section_index < len(info_sections) - 1:
            separator_y = info_y - 6
            elements.append(
                f'<line class="panel-separator" x1="{panel_x + 16}" y1="{separator_y:.2f}" x2="{panel_x + panel_width - 16}" y2="{separator_y:.2f}"/>'
            )
            info_y += 14

    legend_y = max(info_y + 18, panel_y + panel_h - 72)
    panel_h = legend_y + 68 - panel_y
    elements[3] = f'<rect class="panel" x="{panel_x}" y="{panel_y}" width="{panel_width}" height="{panel_h}"/>'
    elements.extend(
        [
            f'<circle class="point point-filled" cx="{panel_x + 24}" cy="{legend_y}" r="5"/>',
            f'<text class="panel-text" x="{panel_x + 40}" y="{legend_y + 4}">Overlaps exon</text>',
            f'<circle class="point point-open" cx="{panel_x + 24}" cy="{legend_y + 24}" r="5"/>',
            f'<text class="panel-text" x="{panel_x + 40}" y="{legend_y + 28}">Outside exon</text>',
            f'<polygon class="exon-uncovered-marker" points="{panel_x + 24:.2f},{legend_y + 45.50:.2f} {panel_x + 30.5:.2f},{legend_y + 52.00:.2f} {panel_x + 17.5:.2f},{legend_y + 52.00:.2f}"/>',
            f'<text class="panel-text" x="{panel_x + 40}" y="{legend_y + 52}">Uncovered exon</text>',
            "</svg>",
        ]
    )

    return "\n".join(elements) + "\n"


def write_svg_plot(
    rows: list[dict[str, object]],
    output: Path,
    region: Interval,
    sample_name: str | None = None,
    highlight: Interval | None = None,
    transcript: TranscriptAnnotation | None = None,
) -> None:
    """Write an SVG plot for the sample and background."""
    svg = _build_plot_svg(
        rows,
        region,
        sample_name=sample_name,
        highlight=highlight,
        transcript=transcript,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
