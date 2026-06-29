"""CLI tests for gcnvplot."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pytest

import gcnvplot
import gcnvplot.cli as cli


def write_read_counts(path: Path, rows: list[tuple[str, int, int, int]]) -> None:
    """Write a tiny CollectReadCounts-style TSV."""
    lines = ["@RG\tID:test", "CONTIG\tSTART\tEND\tCOUNT"]
    lines.extend(f"{contig}\t{start}\t{end}\t{count}" for contig, start, end, count in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_gtf_gz(path: Path, lines: list[str]) -> None:
    """Write a tiny gzipped GTF file."""
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def svg_polygon_x_values(svg: str) -> list[float]:
    """Return all x coordinates from transcript arrowhead polygon point lists."""
    values: list[float] = []
    for points in re.findall(r'<polygon points="([^"]+)" fill="#5b6472"', svg):
        for point in points.split():
            x_text, _y_text = point.split(",", 1)
            values.append(float(x_text))
    return values


def svg_exon_x_ranges(svg: str) -> list[tuple[float, float]]:
    """Return x ranges for transcript exon rectangles."""
    ranges: list[tuple[float, float]] = []
    pattern = r'<rect x="([0-9.]+)" y="[^"]+" width="([0-9.]+)" height="18.00"[^>]+fill="#1f2933"'
    for x_text, width_text in re.findall(pattern, svg):
        start = float(x_text)
        ranges.append((start, start + float(width_text)))
    return ranges


def svg_exon_label_positions(svg: str) -> list[tuple[int, float, float]]:
    """Return exon number labels as (number, x, y) tuples."""
    positions: list[tuple[int, float, float]] = []
    pattern = (
        r'<text x="([0-9.]+)" y="([0-9.]+)" '
        r'text-anchor="middle" font-size="11" font-weight="700">([0-9]+)</text>'
    )
    for x_text, y_text, label_text in re.findall(pattern, svg):
        positions.append((int(label_text), float(x_text), float(y_text)))
    return positions


def svg_x_axis_ticks(svg: str) -> list[float]:
    """Return x positions of genomic coordinate tick marks."""
    ticks: list[float] = []
    pattern = r'<line class="x-axis-tick" x1="([0-9.]+)" y1="([0-9.]+)" x2="([0-9.]+)" y2="([0-9.]+)" stroke="#1f2933" stroke-width="1.1"/>'
    for x1_text, _y1_text, x2_text, _y2_text in re.findall(pattern, svg):
        ticks.append(float(x1_text))
        assert float(x1_text) == float(x2_text)
    return ticks


def svg_highlight_bands(svg: str) -> list[tuple[float, float, float, float]]:
    """Return highlight band rectangles as (x, y, width, height)."""
    bands: list[tuple[float, float, float, float]] = []
    pattern = r'<rect class="highlight-band" x="([0-9.]+)" y="([0-9.]+)" width="([0-9.]+)" height="([0-9.]+)"/>'
    for x_text, y_text, width_text, height_text in re.findall(pattern, svg):
        bands.append((float(x_text), float(y_text), float(width_text), float(height_text)))
    return bands


def test_main_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--version"])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out == f"gcnvplot {gcnvplot.__version__}\n"


def test_main_requires_subcommand(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main([])

    assert exc_info.value.code == 2
    assert "required" in capsys.readouterr().err


def test_build_background_writes_expected_summary(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    sample_a = tmp_path / "sample_a.tsv"
    sample_b = tmp_path / "sample_b.tsv"
    write_read_counts(
        sample_a,
        [
            ("chr1", 100, 199, 10),
            ("chr1", 200, 299, 20),
        ],
    )
    write_read_counts(
        sample_b,
        [
            ("chr1", 100, 199, 15),
            ("chr1", 200, 299, 30),
        ],
    )
    paths_file = tmp_path / "background_inputs.txt"
    paths_file.write_text("sample_a.tsv\nsample_b.tsv\n", encoding="utf-8")
    output_path = tmp_path / "background.tsv"

    assert (
        cli.main(
            [
                "create-background",
                "--read-counts-list",
                str(paths_file),
                "--output",
                str(output_path),
            ]
        )
        == 0
    )

    background_text = output_path.read_text(encoding="utf-8")
    assert "# normalization=median-of-ratios" in background_text
    assert "chr1\t100\t199\t12.5\t2\t12.5\t12.5\t0\t12.5\t12.5" in background_text
    assert "chr1\t200\t299\t25\t2\t25\t25\t0\t25\t25" in background_text

    assert capsys.readouterr().out == (
        "Background samples: 2\n"
        "Background intervals: 2\n"
        f"Wrote: {output_path}\n"
    )


def test_create_background_requires_output(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    paths_file = tmp_path / "background_inputs.txt"
    paths_file.write_text("sample_a.tsv\n", encoding="utf-8")

    with pytest.raises(SystemExit) as exc_info:
        cli.main(["create-background", "--read-counts-list", str(paths_file)])

    assert exc_info.value.code == 2
    assert "the following arguments are required: --output" in capsys.readouterr().err


def test_plot_requires_output(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    background_path = tmp_path / "background.tsv"
    background_path.write_text(
        "\n".join(
            [
                "# normalization=median",
                "# samples=2",
                "# lower_percentile=5",
                "# upper_percentile=95",
                "CONTIG\tSTART\tEND\tN\tBG_NORM_MEAN\tBG_NORM_MEDIAN\tBG_NORM_SD\tBG_NORM_P5\tBG_NORM_P95",
                "chr1\t100\t199\t2\t1\t1\t0.1\t0.8\t1.2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sample_path = tmp_path / "sample.tsv"
    write_read_counts(sample_path, [("chr1", 100, 199, 10)])

    with pytest.raises(SystemExit) as exc_info:
        cli.main(
            [
                "plot",
                "--read-counts",
                str(sample_path),
                "--background",
                str(background_path),
            ]
        )

    assert exc_info.value.code == 2
    assert "the following arguments are required: --output" in capsys.readouterr().err


def test_plot_log2_ratio_writes_svg_with_zero_centered_axis(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    background_path = tmp_path / "background.tsv"
    background_path.write_text(
        "\n".join(
            [
                "# normalization=median-of-ratios",
                "# baseline=median-positive-count",
                "# samples=2",
                "# lower_percentile=5",
                "# upper_percentile=95",
                "CONTIG\tSTART\tEND\tBASELINE_MEDIAN\tN\tBG_NORM_MEAN\tBG_NORM_MEDIAN\tBG_NORM_SD\tBG_NORM_P5\tBG_NORM_P95",
                "chr1\t100\t199\t12.5\t2\t12.5\t12.5\t0\t12.5\t12.5",
                "chr2\t100\t199\t0\t2\t1\t1\t0.1\t0.8\t1.2",
                "chr1\t200\t299\t25\t2\t25\t25\t0\t25\t25",
                "chr2\t200\t299\t0\t2\t1\t1\t0.1\t0.8\t1.2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sample_path = tmp_path / "sample.tsv"
    write_read_counts(
        sample_path,
        [
            ("chr1", 100, 199, 10),
            ("chr2", 100, 199, 0),
            ("chr1", 200, 299, 40),
            ("chr2", 200, 299, 0),
        ],
    )
    output_path = tmp_path / "plot.svg"

    assert (
        cli.main(
            [
                "plot",
                "--read-counts",
                str(sample_path),
                "--background",
                str(background_path),
                "--region",
                "chr1:100-299",
                "--highlight",
                "chr1:150-250",
                "--output",
                str(output_path),
            ]
        )
        == 0
    )

    svg = output_path.read_text(encoding="utf-8")
    assert "chr1:100-199" in svg
    assert "chr1:200-299" in svg
    assert "chr2:100-199" not in svg
    assert "chr2:200-299" not in svg
    assert "Log2(sample/background)" in svg
    assert "SIGNAL=-0.584386" in svg
    assert "SIGNAL=0.41489328" in svg
    highlight_bands = svg_highlight_bands(svg)
    assert len(highlight_bands) == 1

    assert capsys.readouterr().out == (
        "Plotted intervals: 2\n"
        f"Wrote: {output_path}\n"
    )


def test_plot_transcript_writes_exon_track_and_gene_name(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    background_path = tmp_path / "background.tsv"
    background_path.write_text(
        "\n".join(
            [
                "# normalization=median-of-ratios",
                "# baseline=median-positive-count",
                "# samples=2",
                "# lower_percentile=5",
                "# upper_percentile=95",
                "CONTIG\tSTART\tEND\tBASELINE_MEDIAN\tN\tBG_NORM_MEAN\tBG_NORM_MEDIAN\tBG_NORM_SD\tBG_NORM_P5\tBG_NORM_P95",
                "chr1\t100\t199\t12.5\t2\t12.5\t12.5\t0\t12.5\t12.5",
                "chr1\t200\t299\t25\t2\t25\t25\t0\t25\t25",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sample_path = tmp_path / "sample.tsv"
    write_read_counts(
        sample_path,
        [
            ("chr1", 100, 199, 10),
            ("chr1", 200, 299, 40),
        ],
    )
    gtf_path = tmp_path / "transcript.gtf.gz"
    write_gtf_gz(
        gtf_path,
        [
            'chr1\tsource\ttranscript\t100\t220\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "TX1";',
            'chr1\tsource\texon\t100\t120\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "TX1"; exon_number "1";',
            'chr1\tsource\texon\t180\t220\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "TX1"; exon_number "2";',
        ],
    )
    output_path = tmp_path / "plot.svg"

    assert (
        cli.main(
            [
                "plot",
                "--read-counts",
                str(sample_path),
                "--background",
                str(background_path),
                "--transcript",
                "TX1",
                "--gtf",
                str(gtf_path),
                "--highlight",
                "chr1:150-250",
                "--output",
                str(output_path),
            ]
        )
        == 0
    )

    svg = output_path.read_text(encoding="utf-8")
    assert "MYGENE" in svg
    assert ">TX1</text>" in svg
    assert "Genomic coordinate" not in svg
    assert "5&#x27;" not in svg
    assert "transcript-arrow-forward" not in svg
    assert "rotate(-35" not in svg
    assert svg_x_axis_ticks(svg)
    highlight_bands = svg_highlight_bands(svg)
    assert len(highlight_bands) == 2
    assert ">1</text>" in svg
    assert ">2</text>" in svg
    assert "chr1:100-199" in svg
    assert "chr1:200-299" in svg

    assert capsys.readouterr().out == (
        "Plotted intervals: 2\n"
        f"Wrote: {output_path}\n"
    )


def test_plot_reverse_transcript_arrows_do_not_touch_exons(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    background_path = tmp_path / "background.tsv"
    background_path.write_text(
        "\n".join(
            [
                "# normalization=median-of-ratios",
                "# baseline=median-positive-count",
                "# samples=2",
                "# lower_percentile=5",
                "# upper_percentile=95",
                "CONTIG\tSTART\tEND\tBASELINE_MEDIAN\tN\tBG_NORM_MEAN\tBG_NORM_MEDIAN\tBG_NORM_SD\tBG_NORM_P5\tBG_NORM_P95",
                "chr1\t100\t199\t12.5\t2\t12.5\t12.5\t0\t12.5\t12.5",
                "chr1\t200\t299\t25\t2\t25\t25\t0\t25\t25",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sample_path = tmp_path / "sample.tsv"
    write_read_counts(
        sample_path,
        [
            ("chr1", 100, 199, 10),
            ("chr1", 200, 299, 40),
        ],
    )
    gtf_path = tmp_path / "reverse_transcript.gtf.gz"
    write_gtf_gz(
        gtf_path,
        [
            'chr1\tsource\ttranscript\t100\t220\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1";',
            'chr1\tsource\texon\t100\t120\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1"; exon_number "2";',
            'chr1\tsource\texon\t180\t220\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1"; exon_number "1";',
        ],
    )
    output_path = tmp_path / "plot.svg"

    assert (
        cli.main(
            [
                "plot",
                "--read-counts",
                str(sample_path),
                "--background",
                str(background_path),
                "--transcript",
                "RTX1",
                "--gtf",
                str(gtf_path),
                "--output",
                str(output_path),
            ]
        )
        == 0
    )

    svg = output_path.read_text(encoding="utf-8")
    arrow_x_values = svg_polygon_x_values(svg)
    exon_x_ranges = svg_exon_x_ranges(svg)
    assert arrow_x_values
    assert len(exon_x_ranges) == 2
    for arrow_x in arrow_x_values:
        assert all(not start <= arrow_x <= end for start, end in exon_x_ranges)

    assert capsys.readouterr().out == (
        "Plotted intervals: 2\n"
        f"Wrote: {output_path}\n"
    )


def test_plot_transcript_keeps_full_transcript_width_when_exons_are_uncovered(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    background_path = tmp_path / "background.tsv"
    background_path.write_text(
        "\n".join(
            [
                "# normalization=median-of-ratios",
                "# baseline=median-positive-count",
                "# samples=2",
                "# lower_percentile=5",
                "# upper_percentile=95",
                "CONTIG\tSTART\tEND\tBASELINE_MEDIAN\tN\tBG_NORM_MEAN\tBG_NORM_MEDIAN\tBG_NORM_SD\tBG_NORM_P5\tBG_NORM_P95",
                "chr1\t100\t199\t12.5\t2\t12.5\t12.5\t0\t12.5\t12.5",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sample_path = tmp_path / "sample.tsv"
    write_read_counts(
        sample_path,
        [
            ("chr1", 100, 199, 10),
        ],
    )
    gtf_path = tmp_path / "sparse_transcript.gtf.gz"
    write_gtf_gz(
        gtf_path,
        [
            'chr1\tsource\ttranscript\t100\t400\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "SPARSE1";',
            'chr1\tsource\texon\t100\t120\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "SPARSE1"; exon_number "1";',
            'chr1\tsource\texon\t350\t400\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "SPARSE1"; exon_number "2";',
        ],
    )
    output_path = tmp_path / "plot.svg"

    assert (
        cli.main(
            [
                "plot",
                "--read-counts",
                str(sample_path),
                "--background",
                str(background_path),
                "--transcript",
                "SPARSE1",
                "--gtf",
                str(gtf_path),
                "--output",
                str(output_path),
            ]
        )
        == 0
    )

    svg = output_path.read_text(encoding="utf-8")
    exon_x_ranges = svg_exon_x_ranges(svg)
    assert len(exon_x_ranges) == 2
    assert exon_x_ranges[0][0] >= 95.0
    assert exon_x_ranges[1][1] <= 1355.0
    assert exon_x_ranges[1][0] > exon_x_ranges[0][1]

    assert capsys.readouterr().out == (
        "Plotted intervals: 1\n"
        f"Wrote: {output_path}\n"
    )


def test_plot_transcript_prunes_overlapping_exon_labels(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    background_path = tmp_path / "background.tsv"
    background_path.write_text(
        "\n".join(
            [
                "# normalization=median-of-ratios",
                "# baseline=median-positive-count",
                "# samples=2",
                "# lower_percentile=5",
                "# upper_percentile=95",
                "CONTIG\tSTART\tEND\tBASELINE_MEDIAN\tN\tBG_NORM_MEAN\tBG_NORM_MEDIAN\tBG_NORM_SD\tBG_NORM_P5\tBG_NORM_P95",
                "chr1\t100\t199\t12.5\t2\t12.5\t12.5\t0\t12.5\t12.5",
                "chr1\t200\t299\t25\t2\t25\t25\t0\t25\t25",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sample_path = tmp_path / "sample.tsv"
    write_read_counts(
        sample_path,
        [
            ("chr1", 100, 199, 10),
            ("chr1", 200, 299, 40),
        ],
    )
    gtf_path = tmp_path / "close_exons.gtf.gz"
    write_gtf_gz(
        gtf_path,
        [
            'chr1\tsource\ttranscript\t100\t10000\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "CLOSE1";',
            'chr1\tsource\texon\t100\t110\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "CLOSE1"; exon_number "1";',
            'chr1\tsource\texon\t111\t121\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "CLOSE1"; exon_number "2";',
            'chr1\tsource\texon\t122\t132\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "CLOSE1"; exon_number "3";',
            'chr1\tsource\texon\t133\t143\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "CLOSE1"; exon_number "4";',
            'chr1\tsource\texon\t144\t154\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "CLOSE1"; exon_number "5";',
            'chr1\tsource\texon\t155\t165\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "CLOSE1"; exon_number "6";',
            'chr1\tsource\texon\t9950\t10000\t.\t+\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "CLOSE1"; exon_number "7";',
        ],
    )
    output_path = tmp_path / "plot.svg"

    assert (
        cli.main(
            [
                "plot",
                "--read-counts",
                str(sample_path),
                "--background",
                str(background_path),
                "--transcript",
                "CLOSE1",
                "--gtf",
                str(gtf_path),
                "--output",
                str(output_path),
            ]
        )
        == 0
    )

    svg = output_path.read_text(encoding="utf-8")
    exon_label_positions = svg_exon_label_positions(svg)
    assert 1 <= len(exon_label_positions) < 7
    xs = sorted(x for _number, x, _y in exon_label_positions)
    assert all(right - left >= 12.0 for left, right in zip(xs, xs[1:]))

    assert capsys.readouterr().out == (
        "Plotted intervals: 2\n"
        f"Wrote: {output_path}\n"
    )


def test_plot_reverse_transcript_prunes_overlapping_exon_labels(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    background_path = tmp_path / "background.tsv"
    background_path.write_text(
        "\n".join(
            [
                "# normalization=median-of-ratios",
                "# baseline=median-positive-count",
                "# samples=2",
                "# lower_percentile=5",
                "# upper_percentile=95",
                "CONTIG\tSTART\tEND\tBASELINE_MEDIAN\tN\tBG_NORM_MEAN\tBG_NORM_MEDIAN\tBG_NORM_SD\tBG_NORM_P5\tBG_NORM_P95",
                "chr1\t100\t199\t12.5\t2\t12.5\t12.5\t0\t12.5\t12.5",
                "chr1\t200\t299\t25\t2\t25\t25\t0\t25\t25",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    sample_path = tmp_path / "sample.tsv"
    write_read_counts(
        sample_path,
        [
            ("chr1", 100, 199, 10),
            ("chr1", 200, 299, 40),
        ],
    )
    gtf_path = tmp_path / "reverse_pairs.gtf.gz"
    write_gtf_gz(
        gtf_path,
        [
            'chr1\tsource\ttranscript\t100\t10000\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1";',
            'chr1\tsource\texon\t9950\t10000\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1"; exon_number "1";',
            'chr1\tsource\texon\t100\t110\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1"; exon_number "6";',
            'chr1\tsource\texon\t111\t121\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1"; exon_number "5";',
            'chr1\tsource\texon\t122\t132\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1"; exon_number "4";',
            'chr1\tsource\texon\t133\t143\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1"; exon_number "3";',
            'chr1\tsource\texon\t144\t154\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1"; exon_number "2";',
            'chr1\tsource\texon\t155\t165\t.\t-\t.\tgene_id "GENE1"; gene_name "MYGENE"; transcript_id "RTX1"; exon_number "7";',
        ],
    )
    output_path = tmp_path / "plot.svg"

    assert (
        cli.main(
            [
                "plot",
                "--read-counts",
                str(sample_path),
                "--background",
                str(background_path),
                "--transcript",
                "RTX1",
                "--gtf",
                str(gtf_path),
                "--output",
                str(output_path),
            ]
        )
        == 0
    )

    svg = output_path.read_text(encoding="utf-8")
    exon_label_positions = svg_exon_label_positions(svg)
    assert 1 <= len(exon_label_positions) < 7
    xs = sorted(x for _number, x, _y in exon_label_positions)
    assert all(right - left >= 12.0 for left, right in zip(xs, xs[1:]))

    assert capsys.readouterr().out == (
        "Plotted intervals: 2\n"
        f"Wrote: {output_path}\n"
    )
