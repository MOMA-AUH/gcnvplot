"""CLI tests for gcnvplot."""

from __future__ import annotations

from pathlib import Path

import pytest

import gcnvplot
import gcnvplot.cli as cli


def write_read_counts(path: Path, rows: list[tuple[str, int, int, int]]) -> None:
    """Write a tiny CollectReadCounts-style TSV."""
    lines = ["@RG\tID:test", "CONTIG\tSTART\tEND\tCOUNT"]
    lines.extend(f"{contig}\t{start}\t{end}\t{count}" for contig, start, end, count in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
    assert "# normalization=median" in background_text
    assert "chr1\t100\t199\t2\t0.66666667\t0.66666667" in background_text
    assert "chr1\t200\t299\t2\t1.3333333\t1.3333333" in background_text

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
    assert "the following arguments are required: --region" in capsys.readouterr().err


def test_plot_log2_ratio_writes_svg_with_zero_centered_axis(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
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
                "chr2\t100\t199\t2\t1\t1\t0.1\t0.8\t1.2",
                "chr1\t200\t299\t2\t1\t1\t0.1\t0.8\t1.2",
                "chr2\t200\t299\t2\t1\t1\t0.1\t0.8\t1.2",
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
                "--output",
                str(output_path),
            ]
        )
        == 0
    )

    svg = output_path.read_text(encoding="utf-8")
    assert "chr2:100-199" not in svg
    assert "chr2:200-299" not in svg
    assert "Log2(sample/background)" in svg
    assert "SIGNAL=-1.3006595" in svg
    assert "SIGNAL=0.6727054" in svg

    assert capsys.readouterr().out == (
        "Plotted intervals: 2\n"
        f"Wrote: {output_path}\n"
    )
