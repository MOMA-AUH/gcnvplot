"""CLI tests for gatk-germline-cnv-plotter."""

from __future__ import annotations

import pytest

import gatk_germline_cnv_plotter
import gatk_germline_cnv_plotter.cli as cli


def test_main_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--version"])

    assert exc_info.value.code == 0
    assert (
        capsys.readouterr().out
        == f"gatk-germline-cnv-plotter {gatk_germline_cnv_plotter.__version__}\n"
    )


def test_main_requires_subcommand(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main([])

    assert exc_info.value.code == 2
    assert "required" in capsys.readouterr().err


def test_plot_command_reports_placeholder_output(
    tmp_path, capsys: pytest.CaptureFixture[str]
) -> None:
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "plot.png"

    assert (
        cli.main(
            [
                "plot",
                "--input",
                str(input_path),
                "--output",
                str(output_path),
            ]
        )
        == 0
    )

    assert capsys.readouterr().out == (
        f"Input: {input_path}\n"
        f"Output: {output_path}\n"
        "Plot generation is not implemented yet.\n"
    )
