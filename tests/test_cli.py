"""CLI tests for gcnvplot."""

from __future__ import annotations

import pytest

import gcnvplot
import gcnvplot.cli as cli


def test_main_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc_info:
        cli.main(["--version"])

    assert exc_info.value.code == 0
    assert (
        capsys.readouterr().out
        == f"gcnvplot {gcnvplot.__version__}\n"
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
