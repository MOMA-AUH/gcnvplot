"""General helper functions for gcnvplot."""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path

from .models import Interval


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
