"""Read-count file parsing for gcnvplot."""

from __future__ import annotations

import csv
from pathlib import Path

from .models import Interval
from .utils import open_text

READ_COUNT_FIELDS = ["CONTIG", "START", "END", "COUNT"]


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
