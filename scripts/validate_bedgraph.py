#!/usr/bin/env python3
"""Validate non-empty bedGraph output, including scientific-notation values."""

import math
import sys


def validate(path):
    rows = 0
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("track", "browser", "#")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise ValueError(f"{path}:{line_number}: expected four tab-separated fields")
            try:
                start, end = int(fields[1]), int(fields[2])
                value = float(fields[3])
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid coordinate or value") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval {start}-{end}")
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{path}:{line_number}: coverage must be finite and non-negative")
            rows += 1
    if rows == 0:
        raise ValueError(f"{path}: contains no data rows")


if __name__ == "__main__":
    try:
        validate(sys.argv[1])
    except (IndexError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
