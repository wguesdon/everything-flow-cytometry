#!/usr/bin/env python3
"""Print the metadata keywords of an FCS file without loading the events.

An FCS file starts with a HEADER of fixed offsets, followed by a TEXT segment of
delimited key and value pairs. This script reads the TEXT segment only, so it is
fast on a large file and it needs no third party package.

Use it to see the panel, the instrument and the sample name before you decide
whether a dataset is worth loading.

Example:
    uv run python scripts/fcs_header.py data/datasets/flowjo/*/*/LD1_NS+NS_A01_exp.fcs
    uv run python scripts/fcs_header.py --panel data/datasets/**/*.fcs
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def read_text_segment(path: Path) -> dict[str, str]:
    """Read the TEXT segment of an FCS file and return its keywords.

    Args:
        path: Path to the FCS file.

    Returns:
        A dictionary of the FCS keywords. A key keeps the leading ``$`` that the
        standard uses for a required keyword.

    Raises:
        ValueError: If the file does not start with an FCS version string, or if
            the HEADER offsets are not readable.
    """
    with path.open("rb") as handle:
        header = handle.read(58)
        if not header[:3] == b"FCS":
            raise ValueError(f"{path} does not start with an FCS version string")

        try:
            start = int(header[10:18])
            end = int(header[18:26])
        except ValueError as exc:
            raise ValueError(f"{path} has an unreadable HEADER offset") from exc

        handle.seek(start)
        raw = handle.read(end - start + 1)

    text = raw.decode("latin-1")
    if not text:
        return {}

    delimiter = text[0]
    # A doubled delimiter escapes a literal delimiter character inside a value.
    parts = text[1:].split(delimiter)
    fields = [p for p in parts if p != ""]

    return dict(zip(fields[0::2], fields[1::2], strict=False))


def panel_rows(keywords: dict[str, str]) -> list[tuple[str, str, str]]:
    """Return the detector, the marker and the range for every parameter.

    Args:
        keywords: The keyword dictionary from :func:`read_text_segment`.

    Returns:
        One tuple per parameter, holding the ``$PnN`` detector name, the ``$PnS``
        marker name and the ``$PnR`` range. A missing value becomes an empty string.
    """
    count = int(keywords.get("$PAR", 0))
    rows = []
    for i in range(1, count + 1):
        rows.append(
            (
                keywords.get(f"$P{i}N", ""),
                keywords.get(f"$P{i}S", ""),
                keywords.get(f"$P{i}R", ""),
            )
        )
    return rows


def summarise(path: Path, show_panel: bool) -> None:
    """Print a summary of one FCS file.

    Args:
        path: Path to the FCS file.
        show_panel: Print the full parameter table when true.
    """
    try:
        keywords = read_text_segment(path)
    except (ValueError, OSError) as exc:
        print(f"{path}: {exc}", file=sys.stderr)
        return

    print(f"\n=== {path} ===")
    for key in ("$FIL", "$SRC", "$CYT", "$CYTSN", "$DATE", "$TOT", "$PAR", "$SPILLOVER"):
        value = keywords.get(key)
        if value is None:
            continue
        if key == "$SPILLOVER":
            value = f"present, {value.count(',') + 1} fields"
        print(f"{key:<12} {value}")

    if "GUID" in keywords:
        print(f"{'GUID':<12} {keywords['GUID']}")

    markers = [s for _, s, _ in panel_rows(keywords) if s]
    print(f"{'markers':<12} {len(markers)} named of {keywords.get('$PAR', '?')} parameters")

    if show_panel:
        print(f"\n{'detector':<24} {'marker':<24} range")
        for detector, marker, rng in panel_rows(keywords):
            print(f"{detector:<24} {marker:<24} {rng}")


def main() -> None:
    """Parse the arguments and summarise every file given."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", type=Path, help="FCS files to read")
    parser.add_argument(
        "--panel",
        action="store_true",
        help="print the detector and marker for every parameter",
    )
    args = parser.parse_args()

    for path in args.files:
        summarise(path, args.panel)


if __name__ == "__main__":
    main()
