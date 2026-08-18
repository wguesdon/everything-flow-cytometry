#!/usr/bin/env python3
"""Scan every FCS file in data/ and report which ones came from a spectral analyser.

A conventional analyser assigns one detector to one fluorochrome. A spectral
analyser records the whole emission spectrum across many detectors, so a raw
spectral file carries 60 or more channels and an unmixed one names the cytometer
in the ``$CYT`` keyword. Neither fact is visible from a folder name, which is why
this script reads the files.

Only the TEXT segment of each file is read, so no event is loaded and a 5 GB file
costs the same as a 5 MB one.

Example:
    uv run --project python python scripts/find_spectral_datasets.py
    uv run --project python python scripts/find_spectral_datasets.py --csv output/spectral_scan.csv
"""

from __future__ import annotations

import argparse
import collections
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fcs_header import read_text_segment

# Cytek Aurora and Northern Lights name a raw detector by laser and index, such as
# UV1-A, V7-A or B14-A. Sony ID7000 uses an FL<laser>-<detector> form.
CYTEK_DETECTOR = re.compile(r"^(UV|V|B|YG|R)\d{1,2}-[AHW]$")
SONY_DETECTOR = re.compile(r"^FL\d{1,2}-\d{1,2}")

# A cytometer name that identifies a spectral instrument.
SPECTRAL_CYT = ("aurora", "cytek", "northern lights", "id7000", "sp6800", "spectral")

# A file with at least this many parameters is worth a look even when the
# cytometer name is missing.
HIGH_PARAMETER_COUNT = 40

# At least this many detector names must match a raw spectral pattern.
MIN_SPECTRAL_DETECTORS = 20


def classify(keywords: dict[str, str]) -> dict[str, object]:
    """Describe one file from its FCS keywords.

    Args:
        keywords: The keyword dictionary from :func:`fcs_header.read_text_segment`.

    Returns:
        A dictionary with the cytometer name, the serial number, the parameter
        count, the number of named markers, the count of detector names that match
        each raw spectral pattern, and a ``flag`` naming the reason the file looks
        spectral. ``flag`` is an empty string when nothing matched.
    """
    cyt = keywords.get("$CYT", "").strip()
    try:
        par = int(keywords.get("$PAR", 0))
    except ValueError:
        par = 0

    cytek = sony = named = 0
    for i in range(1, par + 1):
        name = keywords.get(f"$P{i}N", "")
        if CYTEK_DETECTOR.match(name):
            cytek += 1
        if SONY_DETECTOR.match(name):
            sony += 1
        if keywords.get(f"$P{i}S", ""):
            named += 1

    if any(token in cyt.lower() for token in SPECTRAL_CYT):
        flag = "cyt_keyword"
    elif cytek >= MIN_SPECTRAL_DETECTORS:
        flag = "cytek_detectors"
    elif sony >= MIN_SPECTRAL_DETECTORS:
        flag = "sony_detectors"
    elif par >= HIGH_PARAMETER_COUNT:
        flag = "high_parameter_count"
    else:
        flag = ""

    return {
        "cyt": cyt,
        "cytsn": keywords.get("$CYTSN", ""),
        "par": par,
        "named_markers": named,
        "cytek_detectors": cytek,
        "sony_detectors": sony,
        "unmixed": "yes" if named > 0 and cytek == 0 else "no",
        "flag": flag,
    }


def scan(data_dir: Path) -> list[dict[str, object]]:
    """Read every FCS file under a folder and classify each one.

    Args:
        data_dir: The folder to walk.

    Returns:
        One record per readable file. A file that cannot be parsed is skipped and
        the reason is written to standard error.
    """
    files = sorted(data_dir.rglob("*.fcs"))
    print(f"Reading the TEXT segment of {len(files)} FCS files.", file=sys.stderr)

    records = []
    for i, path in enumerate(files, start=1):
        if i % 250 == 0:
            print(f"  {i} of {len(files)}", file=sys.stderr)
        try:
            keywords = read_text_segment(path)
        except (ValueError, OSError) as exc:
            print(f"  skipped {path}: {exc}", file=sys.stderr)
            continue

        relative = path.relative_to(data_dir)
        record = {"dataset": "/".join(relative.parts[:3]), "file": relative.name}
        record.update(classify(keywords))
        records.append(record)

    return records


def report(records: list[dict[str, object]]) -> None:
    """Print the spectral datasets and every distinct cytometer name.

    Args:
        records: The output of :func:`scan`.
    """
    flagged = [r for r in records if r["flag"]]
    print(f"\n{len(flagged)} of {len(records)} files look spectral.\n")

    grouped: dict[tuple, list] = collections.defaultdict(list)
    for r in flagged:
        grouped[(r["dataset"], r["cyt"], r["par"], r["unmixed"], r["flag"])].append(r)

    print(f"{'files':>6}  {'dataset':<58} {'cytometer':<30} {'par':>4} {'unmixed':<8} reason")
    for key, rows in sorted(grouped.items(), key=lambda kv: -len(kv[1])):
        dataset, cyt, par, unmixed, flag = key
        print(f"{len(rows):>6}  {dataset:<58} {cyt:<30} {par:>4} {unmixed:<8} {flag}")

    print("\nEvery distinct $CYT value in the archive:")
    for cyt, n in collections.Counter(r["cyt"] for r in records).most_common():
        print(f"  {n:>5}  {cyt or '(empty)'}")


def main() -> None:
    """Parse the arguments, scan the archive and print the report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data",
        help="Folder to scan. Defaults to data/ at the repository root.",
    )
    parser.add_argument("--csv", type=Path, help="Also write every record to this CSV.")
    args = parser.parse_args()

    if not args.data_dir.is_dir():
        raise SystemExit(
            f"Error: {args.data_dir} does not exist. Pull the archive with ./sync.sh pull."
        )

    records = scan(args.data_dir)
    if not records:
        raise SystemExit("Error: no FCS file was read, so nothing can be reported.")

    report(records)

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        print(f"\nWrote {args.csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
