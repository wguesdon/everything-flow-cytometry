#!/usr/bin/env python3
"""Check that every number in a report's prose still appears in its output.

A report states its results as literal numbers in running text, while its tables
are read from CSV files when the report renders. A table therefore cannot go
stale and a sentence can. Re-running an analysis is the moment a sentence goes
stale, so this script reads the numbers out of the prose and looks for each one
in the CSV files that the same analysis wrote.

A number that is not found is a candidate and not a verdict. A report also
carries numbers that no output holds: a year, a PMID, a panel size, a value
quoted from the paper under test. Read the list and judge each entry.

Run it after an analysis re-run and before a commit:

    python scripts/check_report_numbers.py reports/omip51_bcell_dc.qmd
    python scripts/check_report_numbers.py reports/*.qmd
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# The same number form that scripts/check_prose_rewrite.py uses, with one
# addition. The lookbehind drops a number that is part of a name: the 141 of
# CD141, the 780 of B780-A and the 051 of OMIP-051 are labels and not results.
NUMBER_PATTERN = re.compile(
    r"(?<![A-Za-z0-9.-])-?\d+(?:,\d{3})*(?:\.\d+)?")

# Anything a report sets in a code span is a channel name, a file name or an
# argument, so its digits are not results either.
CODE_SPAN_PATTERN = re.compile(r"`[^`]*`")

# A four digit number in this range is a year.
YEAR_RANGE = range(1900, 2101)

# Everything from this heading onward is a citation list, whose volume numbers,
# page ranges, PMIDs and DOIs are not results.
REFERENCES_HEADING = re.compile(r"^#+\s+References\s*$", re.IGNORECASE)

# The line that names the folder an analysis wrote.
OUTPUT_DIR_PATTERN = re.compile(
    r'kOutputDir\s*<-\s*file\.path\((.*?)\)', re.DOTALL)

# A number this small carries no evidence, and every table holds it somewhere.
SMALLEST_CHECKED = 3


def prose_lines(text: str) -> list[tuple[int, str]]:
    """Return the numbered lines of a Quarto file that a reader reads.

    Code chunks, YAML front matter and chunk options are dropped, because a
    number inside them is an argument and not a claim.

    Args:
        text: The whole contents of the .qmd file.

    Returns:
        A list of (line number, line) pairs, counting from 1.
    """
    kept: list[tuple[int, str]] = []
    in_chunk = False
    in_yaml = False
    for number, line in enumerate(text.split("\n"), start=1):
        stripped = line.strip()
        if number == 1 and stripped == "---":
            in_yaml = True
            continue
        if in_yaml:
            if stripped == "---":
                in_yaml = False
            continue
        if stripped.startswith("```"):
            in_chunk = not in_chunk
            continue
        if in_chunk or stripped.startswith("#|"):
            continue
        if REFERENCES_HEADING.match(stripped):
            break
        kept.append((number, line))
    return kept


def output_directory(report: Path) -> Path | None:
    """Find the output folder that a report reads.

    Args:
        report: The path to the .qmd file.

    Returns:
        The folder, or None when the report does not name one.
    """
    match = OUTPUT_DIR_PATTERN.search(report.read_text())
    if match is None:
        return None
    parts = re.findall(r'"([^"]+)"', match.group(1))
    if not parts:
        return None
    return (report.parent / Path(*parts)).resolve()


def numbers_in_outputs(directory: Path) -> set[str]:
    """Collect every number that the CSV files in a folder hold.

    A value is recorded in several forms, because a report rounds. 0.5825 is
    written 0.58 in one sentence and 0.583 in another, and both must match.

    Args:
        directory: The folder of CSV files.

    Returns:
        The set of number strings that the folder supports.
    """
    found: set[str] = set()
    for path in sorted(directory.rglob("*.csv")):
        try:
            with path.open(newline="") as handle:
                for row in csv.reader(handle):
                    for cell in row:
                        text = cell.strip()
                        forms = _forms(text)
                        if forms:
                            found.update(forms)
                            continue
                        # A verdict column carries its evidence as a sentence,
                        # so the numbers inside it support the prose as well.
                        found.update(NUMBER_PATTERN.findall(text))
        except (OSError, csv.Error):
            continue
    return found


def _forms(cell: str) -> set[str]:
    """Return the written forms of one cell value.

    Args:
        cell: The raw cell text.

    Returns:
        Every string a report might use for this value.
    """
    try:
        value = float(cell)
    except ValueError:
        return set()
    forms = {cell}
    if value.is_integer():
        whole = int(value)
        forms.add(str(whole))
        forms.add(f"{whole:,}")
        # A percentage stated as a whole number of events, and the other way.
        forms.add(f"{whole}.0")
    for places in (0, 1, 2, 3, 4):
        rounded = round(value, places)
        forms.add(f"{rounded:.{places}f}")
        if rounded.is_integer():
            forms.add(f"{int(rounded):,}")
    # A fraction written as a percentage, and a percentage written as a fraction.
    for scaled in (value * 100, value / 100):
        for places in (1, 2, 3):
            forms.add(f"{round(scaled, places):.{places}f}")
    return forms


def check(report: Path) -> list[tuple[int, str, str]]:
    """Report every prose number that the analysis output does not hold.

    Args:
        report: The path to the .qmd file.

    Returns:
        A list of (line number, number, line) for each unmatched number.

    Raises:
        FileNotFoundError: When the report names an output folder that is not
            on disk, which means the analysis has not been run.
    """
    directory = output_directory(report)
    if directory is None:
        return []
    if not directory.is_dir():
        raise FileNotFoundError(f"{report}: {directory} does not exist")
    supported = numbers_in_outputs(directory)

    missing: list[tuple[int, str, str]] = []
    for number, line in prose_lines(report.read_text()):
        for match in NUMBER_PATTERN.finditer(CODE_SPAN_PATTERN.sub(" ", line)):
            written = match.group(0)
            plain = written.replace(",", "")
            try:
                value = float(plain)
            except ValueError:
                continue
            if abs(value) < SMALLEST_CHECKED:
                continue
            if value.is_integer() and int(value) in YEAR_RANGE and "," not in written:
                continue
            if written in supported or plain in supported:
                continue
            missing.append((number, written, line.strip()))
    return missing


def main() -> int:
    """Check each report named on the command line.

    Returns:
        0 when every report was read, 1 when one could not be checked.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+", type=Path)
    arguments = parser.parse_args()

    failed = False
    for report in arguments.reports:
        try:
            missing = check(report)
        except FileNotFoundError as error:
            print(f"SKIP {error}", file=sys.stderr)
            failed = True
            continue
        print(f"\n=== {report} : {len(missing)} unmatched ===")
        for line_number, written, line in missing:
            print(f"  line {line_number:>4}  {written:>12}  {line[:96]}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
