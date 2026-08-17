#!/usr/bin/env python3
"""Check a prose rewrite of a Quarto report against the version in git.

The rewrite may change body text and nothing else. This script compares a
working tree file with a git revision of the same file and reports every
difference that the brief in `prompts/rewrite_report_as_prose.md` forbids.

Run it after a rewrite and before a commit:

    uv run --project python python scripts/check_prose_rewrite.py \
        reports/automated_gating_pbmc.qmd
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

# A number the prose may introduce without it being a claim. A section number, a
# figure number and a single digit inside a name are not measurements.
NUMBER_PATTERN = re.compile(r"\d[\d,]*\.?\d*")
DASH_PATTERN = re.compile(r"[–—]")


def read_revision(path: Path, revision: str) -> str:
    """Return the contents of a file at a git revision.

    Args:
        path: Repository relative path to the file.
        revision: Any git revision, for example `HEAD`.

    Returns:
        The file contents as text.

    Raises:
        SystemExit: If git cannot produce the file.
    """
    result = subprocess.run(
        ["git", "show", f"{revision}:{path}"],
        capture_output=True, text=True, check=False,
    )
    if result.returncode != 0:
        raise SystemExit(f"git could not read {revision}:{path}\n{result.stderr}")
    return result.stdout


def split_parts(text: str) -> dict[str, list[str]]:
    """Split a Quarto file into the parts a rewrite may and may not touch.

    Args:
        text: The whole file.

    Returns:
        A dictionary with the keys `yaml`, `chunks`, `quotes`, `headers`,
        `captions` and `body`.
    """
    lines = text.splitlines()
    yaml: list[str] = []
    chunks: list[str] = []
    quotes: list[str] = []
    headers: list[str] = []
    captions: list[str] = []
    body: list[str] = []

    in_yaml = False
    in_chunk = False
    current: list[str] = []

    for index, line in enumerate(lines):
        stripped = line.strip()
        if index == 0 and stripped == "---":
            in_yaml = True
            yaml.append(line)
            continue
        if in_yaml:
            yaml.append(line)
            if stripped == "---":
                in_yaml = False
            continue
        if stripped.startswith("```"):
            current.append(line)
            if in_chunk:
                chunks.append("\n".join(current))
                current = []
            in_chunk = not in_chunk
            continue
        if in_chunk:
            current.append(line)
            continue
        if stripped.startswith(">"):
            quotes.append(stripped)
            continue
        if stripped.startswith("#"):
            headers.append(stripped)
            continue
        if "caption =" in line or "fig-cap" in line:
            captions.append(stripped)
            continue
        if stripped:
            body.append(stripped)

    return {
        "yaml": yaml, "chunks": chunks, "quotes": quotes,
        "headers": headers, "captions": captions, "body": body,
    }


def numbers_in(lines: list[str]) -> set[str]:
    """Return every number that appears in a list of lines.

    Args:
        lines: The lines to scan.

    Returns:
        A set of the numbers as they were written.
    """
    found: set[str] = set()
    for line in lines:
        found.update(NUMBER_PATTERN.findall(line))
    return found


def main() -> int:
    """Compare one file with its git revision and report every violation.

    Returns:
        0 when the rewrite is clean, 1 otherwise.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path, help="The rewritten .qmd file")
    parser.add_argument("--revision", default="HEAD",
                        help="The git revision to compare against")
    args = parser.parse_args()

    new = split_parts(args.path.read_text())
    old = split_parts(read_revision(args.path, args.revision))

    problems: list[str] = []

    for part in ("yaml", "chunks", "quotes", "headers"):
        if new[part] != old[part]:
            problems.append(f"{part} changed")
            for line in set(old[part]) ^ set(new[part]):
                problems.append(f"    {line[:110]}")

    dashes = [line for line in new["body"] if DASH_PATTERN.search(line)]
    if dashes:
        problems.append(f"{len(dashes)} body lines carry an em dash or en dash")
        for line in dashes[:5]:
            problems.append(f"    {line[:110]}")

    bullets = [line for line in new["body"] if line.startswith("- ")]
    if bullets:
        problems.append(f"{len(bullets)} body lines are still bullets")
        for line in bullets[:5]:
            problems.append(f"    {line[:110]}")

    invented = numbers_in(new["body"]) - numbers_in(
        old["body"] + old["chunks"] + old["captions"] + old["headers"]
    )
    if invented:
        problems.append(f"numbers in the prose that are not in the original: "
                        f"{sorted(invented)}")

    print(f"{args.path}")
    print(f"  body lines: {len(old['body'])} before, {len(new['body'])} after")
    print(f"  code chunks: {len(new['chunks'])}, unchanged: "
          f"{new['chunks'] == old['chunks']}")
    if not problems:
        print("  clean")
        return 0
    for problem in problems:
        print(f"  {problem}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
