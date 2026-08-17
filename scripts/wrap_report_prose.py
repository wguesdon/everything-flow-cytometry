#!/usr/bin/env python3
"""Rewrap the prose of a Quarto report at a fixed column.

The rest of the repository wraps at 80 characters. A model asked to rewrap a
paragraph can change a word without meaning to, so this does it mechanically.

Only prose is touched. The YAML header, every code chunk, every header line,
every block quote, every table row and every reference entry are copied through
unchanged, and no word is altered.

    uv run --project python python scripts/wrap_report_prose.py \
        reports/flowcap2_challenges.qmd
"""

from __future__ import annotations

import argparse
import re
import sys
import textwrap
from pathlib import Path

REFERENCE_PATTERN = re.compile(r"^\[\d+\]")


def is_prose(line: str) -> bool:
    """Say whether a line is ordinary prose that may be rewrapped.

    Args:
        line: One line of the file.

    Returns:
        True when the line carries prose and nothing structural.
    """
    stripped = line.strip()
    if not stripped:
        return False
    if stripped.startswith(("#", ">", "|", "!", "-", "*", ":")):
        return False
    if REFERENCE_PATTERN.match(stripped):
        return False
    return not stripped.startswith("```")


def rewrap(text: str, width: int) -> tuple[str, int]:
    """Rewrap every prose paragraph of a Quarto file.

    Args:
        text: The whole file.
        width: The column to wrap at.

    Returns:
        The rewrapped file and the number of paragraphs that changed.
    """
    lines = text.splitlines()
    out: list[str] = []
    paragraph: list[str] = []
    changed = 0

    in_yaml = False
    in_chunk = False

    def flush() -> None:
        nonlocal changed
        if not paragraph:
            return
        joined = " ".join(part.strip() for part in paragraph)
        wrapped = textwrap.wrap(joined, width=width, break_long_words=False,
                                break_on_hyphens=False)
        if wrapped != paragraph:
            changed += 1
        out.extend(wrapped)
        paragraph.clear()

    for index, line in enumerate(lines):
        stripped = line.strip()
        if index == 0 and stripped == "---":
            in_yaml = True
            out.append(line)
            continue
        if in_yaml:
            out.append(line)
            if stripped == "---":
                in_yaml = False
            continue
        if stripped.startswith("```"):
            flush()
            in_chunk = not in_chunk
            out.append(line)
            continue
        if in_chunk:
            out.append(line)
            continue
        if is_prose(line):
            paragraph.append(line)
            continue
        flush()
        out.append(line)

    flush()
    return "\n".join(out) + "\n", changed


def main() -> int:
    """Rewrap one or more files and report what changed.

    Returns:
        0 always, because a file with nothing to rewrap is not an error.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", type=Path, nargs="+", help="The .qmd files")
    parser.add_argument("--width", type=int, default=80,
                        help="The column to wrap at")
    args = parser.parse_args()

    for path in args.paths:
        before = path.read_text()
        after, changed = rewrap(before, args.width)
        words_before = before.split()
        words_after = after.split()
        if words_before != words_after:
            print(f"{path}: REFUSED, the word sequence would change")
            continue
        if after != before:
            path.write_text(after)
        over = sum(
            1 for line in after.splitlines()
            if is_prose(line) and len(line) > args.width
        )
        print(f"{path}: {changed} paragraphs rewrapped, "
              f"{over} prose lines still over {args.width}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
