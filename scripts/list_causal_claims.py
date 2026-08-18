#!/usr/bin/env python3
"""List every sentence in a report that asserts a cause.

`check_prose_rewrite.py` compares a rewritten report with its git revision and
fails on a changed code chunk, a changed block quote, an em dash, a bullet, an
invented number or an invented accession. It cannot fail on a sentence that
states the wrong cause, because a swapped mechanism changes no number and no
chunk.

That gap is not hypothetical. A rewrite of `reports/omip58_pytometry.qmd`
replaced "A rule that looks for a second mode either finds nothing ... or places
the cut inside the negative population" with a version naming the density rule
for the first outcome and the mixture rule for the second. The table printed
directly above that paragraph shows the opposite. The checker reported the file
as clean.

This script does not decide anything. It prints the sentences where that mistake
can hide, so a person reads each one against the table beside it. A sentence that
carries a number as well is the one to check first, because the number gives you
something to look up.

    uv run --project python python scripts/list_causal_claims.py reports/*.qmd
    uv run --project python python scripts/list_causal_claims.py --numbered-only \
        reports/omip58_pytometry.qmd
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

CAUSAL_PATTERN = re.compile(
    r"\b(because|since|therefore|so that|which is why|as a result|"
    r"rather than|instead of|means that|causes?|the reason|explains?|"
    r"due to|leads? to|drives?|so )\b",
    re.IGNORECASE,
)
DIGIT_PATTERN = re.compile(r"\d")


def prose_sentences(path: Path) -> list[str]:
    """Return the prose of a Quarto file, split into sentences.

    The YAML header, every code chunk, every header line, every block quote,
    every table row and every reference entry are dropped, because none of them
    is prose a rewrite was asked to change.

    Args:
        path: Path to the .qmd file.

    Returns:
        The sentences of the body text, in order.
    """
    lines = path.read_text().splitlines()
    body: list[str] = []
    in_chunk = False
    in_yaml = False

    for index, line in enumerate(lines):
        stripped = line.strip()
        if index == 0 and stripped == "---":
            in_yaml = True
            continue
        if in_yaml:
            if stripped == "---":
                in_yaml = False
            continue
        if stripped.startswith("```"):
            in_chunk = not in_chunk
            continue
        if in_chunk or not stripped:
            continue
        if stripped.startswith(("#", ">", "|", "[")):
            continue
        body.append(stripped)

    return re.split(r"(?<=[.]) +", " ".join(body))


def main() -> int:
    """Print the causal sentences of every file given.

    Returns:
        0 always. This script reports, it does not judge.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", type=Path, nargs="+", help="The .qmd files")
    parser.add_argument(
        "--numbered-only", action="store_true",
        help="Print only the sentences that also carry a number"
    )
    args = parser.parse_args()

    total = 0
    for path in args.paths:
        sentences = [
            sentence for sentence in prose_sentences(path)
            if CAUSAL_PATTERN.search(sentence)
            and (not args.numbered_only or DIGIT_PATTERN.search(sentence))
        ]
        total += len(sentences)
        print(f"\n### {path.name}: {len(sentences)} sentences")
        for number, sentence in enumerate(sentences, 1):
            print(f"{number:3}. {sentence}")

    print(f"\ntotal: {total}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
