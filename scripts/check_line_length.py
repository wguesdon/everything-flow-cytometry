"""Assert the line length that CLAUDE.md sets for each language.

The rule is 80 characters for R and 100 for Python. A rule nobody measures is a
rule nobody keeps, and a sweep in August 2026 found 225 R lines over the limit.

A line that cannot be broken without hurting it is exempt. That covers a URL, a
single token longer than the limit, a roxygen `@examples` line that has to stay
runnable, and a line that is one long string on a short piece of code. The last
one is the widest, so it is bounded: the code around the string has to be 20
characters or less, which leaves a call with several arguments and a long
message still reported. Every exemption is recognised from the line rather than listed
by name, so none of them can rot.

Run it from the repository root:

    uv run python scripts/check_line_length.py
    uv run python scripts/check_line_length.py --show
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import NamedTuple

REPO_ROOT = Path(__file__).resolve().parent.parent


class Rule(NamedTuple):
    """One language and the limit CLAUDE.md sets for it.

    Attributes:
        name: The language, for the report.
        limit: The longest line the rule allows.
        paths: Glob patterns, relative to the repository root.
    """

    name: str
    limit: int
    paths: tuple[str, ...]


RULES: tuple[Rule, ...] = (
    Rule("R", 80, ("R/*.R", "scripts/*.R", "scripts/cytokit/*.R", "tests/testthat/*.R",
                   "containers/*.R")),
    Rule("Python", 100, ("scripts/*.py", "python/tests/*.py")),
)

# A virtual environment and a cache hold other people's code, and its line
# length is not this repository's to keep.
SKIPPED_PARTS = frozenset({".venv", "venv", "__pycache__", "site-packages",
                           "node_modules", ".git"})


def is_exempt(line: str, in_examples: bool = False) -> bool:
    """Say whether a long line is one that cannot be broken.

    Args:
        line: The line, without its newline.
        in_examples: Whether the line sits in a roxygen `@examples` block.

    Returns:
        True when the line holds something that a break would damage.
    """
    stripped = line.strip()
    if "http://" in stripped or "https://" in stripped:
        return True
    # A roxygen example has to stay runnable as one expression, and the reader
    # copies it whole.
    if in_examples:
        return True
    # A long single token, such as a path or a checksum, has no break point.
    if max((len(word) for word in stripped.split()), default=0) > 60:
        return True
    # One long string on a short line of code. The code is a wrapper and the
    # string is the content, which covers a test description, a figure title, a
    # sprintf format and a file name from a deposit. Breaking one of those turns
    # a sentence into a concatenation, which reads worse and cannot be grepped.
    # The line still has to be almost all string: 20 characters of code is the
    # ceiling, so a call with several arguments and a long message is reported.
    without_strings = re.sub(r'"[^"]*"|\'[^\']*\'', '""', stripped)
    return len(without_strings) <= 20 and stripped.count('"') >= 2


def main() -> int:
    """Run the check.

    Returns:
        0 when every line is inside its limit, 1 otherwise.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--show", action="store_true", help="print every line that is over")
    options = parser.parse_args()

    status = 0
    for rule in RULES:
        offenders: list[tuple[Path, int, int]] = []
        files = 0
        for pattern in rule.paths:
            for path in sorted(REPO_ROOT.glob(pattern)):
                if SKIPPED_PARTS & set(path.parts):
                    continue
                files += 1
                # A script may hold a byte that is not UTF-8, and a length check
                # is not the place to stop for it.
                text = path.read_text(encoding="utf-8", errors="replace")
                in_examples = False
                for number, line in enumerate(text.splitlines(), start=1):
                    stripped = line.lstrip()
                    if stripped.startswith("#'") and "@" in stripped:
                        in_examples = "@examples" in stripped
                    elif not stripped.startswith("#'"):
                        in_examples = False
                    if len(line) > rule.limit and not is_exempt(line, in_examples):
                        offenders.append((path.relative_to(REPO_ROOT), number, len(line)))

        print(f"{rule.name}: {files} file(s), limit {rule.limit}, "
              f"{len(offenders)} line(s) over")
        if offenders:
            status = 1
            shown = offenders if options.show else offenders[:5]
            for path, number, length in shown:
                print(f"  {path}:{number} is {length}")
            if not options.show and len(offenders) > len(shown):
                print(f"  and {len(offenders) - len(shown)} more. Pass --show for all.")

    if status == 0:
        print()
        print("Every line is inside the limit its language sets.")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
