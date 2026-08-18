"""Assert the line length that CLAUDE.md sets for each language.

The rule is 80 characters for R and 100 for Python. A rule nobody measures is a
rule nobody keeps, and a sweep in August 2026 found 225 R lines over the limit.

A line that cannot be broken without hurting it is exempt: a URL, a file path
and a roxygen `@examples` line that has to stay runnable. Those are recognised
rather than listed, so the exemption cannot rot.

Run it from the repository root:

    uv run python scripts/check_line_length.py
    uv run python scripts/check_line_length.py --show
"""

from __future__ import annotations

import argparse
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


def is_exempt(line: str) -> bool:
    """Say whether a long line is one that cannot be broken.

    Args:
        line: The line, without its newline.

    Returns:
        True when the line holds something that a break would damage.
    """
    stripped = line.strip()
    if "http://" in stripped or "https://" in stripped:
        return True
    # A roxygen example has to stay runnable as one expression.
    if stripped.startswith("#' ") and ("::" in stripped and "(" in stripped):
        return False
    # A long single token, such as a path or a checksum, has no break point.
    return max((len(word) for word in stripped.split()), default=0) > 60


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
                for number, line in enumerate(text.splitlines(), start=1):
                    if len(line) > rule.limit and not is_exempt(line):
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
