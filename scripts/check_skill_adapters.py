#!/usr/bin/env python3
"""Check that the skill adapters agree with the cytokit CLI.

Three adapters teach three host tools how to call `cytokit`. Each one lists the
commands, and each one is a separate file, so they drift. An adapter that sends
an agent to a command the CLI does not have is worse than no adapter, because
the agent reports a failure that looks like the scientist's fault.

This compares three things: the subcommands `cli/cytokit` actually dispatches,
the recipes `cytokit list` calls ready, and the commands each adapter documents.

Run it before a commit that touches the CLI or an adapter:

    python scripts/check_skill_adapters.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

CLI = Path("cli/cytokit")
ADAPTERS = [
    Path("skills/claude-code/SKILL.md"),
    Path("skills/codex/AGENTS.md"),
    Path("skills/opencode/AGENTS.md"),
]

# The line in cli/cytokit that dispatches a recipe to its R script.
DISPATCH_PATTERN = re.compile(r"^\s{4}([a-z|]+)\)\s*$", re.MULTILINE)

# A row of the table that `cytokit list` prints.
LIST_PATTERN = re.compile(r"^([a-z]+)\s+(ready|planned)\s", re.MULTILINE)

# A command line inside an adapter's fenced block.
ADAPTER_PATTERN = re.compile(r"^cytokit\s+([a-z]+)", re.MULTILINE)

# Commands that are not recipes and need no row in the list table.
UTILITY = {"build", "list", "shell", "version", "help"}


def dispatched() -> set[str]:
    """Return the recipes that the CLI routes to an R script.

    Returns:
        The subcommand names in the dispatch case statement.
    """
    text = CLI.read_text()
    block = re.search(r"^\s{4}(inspect[a-z|]*)\)", text, re.MULTILINE)
    return set(block.group(1).split("|")) if block else set()


def listed() -> dict[str, str]:
    """Return the recipe table that `cytokit list` prints.

    Returns:
        A mapping of recipe name to `ready` or `planned`.
    """
    return dict(LIST_PATTERN.findall(CLI.read_text()))


def documented(path: Path) -> set[str]:
    """Return the commands one adapter tells an agent to call.

    Args:
        path: The adapter file.

    Returns:
        The command names it names.
    """
    return set(ADAPTER_PATTERN.findall(path.read_text()))


def main() -> int:
    """Compare the CLI with every adapter.

    Returns:
        0 when they agree, 1 when they do not.
    """
    routes = dispatched()
    table = listed()
    ready = {name for name, state in table.items() if state == "ready"}
    problems: list[str] = []

    if routes != ready:
        problems.append(
            f"cli/cytokit routes {sorted(routes)} but calls {sorted(ready)} ready"
        )

    for adapter in ADAPTERS:
        if not adapter.exists():
            problems.append(f"{adapter} is missing")
            continue
        names = documented(adapter) - UTILITY
        unknown = names - set(table)
        stale = names - ready - UTILITY
        if unknown:
            problems.append(
                f"{adapter} names a command the CLI does not have: "
                f"{sorted(unknown)}"
            )
        elif stale:
            problems.append(
                f"{adapter} tells an agent to call a planned recipe: "
                f"{sorted(stale)}"
            )
        if ready - names:
            problems.append(
                f"{adapter} does not document a ready recipe: "
                f"{sorted(ready - names)}"
            )

    print(f"routes   {sorted(routes)}")
    print(f"ready    {sorted(ready)}")
    print(f"planned  {sorted(n for n, s in table.items() if s == 'planned')}")
    for adapter in ADAPTERS:
        if adapter.exists():
            print(f"{adapter.name:<14} {sorted(documented(adapter) - UTILITY)}")

    if problems:
        print("\nThe adapters and the CLI disagree:", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        return 1
    print("\nThe CLI and all three adapters agree.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
