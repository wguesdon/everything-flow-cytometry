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
DOCUMENT = Path("docs/cytokit.md")
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

    # docs/cytokit.md carries the same recipe table, and a state that lags the
    # CLI sends a reader to a recipe they are told does not exist.
    document = DOCUMENT
    if not document.exists():
        problems.append(f"{document} is missing")
    else:
        stated = document_states(document)
        if not stated:
            problems.append(f"{document} carries no recipe table")
        for name, state in sorted(table.items()):
            if name not in stated:
                problems.append(f"{document} does not list the recipe {name}")
            elif stated[name] != state:
                problems.append(
                    f"{document} calls {name} {stated[name]} and the CLI calls it {state}"
                )
        # A row in a table is a name, not a description. A ready recipe needs a
        # section that says what it does and what its flags mean.
        headings = {
            line[3:].strip().lower()
            for line in document.read_text().splitlines()
            if line.startswith("## ")
        }
        for name in sorted(ready):
            if name not in headings:
                problems.append(f"{document} has no section for the ready recipe {name}")

    # A flag the help names and the recipe rejects strands an agent on an
    # error that reads like its own fault.
    for name, flags in sorted(flags_in_usage().items()):
        recipe = Path("scripts/cytokit") / f"{name}.R"
        if not recipe.exists():
            problems.append(f"the help names {name} and scripts/cytokit/{name}.R is missing")
            continue
        unknown = flags - accepted_flags(recipe)
        if unknown:
            problems.append(
                f"the help gives {name} the flag(s) {sorted(unknown)} which it does not accept"
            )

    # README.md tells a stranger what is built. A sentence that lags the CLI is
    # the first thing they read and the last thing anybody edits.
    readme = Path("README.md")
    if readme.exists():
        text = readme.read_text()
        for name, state in sorted(table.items()):
            if state == "ready" and re.search(rf"`{name}`[^.\n]*\bplanned\b", text):
                problems.append(f"README.md calls {name} planned and the CLI calls it ready")

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


def flags_in_usage() -> dict[str, set[str]]:
    """Read the flags the CLI help gives to each recipe.

    Returns:
        A mapping from recipe name to the flags its usage block names.
    """
    text = CLI.read_text()
    usage = text[text.index("Usage:"):text.index("  cytokit list")]
    claimed: dict[str, set[str]] = {}
    recipe = None
    for line in usage.splitlines():
        match = re.match(r"^  cytokit (\w+)\s", line)
        if match:
            recipe = match.group(1)
            claimed.setdefault(recipe, set())
        if recipe:
            claimed[recipe].update(re.findall(r"--([a-z-]+)", line))
    return claimed


def accepted_flags(path: Path) -> set[str]:
    """Read the arguments and flags a recipe accepts.

    Args:
        path: The recipe under scripts/cytokit/.

    Returns:
        Every name its ParseCytokitArguments call allows.
    """
    source = path.read_text()
    block = source[source.index("ParseCytokitArguments"):]
    return set(re.findall(r'"([a-z-]+)"', block[:block.index(")\n")]))


def document_states(path: Path) -> dict[str, str]:
    """Read the recipe table of docs/cytokit.md.

    Args:
        path: The document.

    Returns:
        A mapping from recipe name to the state the document gives it.
    """
    states: dict[str, str] = {}
    for line in path.read_text().splitlines():
        match = re.match(r"^\|\s*`([a-z]+)`\s*\|\s*(ready|planned)\s*\|", line)
        if match:
            states[match.group(1)] = match.group(2)
    return states


if __name__ == "__main__":
    raise SystemExit(main())
