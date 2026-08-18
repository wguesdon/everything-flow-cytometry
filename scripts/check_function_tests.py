"""Assert that every function in R/ has a test, or a written reason for none.

CLAUDE.md says every function in `R/` has a testthat test in `tests/testthat/`.
A sweep in August 2026 found 56 that had none, and nothing in the repository
said so. A rule nobody measures is a rule nobody keeps.

A handful of functions read a real FCS file or a real FlowJo workspace, and the
repository does not put a 17 MB file in its test suite. Those are listed in
`COVERED_BY_A_SCRIPT` with the script that exercises each one. Adding a name to
that list is a decision a reader can see and argue with, which an empty gap is
not.

Run it from the repository root:

    uv run python scripts/check_function_tests.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SOURCE_DIRECTORY = REPO_ROOT / "R"
TEST_DIRECTORY = REPO_ROOT / "tests" / "testthat"
DEFINITION = re.compile(r"^([A-Za-z][A-Za-z0-9._]*)\s*<-\s*function")

# Each of these reads a real FCS file, a real FlowJo workspace or a whole
# deposit, so a unit test would need a file this repository does not put in its
# test suite. The script named beside each one runs it on the real data.
COVERED_BY_A_SCRIPT: dict[str, str] = {
    "GateBoneMarrowFile": "scripts/10_oetjen2018_bone_marrow.R",
    "GateTcellSubsets": "scripts/10_oetjen2018_bone_marrow.R",
    "GateMassCytometryFile": "scripts/10_oetjen2018_bone_marrow.R",
    "ImportFlowJoGates": "scripts/04_omip39_reproduce_paper.R",
    "GateCovidFile": "scripts/11_vanderbeke2021_covid.R",
    "GateHvtnSubsets": "scripts/12_flowcap2_challenges.R",
    "GateHeuSubsets": "scripts/12_flowcap2_challenges.R",
    "GateNaiveMemoryFile": "scripts/09_z282_harmonisation.R",
    "ComputeOmip58Spillover": "scripts/13_omip58_prepare.R",
    "EstimateSpectralTransform": "scripts/08_yu2021_spectral_mait.R",
    "GateSpectralFile": "scripts/08_yu2021_spectral_mait.R",
    "TestYuClaims": "scripts/08_yu2021_spectral_mait.R",
}


def functions_in(path: Path) -> list[str]:
    """List the functions a source file defines at the top level.

    Args:
        path: The R source file.

    Returns:
        The function names, in the order they are defined.
    """
    names = []
    for line in path.read_text().splitlines():
        match = DEFINITION.match(line)
        if match:
            names.append(match.group(1))
    return names


def main() -> int:
    """Run the check.

    Returns:
        0 when every function is tested or listed, 1 otherwise.
    """
    tests = "\n".join(path.read_text() for path in sorted(TEST_DIRECTORY.glob("*.R")))

    untested: list[tuple[str, str]] = []
    listed: list[tuple[str, str]] = []
    tested = 0
    for source in sorted(SOURCE_DIRECTORY.glob("*.R")):
        for name in functions_in(source):
            if f"{name}(" in tests:
                tested += 1
            elif name in COVERED_BY_A_SCRIPT:
                listed.append((source.name, name))
            else:
                untested.append((source.name, name))

    total = tested + len(listed) + len(untested)
    print(f"{total} functions in R/")
    print(f"  {tested} have a test that names them")
    print(f"  {len(listed)} are exercised by a script, and say which one")
    print(f"  {len(untested)} have neither")

    # A name that is listed but has since gained a test, or been deleted, makes
    # the list a record of what used to be true.
    stale = sorted(set(COVERED_BY_A_SCRIPT) - {name for _, name in listed})
    if stale:
        print()
        print("These names are listed as covered by a script and no longer need to be.")
        print("Either they gained a test or they were deleted. Take them out of the list:")
        for name in stale:
            print(f"  {name}")

    if untested:
        print()
        print("These functions have no test and no reason recorded:", file=sys.stderr)
        for source, name in untested:
            print(f"  {source}: {name}", file=sys.stderr)
        print(file=sys.stderr)
        print("Write a test, or add the name to COVERED_BY_A_SCRIPT with the", file=sys.stderr)
        print("script that runs it.", file=sys.stderr)
        return 1
    if stale:
        return 1
    print()
    print("Every function in R/ has a test or a written reason for none.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
