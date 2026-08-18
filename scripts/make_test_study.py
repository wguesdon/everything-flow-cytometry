"""Build a folder that looks like a scientist's own data, outside this repository.

The corpus check proves that the recipes run. It does not prove that an agent
reading the skill chooses the right one, and no script can prove that. The test
is a session that starts with nothing but the adapter and a folder of FCS files.

This script writes that folder. It copies files out of the archive so that the
session has no repository context to lean on, and by default it writes no
metadata table, because an agent that invents a treatment assignment rather than
asking for one has failed the test.

Run it from the repository root:

    uv run python scripts/make_test_study.py --deposit FR-FCM-ZZCA --to ~/cytokit_test
    uv run python scripts/make_test_study.py --list
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import NamedTuple

REPO_ROOT = Path(__file__).resolve().parent.parent
DEPOSIT_ROOT = REPO_ROOT / "data" / "datasets" / "flowrepository"


class Study(NamedTuple):
    """One deposit that makes a useful test study.

    Attributes:
        name: The name to pass to `--deposit`.
        path: The folder holding the FCS files, under the archive root.
        why: What an agent has to get right on this study.
    """

    name: str
    path: str
    why: str


STUDIES: tuple[Study, ...] = (
    Study("FR-FCM-ZZCA", "FR-FCM-ZZCA/FlowRepository_FR-FCM-ZZCA_files",
          "five small files, a stimulated and a control arm in the file names, "
          "no marker names at all, and an identity spillover matrix"),
    Study("FR-FCM-ZZLV", "FlowRepository_FR-FCM-ZZLV_files",
          "three files, every marker named, a plain T cell panel"),
    Study("FR-FCM-Z244", "FlowRepository_FR-FCM-Z244_files",
          "mass cytometry, so a scatter gate and a compensation step are both wrong"),
    Study("FR-FCM-Z4KT", "FlowRepository_FR-FCM-Z4KT_files",
          "a spectral Aurora panel of 39 detectors"),
)


def copy_study(study: Study, destination: Path, limit: int | None) -> list[Path]:
    """Copy the FCS files of one deposit into a fresh folder.

    Args:
        study: The deposit to copy.
        destination: The folder to write. It must not already hold FCS files.
        limit: How many files to copy, or None for all of them.

    Returns:
        The files that were written.

    Raises:
        SystemExit: When the source is absent or the destination is not empty.
    """
    source = DEPOSIT_ROOT / study.path
    if not source.is_dir():
        raise SystemExit(f"The deposit is not here: {source}")

    existing = sorted(destination.glob("*.fcs")) if destination.exists() else []
    if existing:
        raise SystemExit(
            f"{destination} already holds {len(existing)} FCS file(s). "
            "Name an empty folder, so that a run starts from a known state."
        )

    files = sorted(source.glob("*.fcs")) + sorted(source.glob("*.FCS"))
    if not files:
        raise SystemExit(f"No FCS file is in {source}")
    if limit is not None:
        files = files[:limit]

    destination.mkdir(parents=True, exist_ok=True)
    written = []
    for path in files:
        target = destination / path.name
        shutil.copy2(path, target)
        written.append(target)
    return written


def main() -> int:
    """Build the test study.

    Returns:
        The process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deposit", default="FR-FCM-ZZCA", help="which study to copy")
    parser.add_argument("--to", type=Path, help="the folder to write, outside this repository")
    parser.add_argument("--limit", type=int, help="copy only the first N files")
    parser.add_argument("--list", action="store_true", help="show the studies and stop")
    options = parser.parse_args()

    if options.list:
        for study in STUDIES:
            print(f"{study.name:<14} {study.why}")
        return 0

    if options.to is None:
        parser.error("--to is required, and it has to be outside this repository")

    destination = options.to.expanduser().resolve()
    if destination == REPO_ROOT or REPO_ROOT in destination.parents:
        raise SystemExit(
            f"{destination} is inside this repository. A test study belongs outside it, "
            "so that the session has no repository context to read."
        )

    matches = [study for study in STUDIES if study.name == options.deposit]
    if not matches:
        raise SystemExit(
            f"No study is called {options.deposit}. "
            f"The studies are: {', '.join(study.name for study in STUDIES)}"
        )

    study = matches[0]
    written = copy_study(study, destination, options.limit)
    print(f"Wrote {len(written)} FCS file(s) to {destination}")
    print(f"This study tests: {study.why}")
    print()
    print("No metadata table was written. An agent that assigns a treatment without")
    print("asking has failed the test.")
    print()
    print("Start a session in that folder and ask for the analysis in your own words.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
