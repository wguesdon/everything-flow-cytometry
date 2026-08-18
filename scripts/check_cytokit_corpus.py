"""Run every ready cytokit recipe over FCS deposits that no analysis here uses.

The unit tests in `tests/testthat/test-cytokit.R` check the functions against a
`flowFrame` built in the test. That proves the logic and it proves nothing about
the panels a scientist arrives with. A skill fails on the first deposit whose
convention its author never saw, so this script runs the ready recipes over the
deposits in `data/` that no script and no report reads.

Each deposit in the corpus is there because something about it differs: an empty
`$PnS`, an identity spillover matrix, a spectral instrument, a file per tube
rather than a file per sample.

Run it from the repository root:

    uv run python scripts/check_cytokit_corpus.py
    uv run python scripts/check_cytokit_corpus.py --deposit FR-FCM-ZZCA
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import NamedTuple

REPO_ROOT = Path(__file__).resolve().parent.parent
CLI = REPO_ROOT / "cli" / "cytokit"
DEPOSIT_ROOT = REPO_ROOT / "data" / "datasets" / "flowrepository"

# Every bundle carries these, whatever the recipe wrote on top.
BUNDLE_FILES = ("manifest.json", "session_info.txt", "REPRODUCE.md")
INSPECT_FILES = ("files.csv", "panel.csv", "markers.txt", "naming.csv")


class Deposit(NamedTuple):
    """One deposit in the corpus.

    Attributes:
        name: A short name for the report line.
        path: The folder under `data/datasets/flowrepository/` to pass as `--data`.
        why: What this deposit tests that the others do not.
        recursive: Whether the FCS files sit below `path` rather than in it.
    """

    name: str
    path: str
    why: str
    recursive: bool = False


# No script and no report in this repository reads any of these. That is the
# point: a recipe that only works on the deposits it was written against is not
# a tool a scientist can use on their own data.
CORPUS: tuple[Deposit, ...] = (
    Deposit("FR-FCM-ZZCA", "FR-FCM-ZZCA",
            "no marker names at all, an identity spillover matrix, and the FCS files one "
            "folder down, which is how a deposit arrives when you unzip it",
            recursive=True),
    Deposit("FR-FCM-ZZLV", "FlowRepository_FR-FCM-ZZLV_files",
            "a 2007 FACSAria, three files, every detector named"),
    Deposit("FR-FCM-Z6UG", "FR-FCM-Z6UG", "a recent deposit, eight files"),
    Deposit("OMIP-018", "OMIP-018/FlowRepository_FR-FCM-ZZ36_files", "seventeen files in one folder"),
    Deposit("OMIP-40", "OMIP-40/FlowRepository_FR-FCM-ZY6D_files", "nine files"),
    Deposit("OMIP-60", "OMIP-60/FlowRepository_FR-FCM-ZYRX_files",
            "thirty three files, so the panel check has to agree across many"),
    Deposit("FR-FCM-Z4KT", "FlowRepository_FR-FCM-Z4KT_files", "sixteen files"),
    Deposit("OMIP-47", "OMIP-47/FlowRepository_FR-FCM-ZYFB_files",
            "a large deposit, so reading the header has to stay cheap"),
    Deposit("FR-FCM-Z244", "FlowRepository_FR-FCM-Z244_files", "twenty eight files over a gigabyte"),
)


class Result(NamedTuple):
    """The outcome of one check.

    Attributes:
        what: The deposit or the case the check ran against.
        step: The recipe or the rule that was checked.
        ok: Whether the check passed.
        detail: What was found, or why the check failed.
        seconds: Wall clock time of the command.
    """

    what: str
    step: str
    ok: bool
    detail: str
    seconds: float


def run_cli(arguments: list[str]) -> tuple[int, str, float]:
    """Call the cytokit CLI and capture everything it wrote.

    Args:
        arguments: The arguments after the `cytokit` name.

    Returns:
        The exit status, the combined output, and the elapsed seconds.
    """
    started = time.monotonic()
    completed = subprocess.run(
        [str(CLI), *arguments],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    elapsed = time.monotonic() - started
    return completed.returncode, completed.stdout + completed.stderr, elapsed


def newest_bundle(root: Path) -> Path | None:
    """Return the bundle folder a recipe just wrote, or None when it wrote none.

    Args:
        root: The folder the recipe was given as `--out`.

    Returns:
        The most recently modified sub folder of `root`.
    """
    folders = [item for item in root.iterdir() if item.is_dir()] if root.exists() else []
    return max(folders, key=lambda item: item.stat().st_mtime) if folders else None


def check_inspect(deposit: Deposit, data: Path, work: Path) -> list[Result]:
    """Run inspect and assert the bundle it wrote.

    Args:
        deposit: The corpus entry.
        data: The deposit folder on the host.
        work: A scratch folder outside the repository.

    Returns:
        One result per rule checked.
    """
    out = work / "inspect"
    out.mkdir(parents=True, exist_ok=True)
    arguments = ["inspect", "--data", str(data), "--out", str(out)]
    if deposit.recursive:
        arguments.append("--recursive")
    status, output, seconds = run_cli(arguments)
    if status != 0:
        return [Result(deposit.name, "inspect", False, first_error(output), seconds)]

    bundle = newest_bundle(out)
    if bundle is None:
        return [Result(deposit.name, "inspect", False, "the recipe wrote no bundle", seconds)]

    results = []
    missing = [name for name in (*BUNDLE_FILES, *INSPECT_FILES) if not (bundle / name).exists()]
    results.append(
        Result(
            deposit.name,
            "inspect",
            not missing,
            f"missing {', '.join(missing)}" if missing else summarise_panel(bundle),
            seconds,
        )
    )

    # A manifest that names no input file cannot be traced back to the data.
    manifest_path = bundle / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        inputs = manifest.get("inputs", {})
        results.append(
            Result(
                deposit.name,
                "manifest",
                bool(inputs),
                f"{len(inputs)} input file(s) checksummed" if inputs else "no input recorded",
                0.0,
            )
        )

    # Every panel is either mass or fluorescence, and the answer changes what
    # the next step is allowed to be. Silence here lets an agent reach for a
    # scatter gate on a CyTOF file.
    stated = [line for line in output.splitlines() if line.startswith("Acquisition:")]
    results.append(
        Result(deposit.name, "acquisition", bool(stated),
               stated[0].strip() if stated else "the kind was not stated", 0.0)
    )

    # A container path in the output is a path the scientist cannot type.
    leaked = [line.strip() for line in output.splitlines() if "/indata" in line or "/outdata" in line]
    results.append(
        Result(
            deposit.name,
            "host paths",
            not leaked,
            f"leaked {leaked[0]}" if leaked else "no container path printed",
            0.0,
        )
    )
    return results


def summarise_panel(bundle: Path) -> str:
    """Describe the panel a bundle recorded, for the report line.

    Args:
        bundle: The bundle folder.

    Returns:
        A short description of the panel size and the marker count.
    """
    panel_lines = (bundle / "panel.csv").read_text().splitlines()
    markers = [line for line in (bundle / "markers.txt").read_text().splitlines() if line.strip()]
    return f"{len(panel_lines) - 1} detectors, {len(markers)} named markers"


def check_scaffold(deposit: Deposit, data: Path, work: Path, recipe: str) -> list[Result]:
    """Run template or definitions and assert both files it writes.

    Args:
        deposit: The corpus entry.
        data: The deposit folder on the host.
        work: A scratch folder outside the repository.
        recipe: Either `template` or `definitions`.

    Returns:
        One result per rule checked.
    """
    target = work / recipe / f"{recipe}.csv"
    target.parent.mkdir(parents=True, exist_ok=True)
    arguments = [recipe, "--data", str(data), "--out", str(target)]
    if deposit.recursive:
        arguments.append("--recursive")
    if recipe == "definitions":
        arguments += ["--populations", "T cells,B cells"]
    status, output, seconds = run_cli(arguments)
    if status != 0:
        return [Result(deposit.name, recipe, False, first_error(output), seconds)]

    notes = target.with_name(f"{recipe}_notes.md")
    missing = [path.name for path in (target, notes) if not path.exists()]
    if missing:
        return [Result(deposit.name, recipe, False, f"missing {', '.join(missing)}", seconds)]

    # Both scaffolds read their own output back with the function that consumes
    # it. That read back is the check, so its verdict is what matters here.
    parsed = "parses" in output or "read back" in output
    rows = len(target.read_text().splitlines()) - 1
    # A scaffold with no rows is a fair answer on a panel it cannot start, for
    # example a mass cytometry file with no scatter. It is only a defect when the
    # recipe writes an empty file and says nothing about why.
    if rows == 0 and "Note:" not in output:
        return [Result(deposit.name, recipe, False, "wrote no row and gave no reason", seconds)]
    reason = " (and said why)" if rows == 0 else ""
    return [
        Result(deposit.name, recipe, parsed,
               f"{rows} row(s){reason}, read back {'ok' if parsed else 'silent'}", seconds)
    ]


def first_error(output: str) -> str:
    """Pull the first useful error line out of a failed run.

    Args:
        output: Everything the command wrote.

    Returns:
        One line naming the failure.
    """
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    lines = [line for line in lines if not line.startswith("Execution halted")]
    for index, line in enumerate(lines):
        if not line.startswith("Error"):
            continue
        # R prints `Error in Call(...) :` and puts the message on the next line
        # when the call is long. The message is what a reader needs.
        message = line.split(":", 1)[1].strip() if ":" in line else ""
        if message:
            return message
        return lines[index + 1] if index + 1 < len(lines) else line
    return lines[-1] if lines else "no output"


def check_refusals(work: Path) -> list[Result]:
    """Assert that the CLI refuses what it has to refuse.

    A recipe that runs on a path that does not exist, or that overwrites a file
    the scientist spent an afternoon filling in, costs more than one that stops.

    Args:
        work: A scratch folder outside the repository.

    Returns:
        One result per refusal checked.
    """
    results = []
    empty = work / "empty"
    empty.mkdir(parents=True, exist_ok=True)

    missing_path = work / "no_such_folder"
    status, output, seconds = run_cli(["inspect", "--data", str(missing_path)])
    named_host_path = str(missing_path) in output and "/indata" not in output
    results.append(
        Result(
            "refusals",
            "a path that does not exist",
            status != 0 and named_host_path,
            "names the path the caller typed" if named_host_path else f"says: {first_error(output)}",
            seconds,
        )
    )

    status, output, seconds = run_cli(["inspect", "--data", str(empty)])
    hints = "--recursive" in output
    results.append(
        Result("refusals", "a folder with no FCS file", status != 0 and hints,
               "suggests --recursive" if hints else first_error(output), seconds)
    )

    # A scaffold that overwrites a filled in template destroys work that only the
    # scientist can redo.
    deposit = DEPOSIT_ROOT / CORPUS[1].path
    kept = work / "kept" / "template.csv"
    kept.parent.mkdir(parents=True, exist_ok=True)
    kept.write_text("alias,pop,parent,dims,gating_method\nmine,mine,root,FSC-A,mindensity\n")
    before = kept.read_text()
    status, output, seconds = run_cli(["template", "--data", str(deposit), "--out", str(kept)])
    results.append(
        Result("refusals", "an --out that already exists",
               status != 0 and kept.read_text() == before,
               "the file is untouched" if kept.read_text() == before else "the file was overwritten",
               seconds)
    )

    # --label belongs to a recipe that writes a timestamped bundle. Accepting it
    # elsewhere teaches an agent a flag that does nothing.
    status, output, seconds = run_cli(
        ["template", "--data", str(deposit), "--out", str(work / "labelled.csv"), "--label", "x"]
    )
    results.append(
        Result("refusals", "--label on a recipe without a bundle", status != 0,
               first_error(output) if status != 0 else "the flag was accepted", seconds)
    )

    status, output, seconds = run_cli(["inspect", "--data", str(deposit), "--nonsense"])
    results.append(
        Result("refusals", "an unknown flag", status != 0,
               first_error(output) if status != 0 else "the flag was accepted", seconds)
    )
    return results


def report(results: list[Result]) -> int:
    """Print the table and return the exit status.

    Args:
        results: Every check that ran.

    Returns:
        0 when every check passed, 1 otherwise.
    """
    width = max(len(item.what) for item in results)
    step_width = max(len(item.step) for item in results)
    current = ""
    for item in results:
        if item.what != current:
            print()
            current = item.what
        mark = "ok  " if item.ok else "FAIL"
        timing = f"{item.seconds:5.1f}s" if item.seconds else "      "
        print(f"{mark} {item.what:<{width}}  {item.step:<{step_width}}  {timing}  {item.detail}")

    failed = [item for item in results if not item.ok]
    print()
    if failed:
        print(f"{len(failed)} of {len(results)} checks failed.")
        return 1
    print(f"All {len(results)} checks passed over {len({item.what for item in results}) - 1} deposits.")
    return 0


def main() -> int:
    """Run the corpus check.

    Returns:
        The process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deposit", action="append", help="check one deposit, repeatable")
    parser.add_argument("--skip-refusals", action="store_true", help="skip the refusal cases")
    options = parser.parse_args()

    corpus = CORPUS
    if options.deposit:
        wanted = set(options.deposit)
        corpus = tuple(item for item in CORPUS if item.name in wanted)
        if not corpus:
            print(f"No corpus deposit matches {', '.join(sorted(wanted))}.", file=sys.stderr)
            print(f"The corpus holds: {', '.join(item.name for item in CORPUS)}", file=sys.stderr)
            return 2

    results: list[Result] = []
    work_root = Path(tempfile.mkdtemp(prefix="cytokit_corpus_"))
    try:
        for deposit in corpus:
            data = DEPOSIT_ROOT / deposit.path
            if not data.exists():
                results.append(Result(deposit.name, "data", False, f"{data} is not here", 0.0))
                continue
            work = work_root / deposit.name
            results += check_inspect(deposit, data, work)
            results += check_scaffold(deposit, data, work, "template")
            results += check_scaffold(deposit, data, work, "definitions")
        if not options.skip_refusals:
            results += check_refusals(work_root / "refusals")
    finally:
        shutil.rmtree(work_root, ignore_errors=True)

    return report(results)


if __name__ == "__main__":
    sys.exit(main())
