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

    # An FCS file carries no metadata table, so the identity keywords are the
    # only record of the design that ships with the data.
    lines = output.splitlines()
    heading = next((index for index, line in enumerate(lines)
                    if "What the files say about themselves" in line), None)
    roles = []
    if heading is not None:
        for line in lines[heading + 1:]:
            first = line.strip().split(" ", 1)[0]
            if first in {"grouping", "identifier", "constant", "timing"}:
                roles.append(first)
            elif not line.strip():
                break
    said_nothing = "Nothing. Every identity" in output
    results.append(
        Result(deposit.name, "identity", heading is not None and (bool(roles) or said_nothing),
               ", ".join(sorted(set(roles))) if roles else
               ("no keyword recorded, and it says so" if said_nothing else "said nothing about the files"),
               0.0)
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


def check_input_shapes(work: Path) -> list[Result]:
    """Assert the shapes of `--data` that a scientist actually types.

    A scientist points at one file as often as at a folder, and a deposit
    arrives with the files one level down. All three shapes have to work.

    Args:
        work: A scratch folder outside the repository.

    Returns:
        One result per shape checked.
    """
    results = []
    folder = DEPOSIT_ROOT / "FlowRepository_FR-FCM-ZZLV_files"
    one_file = min(folder.glob("*.fcs"))

    out = work / "one_file"
    out.mkdir(parents=True, exist_ok=True)
    status, output, seconds = run_cli(["inspect", "--data", str(one_file), "--out", str(out)])
    counted = "files  1" in output
    results.append(
        Result("input shapes", "--data names one file", status == 0 and counted,
               "reads the one file" if counted else first_error(output), seconds)
    )

    # A file name with a space and a comma in it survives the mount and the
    # label. FR-FCM-ZZCA carries several.
    awkward = work / "a folder, with a comma"
    awkward.mkdir(parents=True, exist_ok=True)
    shutil.copy2(one_file, awkward / "PBMC Stim_5 hr PMA,2f,Ionomycin.fcs")
    out = work / "awkward_out"
    out.mkdir(parents=True, exist_ok=True)
    status, output, seconds = run_cli(["inspect", "--data", str(awkward), "--out", str(out)])
    results.append(
        Result("input shapes", "a path with a space and a comma", status == 0,
               "reads it" if status == 0 else first_error(output), seconds)
    )

    # A scientist names an output folder that does not exist yet. Podman would
    # otherwise create it owned by root, or the recipe would fail on a mount.
    deep = work / "not" / "here" / "yet" / "template.csv"
    status, output, seconds = run_cli(
        ["template", "--data", str(folder), "--out", str(deep)])
    results.append(
        Result("input shapes", "--out in a folder that is not there", status == 0 and deep.exists(),
               "creates the folder" if deep.exists() else first_error(output), seconds)
    )

    # The bundle has to land where it was asked for, and not inside the image.
    wrote = list(out.iterdir()) if out.exists() else []
    results.append(
        Result("input shapes", "--out outside the repository", bool(wrote),
               f"wrote {wrote[0].name}" if wrote else "wrote nothing to the folder given", 0.0)
    )
    return results


def check_marker_mapping(work: Path) -> list[Result]:
    """Assert what definitions does on a panel that names no antibody.

    FR-FCM-ZZCA leaves `$PnS` empty for every detector. A definitions table
    keyed on `APC-A` invites a cell type label that no antibody supports, so the
    recipe has to say so, and it has to accept the mapping when it is given.

    Args:
        work: A scratch folder outside the repository.

    Returns:
        One result per rule checked.
    """
    data = DEPOSIT_ROOT / "FR-FCM-ZZCA"
    results = []

    plain = work / "unnamed.csv"
    plain.parent.mkdir(parents=True, exist_ok=True)
    status, output, seconds = run_cli(
        ["definitions", "--data", str(data), "--recursive", "--out", str(plain),
         "--populations", "T cells"])
    warned = "detector names, not antibodies" in output
    results.append(
        Result("marker mapping", "a panel that names no antibody", status == 0 and warned,
               "says the columns are detectors" if warned else "wrote the file with no warning",
               seconds)
    )

    mapped = work / "mapped.csv"
    status, output, seconds = run_cli(
        ["definitions", "--data", str(data), "--recursive", "--out", str(mapped),
         "--markers", "APC-A=CD3,FITC-A=CD4", "--populations", "T cells"])
    header = mapped.read_text().splitlines()[0] if mapped.exists() else ""
    correct = header.startswith("cell_type,CD3,CD4")
    results.append(
        Result("marker mapping", "a detector=antibody mapping", status == 0 and correct,
               f"writes {header}" if correct else first_error(output), seconds)
    )

    # Losing the mapping leaves a column called CD3 with nothing to say where
    # the name came from.
    notes = mapped.with_name("mapped_notes.md")
    recorded = notes.exists() and "| APC-A | CD3 |" in notes.read_text()
    results.append(
        Result("marker mapping", "the mapping is recorded", recorded,
               "the notes name each detector" if recorded else "the mapping was not written down", 0.0)
    )

    status, output, seconds = run_cli(
        ["definitions", "--data", str(data), "--recursive", "--out", str(work / "bad.csv"),
         "--markers", "CD3"])
    guides = "give the mapping instead" in output
    results.append(
        Result("marker mapping", "a marker name the panel lacks", status != 0 and guides,
               "the error says what to do instead" if guides else first_error(output), seconds)
    )
    return results


def check_compensate(work: Path) -> list[Result]:
    """Assert what compensate says about three deposits that differ.

    The recipe reads the events of one file, so it is slower than the others.
    It runs on three deposits rather than on all nine: one with a real stored
    matrix, one with an identity matrix, and one mass cytometry run. The report
    says so, because a check that covers three of nine and reads like nine is
    worse than no check.

    Args:
        work: A scratch folder outside the repository.

    Returns:
        One result per rule checked.
    """
    results = []

    # A deposit whose matrix is real. Applying it has to lower the correlation.
    real = DEPOSIT_ROOT / "FlowRepository_FR-FCM-ZZLV_files"
    out = work / "real"
    out.mkdir(parents=True, exist_ok=True)
    status, output, seconds = run_cli(["compensate", "--data", str(real), "--out", str(out)])
    bundle = newest_bundle(out)
    table = bundle / "marker_correlation.csv" if bundle else None
    rows = len(table.read_text().splitlines()) - 1 if table and table.exists() else 0
    lowered = "the matrix lowers the correlation" in output
    results.append(
        Result("compensate", "a deposit with a real matrix", status == 0 and rows == 2 and lowered,
               f"{rows} row(s), and the verdict is earned" if lowered else first_error(output), seconds)
    )
    if bundle:
        drawn = [name.name for name in bundle.iterdir() if name.suffix == ".svg"]
        results.append(
            Result("compensate", "the matrix is drawn", bool(drawn),
                   ", ".join(drawn) if drawn else "no heat map was written", 0.0)
        )

    # A deposit with an identity matrix. There is nothing to apply, so the
    # recipe has to say what to do next rather than print one row and stop.
    identity = DEPOSIT_ROOT / "FR-FCM-ZZCA"
    out = work / "identity"
    out.mkdir(parents=True, exist_ok=True)
    status, output, seconds = run_cli(
        ["compensate", "--data", str(identity), "--recursive", "--out", str(out)])
    guides = "--controls" in output and "no matrix was supplied" in output.lower()
    results.append(
        Result("compensate", "a deposit with an identity matrix", status == 0 and guides,
               "says no matrix was supplied and what to do" if guides else first_error(output), seconds)
    )

    # A mass cytometry run. The refusal has to come before the bundle is made.
    mass = DEPOSIT_ROOT / "FlowRepository_FR-FCM-Z244_files"
    out = work / "mass"
    out.mkdir(parents=True, exist_ok=True)
    status, output, seconds = run_cli(["compensate", "--data", str(mass), "--out", str(out)])
    left_behind = newest_bundle(out)
    results.append(
        Result("compensate", "a mass cytometry run", status != 0 and left_behind is None,
               "refuses, and writes no bundle" if left_behind is None else
               f"refused but left {left_behind.name}", seconds)
    )

    # Two runs of one recipe have to give one answer. flowStats draws a random
    # subset while it gates, so an unseeded run moves between runs.
    digests = []
    for run in ("first", "second"):
        out = work / f"seed_{run}"
        out.mkdir(parents=True, exist_ok=True)
        run_cli(["compensate", "--data", str(real), "--out", str(out)])
        bundle = newest_bundle(out)
        table = bundle / "marker_correlation.csv" if bundle else None
        digests.append(table.read_text() if table and table.exists() else run)
    results.append(
        Result("compensate", "two runs give one answer", digests[0] == digests[1],
               "the correlation table is the same" if digests[0] == digests[1] else
               "the two runs disagree, so a seed is missing", 0.0)
    )
    return results


def check_gate(work: Path) -> list[Result]:
    """Run the scaffolded template through gate, which is the intended chain.

    `template` writes a valid empty file and `gate` runs it. If the two disagree
    about the schema, a scientist meets the failure only after they have spent
    an afternoon filling the file in.

    Args:
        work: A scratch folder outside the repository.

    Returns:
        One result per rule checked.
    """
    data = DEPOSIT_ROOT / "FlowRepository_FR-FCM-ZZLV_files"
    work.mkdir(parents=True, exist_ok=True)
    results = []

    template = work / "template.csv"
    status, output, seconds = run_cli(
        ["template", "--data", str(data), "--out", str(template)])
    if status != 0:
        return [Result("gate", "scaffold a template", False, first_error(output), seconds)]

    out = work / "run"
    out.mkdir(parents=True, exist_ok=True)
    status, output, seconds = run_cli(
        ["gate", "--data", str(data), "--template", str(template), "--out", str(out)])
    bundle = newest_bundle(out)
    wanted = ("population_stats.csv", "gate_tree.csv", "gate_tree.svg")
    missing = [name for name in wanted if not bundle or not (bundle / name).exists()]
    results.append(
        Result("gate", "the scaffolded template runs", status == 0 and not missing,
               f"missing {', '.join(missing)}" if missing else "writes the stats, the tree and the figure",
               seconds)
    )
    if bundle and (bundle / "gate_tree.csv").exists():
        rows = (bundle / "gate_tree.csv").read_text().splitlines()
        # PlotGateTree reads a missing parent as the root of the drawing.
        root_first = len(rows) > 1 and rows[1].startswith('"all_events",NA')
        results.append(
            Result("gate", "the tree carries a root row", root_first,
                   "all_events has no parent" if root_first else f"first row reads {rows[1] if len(rows) > 1 else 'nothing'}",
                   0.0)
        )
        # A vector figure that hides a raster is not a vector figure.
        drawing = (bundle / "gate_tree.svg").read_text()
        clean = "<image" not in drawing and "<filter" not in drawing
        results.append(
            Result("gate", "the figure is vector", clean,
                   "no raster inside the svg" if clean else "the svg embeds a raster", 0.0)
        )

    # A template that does not parse has to be refused with the command that
    # writes a valid one.
    broken = work / "broken.csv"
    broken.write_text("this,is,not,a,template\n1,2,3,4,5\n")
    status, output, seconds = run_cli(
        ["gate", "--data", str(data), "--template", str(broken), "--out", str(work / "broken_out")])
    guides = "cytokit template" in output
    results.append(
        Result("gate", "a template that does not parse", status != 0 and guides,
               "the error names the command that writes a valid one" if guides else first_error(output),
               seconds)
    )

    status, output, seconds = run_cli(["gate", "--data", str(data), "--out", str(work / "no_template")])
    results.append(
        Result("gate", "no --template at all", status != 0,
               first_error(output) if status != 0 else "it ran without a template", seconds)
    )
    return results


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
    deposits = ({item.what for item in results}
                - {"refusals", "input shapes", "marker mapping", "compensate", "gate"})
    noun = "deposit" if len(deposits) == 1 else "deposits"
    print(f"All {len(results)} checks passed over {len(deposits)} {noun}.")
    if any(item.what == "compensate" for item in results):
        print("compensate ran on 3 of those deposits, not on all of them, "
              "because it reads every event of a file.")
    return 0


def main() -> int:
    """Run the corpus check.

    Returns:
        The process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deposit", action="append", help="check one deposit, repeatable")
    parser.add_argument("--skip-refusals", action="store_true", help="skip the refusal cases")
    parser.add_argument("--skip-compensate", action="store_true",
                        help="skip the compensate and gate cases, which read whole files")
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
        results += check_input_shapes(work_root / "shapes")
        results += check_marker_mapping(work_root / "mapping")
        if not options.skip_compensate:
            results += check_compensate(work_root / "compensate")
            results += check_gate(work_root / "gate")
        if not options.skip_refusals:
            results += check_refusals(work_root / "refusals")
    finally:
        shutil.rmtree(work_root, ignore_errors=True)

    return report(results)


if __name__ == "__main__":
    sys.exit(main())
