# cytokit

`cytokit` analyses flow cytometry data that you bring, from your own FCS files.
A library of tested R functions does the compensation, the gating, the
clustering and the statistics. An agent chooses the step with you and calls the
recipe, so nobody writes flowCore or openCyto code by hand.

Every recipe runs inside the pinned Podman image. Your data path is mounted
read only, because an FCS file is a measurement that cannot be taken again.

## The command surface

```bash
./cli/cytokit build                                # once, builds the image
./cli/cytokit list                                 # the recipes and their state
./cli/cytokit inspect     --data PATH              # panel, markers, compensation state
./cli/cytokit template    --data PATH --out FILE   # an empty gating template
./cli/cytokit definitions --data PATH --out FILE   # an empty cell type table
```

Add `--recursive` when the FCS files sit below the folder you name. That is how
a deposit arrives when you unzip it.

Start with `inspect`. Every later recipe reads what it prints.

## The recipes

| Recipe | State | What it does |
|---|---|---|
| `inspect` | ready | Panel, markers, event counts, acquisition kind, compensation state |
| `template` | ready | An empty openCyto gating template for this panel |
| `definitions` | ready | An empty cell type definitions table for this panel |
| `compensate` | planned | Spillover from single stains, and the correlation check |
| `gate` | planned | Run a template, write counts, a gate tree and the flow |
| `cluster` | planned | FlowSOM and UMAP, inside a gate when asked |
| `annotate` | planned | Label each cluster from a definitions table |
| `proportions` | planned | Per sample frequencies joined to a metadata table |
| `compare` | planned | Box plot of proportion by treatment, with the test |
| `claims` | planned | A verdict for every claim in a claims table |
| `reproduce` | planned | Rebuild one of the analyses in this repository |

A recipe that `cytokit list` calls planned does not exist. It is not a command
you can run, and an agent that invents one sends you to a failure that reads
like your own fault.

## The build order

Each step depends on the one before it.

1. `inspect`, and the shared parts of a recipe.
2. `template` and `definitions`, because a panel the tool has not seen needs
   both and nobody arrives with them.
3. `compensate`.
4. `gate`.
5. `cluster` and `annotate`.
6. `proportions` and `compare`.
7. `claims` and `reproduce`.

## What a recipe writes

A recipe that runs an analysis writes one timestamped bundle. The bundle holds
the tables and the figures, and four files that let the result be rebuilt.

| File | What it holds |
|---|---|
| `manifest.json` | The recipe, every argument, the input checksums, the image digest |
| `session_info.txt` | Every package version |
| `REPRODUCE.md` | The command to run again |
| `reader_notes.csv` | What the FCS reader reported, when it reported anything |

`template` and `definitions` write one file that you name, plus a `_notes.md`
beside it. The notes are a separate file because openCyto reads the template
with `read.csv`, whose comment character is empty, so a comment in the CSV
becomes the header and the template stops parsing.

## Two traps this tool watches for

An identity spillover matrix with `APPLY COMPENSATION = TRUE` means that no
matrix was supplied. It does not mean that the values are compensated. Three
deposits in this archive carry that combination.

A mass cytometer counts an isotope, so there is no spillover to compensate and
no forward or side scatter to gate on. `inspect` reads the instrument keyword
and the detector names, and states the acquisition kind before you choose a
step.

A third condition is not a trap but it stops the work. Many deposits leave
`$PnS` empty, and the tool then falls back to the detector name. A detector
name is not an antibody, so supply the mapping before you read a cell type
label.

## How it is tested

| Layer | Command | What it proves |
|---|---|---|
| The functions | `Rscript tests/testthat.R` in the container | Each function returns the right value |
| The surface | `uv run python scripts/check_cytokit_corpus.py` | Every ready recipe runs on nine deposits no analysis reads |
| The adapters | `python3 scripts/check_skill_adapters.py` | The adapters name the commands the CLI dispatches |

The agent is the fourth layer and no script runs it.
`scripts/make_test_study.py` copies a deposit to a folder outside this
repository, for a session that starts with nothing but the skill and the data.

## The agent adapters

`skills/` carries one adapter for Claude Code, one for Codex and one for
opencode. See `skills/README.md` to install them.
