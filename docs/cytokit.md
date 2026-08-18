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
./cli/cytokit compensate  --data PATH [--controls DIR] [--out DIR]
                                      [--cofactor 150] [--seed INTEGER]
```

Add `--recursive` when the FCS files sit below the data folder you name. That
is how a deposit arrives when you unzip it.

Start with `inspect`. Every later recipe reads what it prints.

## The recipes

| Recipe | State | What it does |
|---|---|---|
| `inspect` | ready | Panel, markers, event counts, acquisition kind, compensation state |
| `template` | ready | An empty openCyto gating template for this panel |
| `definitions` | ready | An empty cell type definitions table for this panel |
| `compensate` | ready | The stored matrix, and whether it lowers the correlation |
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

## Compensate

`compensate` answers whether the data is compensated. The `APPLY COMPENSATION`
keyword does not answer that question. The recipe correlates every marker pair
before and after it applies the stored matrix, on arcsinh transformed values.
Read the change between the rows and not the value on its own, because the
recipe gates no events. A matrix that works lowers the median absolute
correlation. A pair driven far negative is over compensated, so the count uses
the absolute correlation.

The recipe reads every event of one file and it names that file. Measure a
second file when the run spans more than one day or one instrument.

Use `--controls DIR` to compute a matrix from the single stain controls. The
recipe compares that matrix with the stored one and draws both. The bundle
checksums the control files as well as the samples, because the matrix depends
on both.

The recipe refuses a mass cytometry run. A mass cytometer counts an isotope, so
there is no spillover to compensate.

| Flag | What it sets | Default |
|---|---|---|
| `--controls DIR` | The folder of single stain controls | None |
| `--cofactor N` | The arcsinh cofactor | 150 |
| `--unstained PATTERN` | The name pattern of the unstained control | `unstain` |
| `--seed N` | The seed for every random draw | 42 |

Set the seed once and the answer repeats. `flowStats::norm2Filter` draws a
random subset while it gates each control, so an unseeded run moves the matrix
between runs. The seed is recorded in `manifest.json`.

## The build order

`inspect`, `template`, `definitions` and `compensate` are ready. Build the
remaining recipes in this order.

1. `gate`.
2. `cluster` and `annotate`.
3. `proportions` and `compare`.
4. `claims` and `reproduce`.

## What a recipe writes

A recipe that runs an analysis writes one timestamped bundle. The bundle holds
the tables and the figures, and four files that let the result be rebuilt.

| File | What it holds |
|---|---|
| `manifest.json` | The recipe, every argument, the input checksums, the image digest |
| `session_info.txt` | Every package version |
| `REPRODUCE.md` | The command to run again |
| `reader_notes.csv` | What the FCS reader reported, when it reported anything |

`compensate` adds `marker_correlation.csv`, which holds one row per state, and
`stored_matrix.csv` and `stored_matrix.svg` when the file carries a matrix.
With `--controls` it also writes `computed_matrix.csv`, `computed_matrix.svg`,
`control_match.csv` and `matrix_comparison.csv`.

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
