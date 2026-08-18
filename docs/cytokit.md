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
./cli/cytokit gate        --data PATH --template FILE [--out DIR] [--cores N]
                                      [--seed N] [--no-compensate] [--no-transform]
                                      [--no-save-gates]
./cli/cytokit cluster     --gates BUNDLE [--parent POP] | --data PATH
                                      [--markers "CD3,CD4"] [--metaclusters N]
                                      [--grid N] [--events N] [--no-umap] [--seed N]
                                      [--out DIR]
./cli/cytokit annotate    --clusters BUNDLE --definitions FILE [--margin N] [--out DIR]
./cli/cytokit proportions --counts BUNDLE --metadata FILE [--out DIR]
                                      [--sample-column sample]
./cli/cytokit compare     --proportions BUNDLE --group COL [--population NAME]
                                      [--all-populations] [--value COLUMN] [--out DIR]
./cli/cytokit claims      --claims FILE --results FILE [--out DIR]
                                      [--tolerance 0.05]
./cli/cytokit reproduce   [--analysis NAME] [--out DIR]
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
| `gate` | ready | Run a template, write counts and a gate tree |
| `cluster` | ready | FlowSOM and UMAP, inside a gate when asked |
| `annotate` | ready | Label each cluster from a definitions table |
| `proportions` | ready | Per sample frequencies joined to a metadata table |
| `compare` | ready | Box plot of proportion by treatment, with the test |
| `claims` | ready | A verdict for every claim in a claims table |
| `reproduce` | ready | Rebuild one of the analyses in this repository |

A recipe that `cytokit list` calls planned does not exist. It is not a command
you can run, and an agent that invents one sends you to a failure that reads
like your own fault.

## Inspect

`inspect` reports the panel, the markers, the event counts, the acquisition
kind and the compensation state. Every later recipe reads what it prints, so
run it first.

It also reads nine keywords that identify the specimen and the run, and gives
each one a role. A keyword with one value across every file describes the run.
A keyword with a value per file is an identifier. A keyword between the two
splits the folder, and that split is the only candidate grouping that ships
with the data. A timestamp never groups anything. Five of the nine deposits in
the corpus check carry a grouping that no file name reveals.

It states whether the run is mass cytometry or fluorescence, because the answer
changes what the next step is allowed to be. It warns when a file leaves `$PnS`
empty, because the tool then falls back to the detector name and a detector
name is not an antibody.

| Flag | What it sets | Default |
|---|---|---|
| `--recursive` | Look for FCS files below the folder | Off |
| `--out DIR` | Where the bundle goes | `output` |
| `--label NAME` | The name inside the bundle folder | The data folder |

The bundle holds `files.csv`, `samples.csv`, `panel.csv`, `naming.csv`,
`markers.txt` and `reader_notes.csv`.

## Template

`template` writes a valid empty openCyto gating template for your panel. When
the panel carries FSC-A, SSC-A and FSC-H, the file already holds the two
scatter rows, so you can see the shape before you add your own populations. The
recipe reads the file back with the function that will consume it and reports
whether it parses.

The guidance goes in a `_notes.md` file beside the CSV and never inside it.
openCyto reads the template with `read.csv`, whose comment character is empty,
so a comment line would become the header and the template would stop parsing.

| Flag | What it sets | Default |
|---|---|---|
| `--out FILE` | The file to write | Required |
| `--recursive` | Look for FCS files below the folder | Off |
| `--force` | Overwrite a file that is already there | Off |

The recipe refuses to overwrite an existing file without `--force`. A template
holds decisions that took work.

## Definitions

`definitions` writes a valid empty cell type definitions table whose columns are
the markers of your panel. Fill a cell with `pos`, `neg` or `high`, and leave it
empty when the marker does not decide that population.

`--markers` takes a selection such as `CD3,CD4`. When the file names no
antibody, it takes a mapping such as `APC-A=CD3,PE-A=CD4` instead, and the
recipe writes the antibody names. It records that mapping in the notes, because
a column called CD3 with nothing to say which detector it came from is not
provenance. It warns whenever a column is still a detector name.

| Flag | What it sets | Default |
|---|---|---|
| `--out FILE` | The file to write | Required |
| `--markers` | A selection, or a detector to antibody mapping | Every marker |
| `--populations` | The rows to write | None |
| `--recursive` | Look for FCS files below the folder | Off |
| `--force` | Overwrite a file that is already there | Off |

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

## Gate

`gate` runs an openCyto gating template that the scientist fills in. The
template records the scientist's decision. The recipe applies the stored
compensation matrix before it applies the logicle transform. A gate drawn on
uncompensated values counts spillover as signal. openCyto cuts values on a
transformed scale.

Use `--no-compensate` to skip the stored matrix. Use `--no-transform` when the
values already sit on a gating scale. Run `cytokit template` first to get a
valid empty template.

| Flag | What it sets | Default |
|---|---|---|
| `--template FILE` | The openCyto gating template | Required |
| `--cores N` | The number of cores | 1 |
| `--seed N` | The seed | 42 |
| `--no-compensate` | Skip the stored matrix | False |
| `--no-transform` | Skip the logicle transform | False |
| `--no-save-gates` | Do not save the gated hierarchy | False |

The recipe writes `population_stats.csv` with one row for each sample and
population. It writes `gate_tree.csv` with the parent of each population. It
also writes `gate_tree.svg`.

The recipe names every gate that keeps more than 99.5 percent of its parent. It
also names every gate that keeps less than 1 percent of its parent. A gate that
keeps every event did not cut. A gate that keeps almost no events cut in the
wrong place. Both cases can look like a result in a table.

When the input has more than one file, the recipe reports the frequency spread
of each population across samples as a coefficient of variation.

`gate` saves the gated hierarchy in its bundle as a `gating_set` folder.
`cluster` reads this folder and clusters events inside one of its gates.
`--no-save-gates` stops this save.

## Cluster

`cluster` groups events with FlowSOM and draws them with UMAP. It does not name
a cluster, because a name is a claim about an antibody. The scientist owns that
claim.

Use `--gates BUNDLE` and `--parent POP` to use events from a saved hierarchy.
Use `--data PATH` to use events from FCS files. Clustering inside a gate is the
usual case, because a clustering over debris spends its clusters on debris.

| Flag | What it sets | Default |
|---|---|---|
| `--gates BUNDLE` | The gate bundle that holds `gating_set` | None |
| `--parent POP` | The population in the saved hierarchy | The last population |
| `--data PATH` | The FCS file or folder | None |
| `--markers LIST` | The channels for clustering | All marker channels |
| `--metaclusters N` | The number of clusters | 12 |
| `--grid N` | The FlowSOM grid size | 10 |
| `--events N` | The maximum event subsample | 50000 |
| `--no-umap` | Skip the UMAP embedding | False |
| `--seed N` | The seed | 42 |

The recipe replaces a detector name with the marker name when the file supplies
one. It names every cluster with fewer than 20 events. A median from that few
events is noise.

The recipe writes `cluster_medians.csv` and `cluster_medians.svg`. It writes
`umap.csv` and `umap.png` unless `--no-umap` is set. The UMAP uses a raster,
because a vector records every point.

## Annotate

`annotate` puts a cell type name on every cluster. The name comes from the
scientist's definitions table.

`--clusters BUNDLE` identifies the bundle that `cluster` writes.
`--definitions FILE` identifies the definitions table. Both flags are
required.

| Flag | What it sets | Default |
|---|---|---|
| `--clusters BUNDLE` | The cluster bundle | Required |
| `--definitions FILE` | The cell type definitions table | Required |
| `--margin N` | The close-call margin | 0.1 |

The recipe scores every cluster against every definition. The margin between
the best score and the second-best score decides whether a label is a fact or a
close call. It names every label below the threshold. Report one of these
labels as a candidate and not as a cell type.

The recipe names every definition that wins no cluster. The population is
absent, or the definition does not suit the panel. It reports how many markers
from each definition occur in the clustering. A definition scores on the other
markers when the clustering lacks one, so its label is weaker than it looks.

The recipe writes `cluster_labels.csv`, `cell_type_summary.csv` and
`cell_type_summary.svg`. It writes `close_calls.csv` when a close call occurs.

## Proportions

`proportions` joins a population count table to a metadata table that the
scientist writes. An FCS file carries no treatment, so the metadata table comes
from the scientist.

`--counts BUNDLE` identifies the bundle that `gate` writes. It can also
identify a CSV with `sample`, `population` and `percent_of_parent` columns.
`--metadata FILE` identifies the metadata table. Both flags are required.

| Flag | What it sets | Default |
|---|---|---|
| `--counts BUNDLE` | The gate bundle or count CSV | Required |
| `--metadata FILE` | The metadata table | Required |
| `--sample-column COLUMN` | The metadata column that holds the file name | `sample` |

The recipe matches sample names on the base name. A path or extension on one
side only still matches. The recipe keeps and reports a sample that the
metadata does not name. A dropped sample can remove a replicate from a group.

The recipe reports which metadata columns can carry a test. A test needs two
groups and two samples in each group. The recipe names each column that does
not meet this condition.

The recipe writes `proportions.csv` and `design.csv`.

## Compare

`compare` draws one population frequency by group and puts the test on the
figure. The recipe reads the bundle that `proportions` writes.

`--proportions BUNDLE` identifies the proportions bundle. `--group COLUMN`
identifies the metadata column that splits the samples. Both flags are
required. Use `--population NAME` for one population. Use `--all-populations`
for every population.

| Flag | What it sets | Default |
|---|---|---|
| `--proportions BUNDLE` | The proportions bundle | Required |
| `--group COLUMN` | The metadata column that splits samples | Required |
| `--population NAME` | One population to compare | None |
| `--all-populations` | Compare every population | False |
| `--value COLUMN` | The frequency column | `percent_of_parent` |

The design chooses the test. Two groups use a Wilcoxon rank sum test. More
than two groups use a Kruskal Wallis test. A cell frequency is bounded, skewed
and measured on a small sample count, so the recipe does not use a t test.

A group with one sample carries no spread, so no test runs. The recipe still
draws the figure and shows one point for each sample. Report this figure as a
picture of the samples and not as a difference.

The figure always draws the points on top of the box. A box over three points
draws a quartile that three points cannot support.

When the recipe runs more than one test from one table, `tests.csv` adds the
Benjamini and Hochberg adjusted column. Read the adjusted column and not the
raw column.

The recipe writes `tests.csv`, one group summary CSV for each population and
one SVG for each population.

## Claims

`claims` puts one verdict on every claim, which is supported, contradicted or
unresolved. Unresolved is a real answer and not a failure to answer. A claim
whose measure nobody computed has not been contradicted, and recording it as
anything else overstates what was checked.

The claims table needs the columns `claim_id`, `claim`, `measure`, `test` and
`expected`. The results table needs the columns `measure` and `value`.

| Test | What it asks | What `expected` holds |
|---|---|---|
| `at_least` | Is the value at or above a floor | One number |
| `at_most` | Is the value at or below a ceiling | One number |
| `between` | Is the value inside a range | Two numbers, separated by a comma |
| `equals` | Is the value near a target | One number |
| `greater_than` | Is the value above another measure | The name of that measure |
| `present` | Was the measure computed at all | Nothing |

| Flag | What it sets | Default |
|---|---|---|
| `--claims FILE` | The claims table | Required |
| `--results FILE` | The measured values | Required |
| `--tolerance N` | How far an `equals` claim may sit from its target | 0.05 |
| `--out DIR` | Where the bundle goes | `output` |

The bundle holds `verdicts.csv`, `verdict_tally.csv` and `verdicts.svg`.

## Reproduce

`reproduce` lists the analyses in this repository and runs one. This is the
secondary use of the tool. The first use is your own data.

With no `--analysis` it lists every analysis with its script, its deposit and
its report. `docs/analyses.csv` holds that map.

| Flag | What it sets | Default |
|---|---|---|
| `--analysis NAME` | The analysis to run | List them and stop |
| `--out DIR` | Where the bundle goes | `output` |

The bundle holds `outcome.csv` and `run_log.txt`. An analysis reads from
`data/` and writes into `output/`.

## Everything is built

Every recipe in the table above is ready. `cytokit list` is the authority, and
`scripts/check_skill_adapters.py` fails when this document, the three adapters
or `README.md` disagree with it.

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

`gate` adds `population_stats.csv`, `gate_tree.csv` and `gate_tree.svg`. It
writes `gates_to_check.csv` when a gate keeps more than 99.5 percent or less
than 1 percent of its parent. With more than one file, it also writes
`population_spread.csv`.

`cluster` adds `cluster_medians.csv` and `cluster_medians.svg`. It adds
`umap.csv` and `umap.png` unless `--no-umap` is set.

`annotate` adds `cluster_labels.csv`, `cell_type_summary.csv` and
`cell_type_summary.svg`. It adds `close_calls.csv` when a close call occurs.

`proportions` adds `proportions.csv` and `design.csv`.

`compare` adds `tests.csv`, one `_by_group.csv` file for each population and
one SVG for each population.

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
