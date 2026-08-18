# cytokit

Analyse flow cytometry data the scientist brings, from their own FCS files. A
library of tested R functions does the compensation, the gating, the clustering
and the statistics. Your job is to choose the step with the scientist and call
the recipe. Do not write flowCore, openCyto or ggplot code by hand. Write
original code only when no recipe fits, and say so when you do.

## The one command surface

Everything goes through `./cli/cytokit`, which runs the pinned Podman image.
Never call R or Rscript directly.

```
cytokit build                                   # once, builds the image
cytokit list                                    # the recipes and their state
cytokit inspect     --data PATH                 # panel, markers, compensation state
cytokit template    --data PATH --out FILE      # an empty gating template for this panel
cytokit definitions --data PATH --out FILE      # an empty cell type table for this panel
cytokit compensate  --data PATH [--controls DIR] [--out DIR]
                                    [--cofactor 150] [--unstained PATTERN]
                                    [--seed INTEGER]
cytokit gate        --data PATH --template FILE [--out DIR] [--cores N]
                                    [--seed N] [--no-compensate] [--no-transform]
                                    [--no-save-gates]
cytokit cluster     --gates BUNDLE [--parent POP] | --data PATH
                                    [--markers "CD3,CD4"] [--metaclusters N]
                                    [--grid N] [--events N] [--no-umap] [--seed N]
                                    [--out DIR]
cytokit annotate    --clusters BUNDLE --definitions FILE [--margin N] [--out DIR]
cytokit proportions --counts BUNDLE --metadata FILE [--out DIR]
                                    [--sample-column sample]
cytokit compare     --proportions BUNDLE --group COL [--population NAME]
                                    [--all-populations] [--out DIR]
cytokit claims      --claims FILE --results FILE [--out DIR]
                                    [--tolerance 0.05]
cytokit reproduce   [--analysis NAME] [--out DIR]

# Add --recursive to an FCS data command when the files sit below the folder
# you name, which is how a deposit arrives when you unzip it.
```

The data path is mounted read only. An FCS file is a measurement that cannot be
taken again.

`compensate` answers whether the data are compensated. The `APPLY
COMPENSATION` keyword does not answer this question. It correlates every marker
pair before and after it applies the stored matrix to arcsinh transformed
values. Read the change between rows because the recipe does not gate events.
A good matrix lowers the median absolute correlation. A pair driven far
negative is over compensated, so the count uses the absolute correlation.

Use `--controls DIR` to compute a matrix from the single stain controls. The
recipe compares that matrix with the stored matrix. `--cofactor` sets the
arcsinh cofactor, which defaults to `150`. `--unstained PATTERN` sets the
pattern for the unstained control, which defaults to `unstain`. `--seed
INTEGER` changes the seed, which defaults to `42`. `flowStats::norm2Filter`
draws a random subset while it gates each control, so a different seed
changes the computed matrix. Do not run `compensate` on mass cytometry data. A
mass cytometer counts an isotope, so there is no spillover to compensate.

`gate` runs an openCyto gating template that the scientist fills in. The
template records the scientist's decision. The recipe applies the stored
compensation matrix before it applies the logicle transform. A gate drawn on
uncompensated values counts spillover as signal. openCyto cuts values on a
transformed scale.

Use `--no-compensate` to skip the stored matrix. Use `--no-transform` when the
values already sit on a gating scale. The `--template FILE` flag is required.
The `--cores N` flag sets the number of cores. The `--seed N` flag sets the
seed.

Run `cytokit template` first to get a valid empty template. The recipe writes
`population_stats.csv` with one row for each sample and population. It writes
`gate_tree.csv` with the parent of each population. It also writes
`gate_tree.svg`.

The recipe names every gate that keeps more than 99.5 percent of its parent. It
also names every gate that keeps less than 1 percent of its parent. A gate that
keeps every event did not cut. A gate that keeps almost no events cut in the
wrong place. Both cases can look like a result in a table.

When the input has more than one file, the recipe reports the frequency spread
of each population across samples as a coefficient of variation.

`gate` saves the gated hierarchy in the bundle as a `gating_set` folder. Use
that folder with `cluster` to cluster events inside a gate. Use
`--no-save-gates` to stop this save.

`cluster` groups events with FlowSOM and draws them with UMAP. It does not name
a cluster, because a name is a claim about an antibody. Use `--gates BUNDLE`
and `--parent POP` to cluster events from a saved hierarchy. Use `--data PATH`
to cluster events from FCS files. Clustering inside a gate is the usual case,
because a clustering over debris spends its clusters on debris.

`--metaclusters N` sets the number of clusters and defaults to `12`. `--grid
N` sets the FlowSOM grid and defaults to `10`. `--events N` caps the subsample
and defaults to `50000`. `--markers` narrows the channels. `--no-umap` skips
the embedding. `--seed N` sets the seed.

The recipe replaces a detector name with the marker name when the file supplies
one. It names every cluster with fewer than 20 events. A median from that few
events is noise. The recipe writes `cluster_medians.csv` and
`cluster_medians.svg`. It writes `umap.csv` and `umap.png` unless `--no-umap`
skips the embedding. The UMAP uses a raster, because a vector records every
point.

`annotate` puts a cell type name on every cluster. The name comes from the
scientist's definitions table. `--clusters BUNDLE` and `--definitions FILE`
are required. `--margin N` sets the close-call margin and defaults to `0.1`.

The recipe scores every cluster against every definition. The margin between
the best score and the second-best score decides whether a label is a fact or a
close call. It names every label below the margin threshold. Report one of
those labels as a candidate and not as a cell type.

The recipe names every definition that wins no cluster. The population is
absent, or the definition does not suit the panel. It reports how many markers
from each definition occur in the clustering. A definition scores on its other
markers when the clustering lacks one, so its label is weaker than it looks.
The recipe writes `cluster_labels.csv`, `cell_type_summary.csv`,
`cell_type_summary.svg` and `close_calls.csv` when a close call occurs.

## Workflow

1. Run `cytokit inspect` first. Read the real markers, event counts and
   compensation state. Do not guess them from a file name.
2. Read what it says about marker names. Many deposits leave `$PnS` empty, and
   the tool then falls back to the detector name. A detector name is not an
   antibody, so ask the scientist for the mapping before any cell type label
   means anything. Pass what they give you to `cytokit definitions` as
   `--markers "APC-A=CD3,PE-A=CD4"`, which writes the antibody names and
   records in the notes that they came from the command line and not from the
   data.
3. Read what it says about the files themselves. `inspect` reads the keywords
   that identify the specimen and the run, and reports which ones split the
   folder. A keyword that splits it is the only candidate grouping that ships
   with the data. A file name is a guess, so confirm any grouping you read from
   one.
4. Read the acquisition line. A mass cytometry file carries no scatter and no
   spillover, so a scatter gate and a compensation step are both wrong on it.
   Gate it on the DNA channel and on event length instead.
5. Read what it says about the compensation state. An identity matrix with
   `APPLY COMPENSATION = TRUE` means no matrix was supplied. It does not mean
   the values are compensated. Three deposits in this repository carry that
   trap.
6. Run `cytokit compensate` after you read the compensation state. Read the
   verdict it prints. A matrix that lowers the median absolute correlation is
   doing work, and a matrix that changes almost nothing means the values are
   already compensated. The recipe measures one file of the folder and names
   it, so measure a second file when the run spans more than one day or one
   instrument.
7. Confirm the gate design with the scientist. Record how many populations the
   design has, which marker separates each one and the negative control.
8. Run `cytokit template` first to get a valid empty template. Fill in the
   template with the scientist. Run `cytokit gate`. Read the gates that the
   recipe names before you report a population.
9. Run `cytokit cluster` inside the gate that holds the events of interest. Read every cluster with fewer than 20 events as noise. Do not name a cluster before the scientist supplies a cell type definition.
10. A cell type definitions table is per panel, and the scientist will not have one. Run `cytokit definitions` to get a valid empty file. Fill it in with the scientist. The recipe reads the file back and reports whether it parses.
11. Run `cytokit annotate` with the cluster bundle and the definitions table. Report every close call as a candidate and not as a cell type.
12. Run `cytokit proportions` with the gate bundle and the metadata table. The
    metadata table comes from the scientist, because an FCS file carries no
    treatment. Use `--sample-column` when its file name column is not
    `sample`. Read every unmatched sample that the recipe reports. The recipe
    keeps that sample, because a dropped sample can remove a replicate from a
    group. Read `design.csv` before a comparison. A test needs two groups and
    two samples in each group.
13. Run `cytokit compare` with the proportions bundle and the metadata group
    column. Use `--population` for one population or `--all-populations` for
    every population. Use `--value` when the frequency column is not
    `percent_of_parent`. Two groups use a Wilcoxon rank sum test. More than two
    groups use a Kruskal Wallis test. A group with one sample has no test, but
    the figure still shows each sample. Report that figure as a picture of the
    samples and not as a difference. Read `p_value_adjusted` and not
    `p_value` when more than one test runs.
## The rule this repository runs on

Measure before you decide, and record what you measured.

Every recipe writes one bundle holding its tables, its figures, a
`manifest.json` with every argument and the checksum of every input, a
`session_info.txt` with every package version, and a `REPRODUCE.md`. A number
you report to the scientist has to be traceable to one of those folders.

A threshold that a rule fits is not a threshold that is right. When a recipe
records the alternatives it did not use, show them to the scientist rather than
reporting only the one it chose.

## Read more

- `docs/cytokit.md` — the plan, the recipe table and the build order.
- `README.md` — the ten worked analyses in this repository and what each asks.
- `docs/literature.md` — the paper behind each package.

## Install

Add this file to `instructions` in `opencode.json`, and link the command:

```bash
mkdir -p .opencode/command
ln -s ../../skills/opencode/command/cytokit.md .opencode/command/cytokit.md
```

`opencode.json` at the repository root:

```json
{
  "$schema": "https://opencode.ai/config.json",
  "instructions": ["skills/opencode/AGENTS.md"]
}
```
