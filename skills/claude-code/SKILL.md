---
name: cytokit
description: >-
  Analyse flow cytometry FCS files: read a panel, compute a compensation matrix,
  gate a hierarchy, cluster and identify cell types, and compare cell proportions
  between treatments. Use when the user has FCS files and wants any of those
  steps. Drives the cytokit CLI. Do not hand-write flowCore or openCyto code.
allowed-tools: Bash, Read, Write, Edit, Glob, Grep
---

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

# Add --recursive to any of the four when the FCS files sit below the folder
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
7. A gating template and a cell type definitions table are per panel, and the
   scientist will not have one. Run `cytokit template` and `cytokit
   definitions` to get a valid empty file, then fill it in with them. Both
   recipes read the file back with the function that will consume it and report
   whether it parses.
8. Confirm the design with the scientist before you draw a gate. How many
   populations, which marker separates each one, and what the negative control
   is.

## What is not built yet

`cytokit list` marks a recipe `ready` or `planned`. The planned recipes are
`gate`, `cluster`, `annotate`, `proportions`, `compare`, `claims` and
`reproduce`. Do not invent a command that `cytokit list` does not show. If the
scientist asks for one, say it is not built and point at `docs/cytokit.md`,
which holds the plan and the build order.

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

The skill source lives in this repository. To activate it as a project skill:

```bash
mkdir -p .claude/skills
ln -s ../../skills/claude-code .claude/skills/cytokit
```

For a personal skill in every project, link it under your home config instead:
`ln -s "$PWD/skills/claude-code" ~/.claude/skills/cytokit`.
