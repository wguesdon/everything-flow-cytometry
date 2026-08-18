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
```

The data path is mounted read only. An FCS file is a measurement that cannot be
taken again.

## Workflow

1. Run `cytokit inspect` first. Read the real markers, event counts and
   compensation state. Do not guess them from a file name.
2. Read what it says about marker names. Many deposits leave `$PnS` empty, and
   the tool then falls back to the detector name. A detector name is not an
   antibody, so ask the scientist for the mapping before any cell type label
   means anything.
3. Read what it says about the compensation state. An identity matrix with
   `APPLY COMPENSATION = TRUE` means no matrix was supplied. It does not mean
   the values are compensated. Three deposits in this repository carry that
   trap.
4. A gating template and a cell type definitions table are per panel, and the
   scientist will not have one. Run `cytokit template` and `cytokit
   definitions` to get a valid empty file, then fill it in with them. Both
   recipes read the file back with the function that will consume it and report
   whether it parses.
5. Confirm the design with the scientist before you draw a gate. How many
   populations, which marker separates each one, and what the negative control
   is.

## What is not built yet

`cytokit list` marks a recipe `ready` or `planned`. The planned recipes are
`compensate`, `gate`, `cluster`, `annotate`, `proportions`, `compare`, `claims`
and `reproduce`. Do not invent a command that `cytokit list` does not show. If
the scientist asks for one, say it is not built and point at
`docs/cytokit_prd.md`, which holds the plan and the build order.

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

- `docs/cytokit_prd.md` — the plan, the recipe table and the build order.
- `README.md` — the ten worked analyses in this repository and what each asks.
- `docs/literature.md` — the paper behind each package.

## Install

Codex reads the nearest `AGENTS.md`. This repository already has one at the
root, which carries the coding rules. Append this file to it:

```bash
cat skills/codex/AGENTS.md >> AGENTS.md
```

Keep the two in step. When the `cytokit` surface changes, update every adapter
in `skills/` together, because an adapter that has drifted is worse than none.
