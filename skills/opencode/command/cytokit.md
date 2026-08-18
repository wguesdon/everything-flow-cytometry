---
description: Analyse flow cytometry FCS files with cytokit.
agent: build
---

The user wants to analyse flow cytometry data. Use cytokit. Follow the adapter in
`skills/opencode/AGENTS.md`: read the panel first, confirm the design with the
user, then call `./cli/cytokit`. Do not hand-write flowCore or openCyto code.

Data path: $ARGUMENTS

Start by reading the panel:

`!./cli/cytokit inspect --data $ARGUMENTS`

Then read three things out of that output before you go further. The markers,
because many deposits leave `$PnS` empty and the tool falls back to the detector
name, which is not an antibody. The compensation state, because an identity
matrix with `APPLY COMPENSATION = TRUE` means no matrix was supplied rather than
compensated values. And whether every file carries the same detectors, because
one gating template cannot cover two panels.

A gating template and a cell type definitions table are per panel and the user
will not have one. Draft both with them using `cytokit template` and
`cytokit definitions`, which write a valid empty file and read it back to check
that it parses.

Run `./cli/cytokit list` to see which recipes are built. Do not invent a command
it does not show.
