# everything-flow-cytometry

This repository collects the literature, the software packages and the example code for computational
flow cytometry analysis in R and Python.

## Why this repository exists

Panels grew faster than the analysis method. A three colour assay gives 3 possible bi-axial plots. A
ten colour panel gives 45. A 45 marker spectral panel gives 990. Manual gating does not scale to that
number, and the operator who draws the gates is a large source of variation between laboratories.

The software to replace manual gating exists. Most of it is command line R, and most clinical
laboratories do not use it. This repository puts the reading list, the package list and runnable
example code in one place.

## Status

The repository is at an early stage. The table below states what is ready and what is not.

| Item | State |
|---|---|
| `docs/literature.md` | Ready. Every citation was verified against a retrieved record on 2026-08-13. |
| `docs/packages.md` | Partial. The package list is complete for the tools in daily use. The version and maintenance status is checked for some entries only. |
| `docs/datasets.md` | Ready. Every repository count was checked on 2026-08-13. |
| `data/` and `sync.sh` | Ready. About 103 GB of FCS files and reference code, held in S3. |
| `examples/` | Empty. It needs the analysis scripts and a container. |
| Agent skills | Not started. They come after the example code runs. |

## Contents

| Path | Content |
|---|---|
| `docs/literature.md` | Key reviews, the gating variability studies and the data sharing papers |
| `docs/packages.md` | R and Python packages, grouped by the step they perform |
| `docs/datasets.md` | Public repositories that hold FCS files, with their current state |
| `docs/data_catalog.md` | Every folder in `data/` with its size, so you can choose what to pull |
| `sync.sh` | Push and pull `data/` to and from S3, in whole or in part |
| `config.sh` | The bucket URI and the storage class |
| `scripts/import_from_wd1.sh` | Copy the archive from the WD1 external drive |
| `scripts/make_data_catalog.sh` | Rewrite `docs/data_catalog.md` from the local `data/` folder |

## The data

`data/` holds about 103 GB and it is gitignored. Git holds the code and the
documentation. S3 holds the data, in `s3://wguesdon-flow-cytometry`, with versioning
enabled and all four public access blocks on.

Most of the archive is FlowRepository downloads. That site stopped accepting new
experiments in 2025 and its TLS certificate expired on 18 March 2023, so several of
these accessions are hard to download again. Treat this archive as the working copy.

A full pull transfers about 100 GB. Read `docs/data_catalog.md` first, then pull the
folders you need.

```bash
./sync.sh catalog                                  # size of each folder, read from S3
./sync.sh pull datasets/flowrepository/FR-FCM-ZZZU  # one accession
./sync.sh pull literature repositories             # two folders at once
./sync.sh push datasets/flowrepository/FR-FCM-Z282  # upload one folder
./sync.sh push                                     # upload everything
```

| Flag | Effect |
|---|---|
| `--dry-run` | Show the changes and transfer nothing |
| `--delete` | Remove remote files that are missing locally. Push only, one folder at a time, and it asks first. |

`--delete` is refused on `pull`. It would erase local files that are absent from the
bucket, and `data/` is gitignored, so this working copy can be the only copy. The
bucket has versioning enabled, so a wrong `push --delete` can be undone. A wrong
`pull --delete` cannot.

### Run a long transfer in tmux

A full push moves about 100 GB. Start it in a tmux session so the transfer survives
a closed terminal or a dropped SSH connection.

```bash
cd /mnt/data/Github/everything-flow-cytometry
tmux new -s flow_push
./sync.sh push 2>&1 | tee logs/push_$(date +%Y_%m_%d).log
```

Press `Ctrl-b` then `d` to detach. The transfer continues.

```bash
tmux attach -t flow_push     # go back to it
tmux ls                      # list the sessions
tail -f logs/push_*.log      # watch the log without attaching
```

To start it detached in one command:

```bash
tmux new -d -s flow_push -c /mnt/data/Github/everything-flow-cytometry \
  './sync.sh push 2>&1 | tee logs/push_'"$(date +%Y_%m_%d)"'.log'
```

`aws s3 sync` compares each object before it transfers, so an interrupted push is
safe. Run the same command again and it continues from where it stopped.

## Planned work

1. Add one worked example on a public dataset. The example reads FCS files, applies compensation and
   transformation, runs quality control, gates the data and finds populations.
2. Put the example in a container so the result is the same on every machine. Use Podman.
3. Add agent skills for Claude Code, Codex and OpenCode. Each skill runs the example code and reports
   the result.

## Provenance

The citations and the repository counts come from a literature search that was run on 2026-08-13.
Each claim was checked against Europe PMC, PubMed, Crossref, Semantic Scholar or the PMC full text. A
claim that failed verification is not in this repository.
