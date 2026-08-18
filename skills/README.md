# skills

Thin adapters that teach each host tool how to call `cytokit` and where the
reference documents live. The reasoning and the code stay in `R/` and
`scripts/cytokit/`, so the three adapters do not drift.

- `claude-code/SKILL.md` — a Claude Code Agent Skill, with `name` and
  `description` frontmatter, that triggers on a flow cytometry request.
- `codex/AGENTS.md` — Codex instructions, read from the nearest `AGENTS.md`.
- `opencode/` — opencode instructions plus a `/cytokit` command in
  `command/cytokit.md`.

Each adapter carries the same operating guide: the command surface, the
workflow, what is not built yet, and the rule the repository runs on. The detail
lives in `docs/cytokit_prd.md` and `README.md`.

Update the three together when the `cytokit` surface changes. An adapter that has
drifted from the CLI is worse than no adapter, because it sends an agent to a
command that does not exist.

## Install

| Host | Activate |
|---|---|
| Claude Code | `mkdir -p .claude/skills && ln -s ../../skills/claude-code .claude/skills/cytokit` |
| Codex | `cat skills/codex/AGENTS.md >> AGENTS.md` |
| opencode | Add `skills/opencode/AGENTS.md` to `instructions` in `opencode.json`, and link `command/cytokit.md` into `.opencode/command/` |

See each adapter for its own install note.

## How the skill is tested

Three layers test different things, and only the third one tests the skill.

| Layer | Command | What it proves |
|---|---|---|
| The functions | `podman run --rm -v "$PWD:/work:z" -w /work everything-flow-cytometry:latest Rscript tests/testthat.R` | Each function returns the value it has to return, on a `flowFrame` built in the test |
| The surface | `uv run python scripts/check_cytokit_corpus.py` | Every ready recipe runs on nine deposits that no analysis here reads, and the CLI refuses what it has to refuse |
| The adapters | `python3 scripts/check_skill_adapters.py` | The three adapters name the same commands the CLI dispatches |

The corpus check is the one that finds the faults a unit test cannot. Its nine
deposits carry different instruments, different naming conventions and different
acquisition kinds, because a skill that only works on the data its author used
is not a tool a scientist can use on their own data. Each entry in `CORPUS`
records what it tests that the others do not.

The third layer is the agent, and no script runs it. Copy a few FCS files to a
folder outside this repository, start a session there with only the adapter
loaded, and ask for the analysis in the words a scientist uses. Then read what
the agent did. A skill fails at this layer by writing flowCore code by hand, by
inventing a command that `cytokit list` calls planned, or by naming a cell type
from a detector name that no antibody was mapped to.
