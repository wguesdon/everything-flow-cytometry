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
