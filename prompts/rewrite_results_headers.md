# Rewrite the results headers of one report in the register of a paper

You are editing one report in `reports/`. Its Introduction, Methods, Discussion
and References are in the right register already. Some of its results headers are
not. They pose a question, promise a contrast, or announce a verdict.

Change those headers and the sentences named below, and nothing else.

Read this whole brief before you change a line.

## The headers

Each header keeps its prefix, whether that is `## Part N:` or a plain `## `. What
changes is the phrase that follows. It becomes a noun phrase naming the subject
of the section. It does not pose a question, promise a contrast or announce a
verdict, and it does not argue with an imagined sceptic.

Write the phrase in lower case after the colon, because that is what every report
in this repository does.

These are the headers to change, by file. A header not listed here does not
change.

| File | Header | The fault |
|---|---|---|
| `flowcap2_challenges.qmd` | `## Part 4: every claim, with a verdict` | Announces the verdict as the point of the section. |
| `oetjen2018_bone_marrow.qmd` | `## Part 1: the finding that changes the whole pipeline` | Promises a finding and withholds it. |
| `oetjen2018_bone_marrow.qmd` | `## Part 3: whether an automated cut can be trusted here` | Shaped as a question, and "trusted" is a verdict. |
| `oetjen2018_bone_marrow.qmd` | `## Part 6: every claim, with a verdict` | As above. |
| `omip43_asc_analysis.qmd` | `## Part 3: automated gating, and why it fails` | Promises an explanation as a hook. |
| `omip43_asc_analysis.qmd` | `## Part 3b: where a density method actually cuts` | "actually" argues with an imagined sceptic. |
| `omip43_asc_analysis.qmd` | `## Part 4: clustering recovers what gating missed` | States the result in the header. |
| `omip43_asc_analysis.qmd` | `## Part 6: every claim, with a verdict` | As above. |
| `yu2021_spectral_mait.qmd` | `## Part 1: what spectral changes about the pipeline` | Shaped as a question. |
| `yu2021_spectral_mait.qmd` | `## Part 2: the cohort, and where the grouping comes from` | Second half shaped as a question. |
| `yu2021_spectral_mait.qmd` | `## Part 3: one cut for the cohort, and the marker where that fails` | Promises a failure as a hook. |
| `yu2021_spectral_mait.qmd` | `## Part 5: what the CD45RA choice changes` | Shaped as a question. |
| `z282_harmonisation.qmd` | `## Part 1: what one pipeline has to survive` | Shaped as a question. |
| `z282_harmonisation.qmd` | `## Part 4: what the four arms say` | Shaped as a question. |

Do not renumber a part. `## Part 3b` keeps its `3b`.

## The prose

In the sections whose headers you change, look for these four patterns and fix
each one you find.

1. A short sentence or a one line paragraph used as a payoff after a longer one.
   Join it to the sentence it depends on.
2. The "not X, it is Y" and "it is not about X, it is about Y" construction.
   State what the thing is.
3. A sentence whose only work is to introduce the next sentence. Delete it.
4. A sentence that tells the reader how to read a number, how to feel about a
   result, or where to look. Give the number and what it measures.

Leave every other paragraph of the file alone. This is a narrow pass.

Keep every admission of failure and every stated limit. A paper reports what did
not work. Do not soften a sentence that says a result is unresolved or that a
method did not work.

## Do not invent a mechanism

This is the rule that the last pass broke, and it is the one that matters most,
because no script catches it.

One report said: "A rule that looks for a second mode either finds nothing, which
is visible, or places the cut inside the negative population and returns a
number, which is not." The rewrite made that specific: "The density minimum rule
finds no second mode on a unimodal distribution. The mixture rule places the cut
inside the negative population and returns a percentage that exceeds the subset
frequency."

It reads better and it is wrong. The table printed directly above that paragraph
shows the density rule placed every cut that selected 38.5 to 62.7 percent of T
cells, and the mixture rule placed one cut only. The rewrite swapped the two
rules. It introduced no new number, no new accession and no changed chunk, so
`check_prose_rewrite.py` passed it as clean.

The rule follows from that. You may cut a sentence, shorten it, or split it. You
may not add specificity that the original did not carry. If the original says
"a rule", do not name which rule. If the original says "a marker", do not name
which marker. If the original gives no cause, do not supply one.

Where a sentence states a cause and you cannot tell from the file whether it is
right, leave that sentence exactly as it is. An unchanged sentence is never a
fault in this pass.

## What you must not change

1. The YAML header. Not one line, not the title, not the subtitle.
2. Any ```{r} chunk. Not the code, not the label, not the chunk options. Do not
   move a chunk to another position in the file.
3. Any `caption =` string inside a chunk, and any `fig-cap:` line.
4. Any block quote that starts with `>`. Those are the words of the paper's
   authors. Never edit quoted material.
5. The `## Introduction`, `## Methods` and `## Discussion` sections and every
   heading inside them.
6. The `## Reproducing this` bash block and the `## References` list.
7. Every number, every percentage, every accession, every file name, every
   marker name, every package name.

## The numbers rule

You may not introduce a number that is not already in the file. You may not round
a number that is already there. If the file says the correlation is 0.972, you
write 0.972, not "above 0.97".

Never write an accession that is not already in the file. An earlier pass
invented `FR-FCM-Z2ZG`, which points at nothing and reads as a fact.

## Rules that still bind you

Everything under "Rules for both modes" in `AGENTS.md`. No em dash and no en
dash. No bold phrase followed by a full stop. No banned opener and no banned
closer. Every sentence is a complete sentence with a subject and a verb. Wrap
every prose line at 80 characters, counting a Greek letter as one character.

## Work in passes

Reading the whole file and then planning the whole rewrite fails on a long file.
The session ends after the reads and it writes nothing. Edit as you go. Take one
header and its section, rewrite it, call the edit tool, and move on.

## Tools

`rg` is not installed. Use `grep`. `grep -nP` is available for a Unicode range,
which is how you find an em dash.

## Acceptance test

Before you report that you are done, confirm all of these and say how you checked.

1. Every header listed for your file is now a noun phrase in lower case.
2. `grep -c '^```'` returns the same count as before you started.
3. No line contains an em dash or an en dash.
4. The YAML header is byte identical to the original.
5. No sentence you wrote names a rule, a marker or a cause that the sentence it
   replaced did not name.

## After a rewrite, read the causal sentences

`check_prose_rewrite.py` cannot fail on a sentence that states the wrong cause,
because a swapped mechanism changes no number and no chunk. Run
`scripts/list_causal_claims.py --numbered-only` on the file and read each
sentence against the table beside it. Sixty sentences across the nine reports
assert a cause, and an audit of them found six that no output supports.

## Do not

Do not run the analysis scripts. Do not render a report. Do not touch any file
other than the single report you were given. Do not commit.
