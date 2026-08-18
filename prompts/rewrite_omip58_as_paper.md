# Rewrite the results of `reports/omip58_pytometry.qmd` in the register of a paper

You are editing one file, `reports/omip58_pytometry.qmd`. Its Introduction,
Methods, Discussion and References are in the right register already. The results
sections are not. They carry headers that pose a question or announce a verdict,
and paragraphs that build to a payoff instead of stating a measurement.

Change the register of the six results sections and nothing else.

Read this whole brief before you change a line.

## The headers

Six headers name a part of the results. Each keeps its `## Part N:` prefix,
because the other eight reports in this repository share that form. What changes
is the phrase after the colon. It becomes a noun phrase that names the subject of
the section. It does not pose a question, promise a contrast or announce a
verdict.

| Current header | The fault |
|---|---|
| `## Part 1: what the deposit says about its compensation, and what the data says` | Promises a contrast and withholds it. |
| `## Part 2: the matrix computed from the single stains` | Already a noun phrase. Leave it or tighten it. |
| `## Part 3: the gate and the handoff` | Already a noun phrase. Leave it. |
| `## Part 4: what a one dimensional cut does below CD3` | Shaped as a question. |
| `## Part 5: the clustering` | Already a noun phrase. Leave it. |
| `## Part 6: every claim, with a verdict` | Announces the verdict as the point of the section. |

A paper names the measurement. "Compensation state of the deposited files" names
a subject. "What the deposit says about its compensation, and what the data says"
sets up a reveal.

Do not change any other header. The `## Discussion` subheadings stay as they are.

## The prose

These sentences are in the file now, and each one has to go or change. They are
the pattern to look for, not the whole list.

| Sentence | The fault |
|---|---|
| "Read on its own, that pair says the values were compensated at export and nothing remains to apply." followed by "The events say otherwise." | A three word paragraph as a dramatic payoff. Join the measurement to the statement. |
| "The identity matrix therefore records that no matrix was supplied, not that none is needed." | The "not X, but Y" pattern. State what the matrix is. |
| "The difference shows in the result." | A sentence whose only work is to introduce the next one. |
| "which is the sign a mutually exclusive pair should carry" | Tells the reader how to read a number. Give the number and what CD4 and CD8 are. |
| "A rule that looks for a second mode either finds nothing, which is visible, or places the cut inside the negative population and returns a number, which is not." | Rhetorical balance in place of a statement. Say what each rule returned. |
| "This table is reported and not used." | Reads as an aside to the reader. State it as a method fact. |

Keep every admission of failure and every stated limit. A paper reports what did
not work. Do not soften a sentence that says a result is unresolved, that a
threshold was unstable, or that a label does not name the population the paper
draws. Those are results.

Keep the first person out of it. The passive voice is allowed in a results
section, which is the one place these documents depart from the usual rule.

## What you must not change

1. The YAML header. Not one line, not the title, not the subtitle. They are
   already in the right register.
2. Any ```{r} chunk. Not the code, not the label, not the chunk options. Do not
   move a chunk to another position in the file.
3. Any `caption =` string inside a chunk, and any `fig-cap:` line.
4. The `## Introduction`, `## Methods` and `## Discussion` sections and every
   heading inside them.
5. The `## Reproducing this` bash block.
6. The `## References` list.
7. Every number, every percentage, every correlation, every accession, every
   file name, every marker name, every package name.

## The numbers rule

You may not introduce a number that is not already in the file. You may not round
a number that is already there. If the file says the correlation is 0.687, you
write 0.687, not "about 0.69" and not "close to 0.7".

Never write an accession that is not already in the file. This deposit is
`FR-FCM-ZYRN`. An earlier pass on another report invented `FR-FCM-Z2ZG`, which
points at nothing and reads as a fact.

## Rules that still bind you

Everything under "Rules for both modes" in `AGENTS.md`. No em dash and no en
dash. No bold phrase followed by a full stop. No label, colon and value inside
running text. No banned opener and no banned closer.

Every sentence is a complete sentence with a subject and a verb.

Wrap every prose line at 80 characters. A Greek letter counts as one character.

Vary sentence length. A results section is plain, and plain is not uniform.

## Work in passes

Reading the whole file and then planning the whole rewrite fails on a file this
size. The session ends after the reads and it writes nothing.

Edit as you go. Take one part, rewrite it, call the edit tool, and move on. There
are six parts. That is six edits, or more.

## Tools

`rg` is not installed. Use `grep`. `grep -nP` is available for a Unicode range,
which is how you find an em dash.

## Acceptance test

Before you report that you are done, confirm all of these and say how you
checked.

1. The six `## Part N:` headers are noun phrases and none poses a question or
   names a verdict.
2. `grep -c '^```' reports/omip58_pytometry.qmd` returns the same count as
   before you started.
3. No line contains an em dash or an en dash.
4. The YAML header is byte identical to the original.
5. Every number in the prose appears in the original file.
6. The three sentences quoted in the table above no longer appear.

The author then runs `scripts/check_prose_rewrite.py` on your output. It compares
the file with its git revision and it fails on a changed chunk, a changed block
quote, an em dash, a bullet or a number that is not in the original.

## Do not invent a mechanism

This is the rule that the last pass broke, and it is the one that matters most,
because no script catches it.

The report said: "A rule that looks for a second mode either finds nothing, which
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

## After a rewrite, read the causal sentences

`check_prose_rewrite.py` cannot fail on a sentence that states the wrong cause,
because a swapped mechanism changes no number and no chunk. Run
`scripts/list_causal_claims.py --numbered-only` on the file and read each
sentence against the table beside it. Sixty sentences across the nine reports
assert a cause, and an audit of them found six that no output supports.

## Do not

Do not run the analysis scripts. Do not render the report. Do not touch anything
outside `reports/omip58_pytometry.qmd`. Do not commit.
