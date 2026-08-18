# Rewrite one report in the register of a scientific paper

You are editing the prose of one Quarto report in `reports/`. The prose is
already written. Your job is to change its register and its structure, not to
write it again from nothing.

Three things change. The register becomes that of a methods and results section
in a journal article. A Methods section is added. Every reference to another
report in this repository is removed, so the document stands alone.

Read this whole brief before you change a line.

## Register

Write as a methods and results section reads. State what was done, what was
measured, and what the measurement showed. Nothing is sold to the reader.

Remove every trace of salesmanship. These are the ones that appear in the
current prose.

1. A superlative about the work itself. "This is the closest agreement anywhere
   in this repository" becomes a statement of the number.
2. A sentence that tells the reader how to feel about a result. "That is the
   useful result", "That is worth stating plainly", "The margin is not a close
   call", "worth naming rather than smoothing over".
3. A rhetorical question. "How close does one automated rule get?" becomes a
   statement of what was tested.
4. An instruction to the reader. "Read the last column", "Read the passes as a
   ladder of tolerance", "Compare the two tables". Say what the column shows
   instead.
5. A sentence whose only work is to introduce the next sentence.
6. A verdict word in place of a number. "impossible", "not credible",
   "catastrophic", "beautiful", "elegant".

Keep the first person plural out of it unless the original already had it. A
methods section says what was done, and it does not need an actor in every
sentence. The passive voice is allowed here, which is the one place these
documents depart from the usual rule, because a methods section is exactly the
case where the actor does not matter.

Keep every admission of failure and every stated limit. A paper reports what did
not work. Do not soften "the template is wrong in both directions" into anything
milder, and do not delete a sentence that says a question is unresolved.

## The Methods section

Add one `## Methods` section to each report. Put it after the introductory
section that explains the dataset and before the first section that reports a
result.

It gathers the method that is currently scattered through the results. Write it
as continuous prose under these subheadings, and drop a subheading when the
report has nothing for it.

```
## Methods

### The dataset
### Preprocessing
### Gating
### Statistics
### Software
```

The Methods section describes. It reports no result and it carries no code
chunk. Take the method sentences out of the results sections when you move them,
so the same sentence does not appear twice. A results section may still name the
cut it used where the reader needs it to follow the number.

Under `### Software`, name the packages the analysis used and say that every
version is pinned in the container. Do not invent a version number. If you want
to name a version and the file does not carry one, write the package name alone.

## Every report stands alone

Delete every reference to another report in this repository. The reader of one
report has not read the others.

These are the references that exist, and every one of them has to go.

| File | What to remove |
|---|---|
| `omip39_automated_gating.qmd` | "The PBMC report showed that..." |
| `yu2021_spectral_mait.qmd` | Four references to the OMIP-039 and OMIP-043 reports |
| `z282_harmonisation.qmd` | "The three earlier reports each ended on the same admission", "The earlier reports named a channel once", "the lesson of the OMIP-043 report", "the pattern every earlier report in this repository found", "The gap that three earlier reports named" |
| `oetjen2018_bone_marrow.qmd` | "the rule this repository uses in every other report", "the third dataset in this repository where an automated density cut fails" |
| `vanderbeke2021_covid.qmd` | "Every earlier report in this repository", "anywhere in this repository", "the least reliable cut in every report in this repository" |
| `flowcap2_challenges.qmd` | "Every other report in this repository", "every other gate in this repository" |

Removing a reference does not mean deleting the point it carried. When a
sentence says that a marker behaved the same way in another report, keep the
statement about the marker and drop the comparison. "CD45RA is a continuum, and
the OMIP-039 report met the same problem" becomes "CD45RA is a continuum".

You may still name a published paper, an accession, a package or a script in this
repository. What you may not do is point the reader at another document in
`reports/`.

## What you must not change

1. The YAML header, apart from the `title` and the `subtitle`. Every other line
   of it stays as it is.
2. Any ```{r} chunk. Not the code, not the label, not the chunk options. Do not
   move a chunk to another position in the file.
3. Any `caption =` string inside a chunk, and any `fig-cap:` line.
4. Any block quote that starts with `>`. Those are the words of the paper's
   authors. Never edit quoted material.
5. The `## Reproducing this` bash block.
6. The `## References` list.
7. Every number, every p value, every accession, every file name, every marker
   name, every package name.

You may add the `## Methods` header and its subheadings. You may not rename or
delete any header that is already there, and you may not reorder the sections
that hold code chunks.

## The numbers rule

You may not introduce a number that is not already in the file. You may not round
a number that is already there. If the file says the correlation is 0.972, you
write 0.972, not "above 0.97" and not "almost perfect".

## Rules that still bind you

Everything under "Rules for both modes" in `AGENTS.md`. No em dash and no en
dash. No bold phrase followed by a full stop. No label, colon and value inside
running text. No banned opener and no banned closer.

Every sentence is a complete sentence with a subject and a verb.

Vary sentence length. A methods section is plain, and plain is not uniform.

## Two things the first pass got wrong

Wrap every prose line at 80 characters. The rest of the repository does, and an
unwrapped paragraph produces a diff that cannot be reviewed line by line.

Leave no table and no figure without a sentence, and do not pad one either. A
sentence next to a table must add something the caption does not say. The
finding, the magnitude, or the reason the table is in the report all qualify. A
sentence that restates the caption is filler, and filler is worse than silence.

Most of these sentences already exist in the bullet or prose version of the file.
Where the original said something like "The largest value names the pair of
detectors that compensation actually moves. Inspect that pair on a plot", the
first sentence is content and the second is an instruction to the reader. Keep
the content and drop the instruction. Do not delete the whole pair and write a
new sentence that repeats the caption.

Read the previous revision with `git show HEAD:<path>` when you are unsure what
a table was there to show.

Recovering content this way can put a method sentence back into a results
section that you already moved into Methods. Check for that before you finish. A
short reminder of the method is allowed where the reader needs it to follow a
number. The same sentence appearing twice is not.

## The title names the study, it does not advertise it

Rewrite the `title` and the `subtitle` in the YAML header. Nothing else in the
YAML changes.

A title names the material and the question. It does not carry a verdict, a
promise or a hook. These are the current titles and the fault in each.

| Current title | Fault |
|---|---|
| "A viability channel that runs backwards, and what it costs an automated pipeline" | Gives away a finding as a hook, and "what it costs" promises drama. |
| "Thirteen operators, one panel, and what centralising the analysis actually fixes" | The word "actually" argues with an imagined sceptic. |
| "Entering a benchmark with one feature" | A story about the author rather than the study. |
| "A deposited workspace that gates every sample, and an automated rule that keeps up" | "Keeps up" is a verdict on the work. |
| "OMIP-39 as a positive control" | States the role the author assigned the dataset rather than what was measured. |
| "A 35 marker spectral panel, and a published finding tested against it" | The second half describes the exercise rather than the subject. |

Write the title as a paper would. Name the measurement, the material and the
accession where it fits. "Automated gating of a PBMC panel in R" is already
close to right and needs little.

The subtitle carries the material, the accession and the design. It does not
carry a claim about the result and it does not tell the reader what is
interesting.

Do not put a finding in the title. A finding belongs in the results.

The `title` and the `subtitle` each stay on one line, however long. The 80
character rule applies to prose and not to YAML, and a wrapped YAML value breaks
the header.

Never write an accession that is not already in the file. The first pass on the
OMIP-43 report invented `FR-FCM-Z2ZG` for its Methods section. The real
accession is `FR-FCM-ZYBP`, and the invented one points at nothing while reading
as a fact. When a Methods section needs an accession and the prose does not
carry one, take it from the `## Reproducing this` block, from a path under
`gating/`, or from the folder name under `data/datasets/flowrepository/`. If you
cannot find it, write the study name and leave the accession out.

## The introduction is not an argument for the dataset

Rename the `## Why this dataset` header to `## Introduction`. The question in
that header is the problem. A paper does not justify its choice of material to
the reader; it states what is known, what the material is, and what was tested.

The introduction says three things and stops.

1. What the source study reported, with its citation.
2. What the deposit contains.
3. What this report measured.

It does not say that the dataset was well chosen, that the data are clean, that
the claim can be tested, or what a reader should count as success.

These sentences are in the current files and every one of them has to go.

| Sentence | Why it goes |
|---|---|
| "The analysis reproduces the paper's findings as closely as the deposited data allows." | A claim about the quality of the work rather than a result. |
| "A report that reproduces only one of those two results has not reproduced the paper." | Tells the reader what counts as success. |
| "The benchmark answered its own question in two opposite ways." | Editorial framing of a result that the results section reports. |
| "The input is clean and the human boundaries are published." | An evaluation of the material. State what the files carry instead. |
| "The dataset was chosen because the claim can be tested." | Justifies the choice. |
| "That combination makes this dataset the fair test of the approach." | Same, with a superlative attached. |

## Do not editorialise inside a results section

A results section states what was measured and what the measurement was. Three
habits in the current text break that.

The first is defending the source paper. The OMIP-43 report currently writes
that the word "regularly" in a quoted sentence "describes the authors' usual
practice across many experiments rather than a promise about this deposit".
That is advocacy on the authors' behalf. Report the count against the stated
range, and if the wording matters, say that the sentence does not quantify
"regularly".

The second is teaching. "The two numbers in that sentence are linked by
arithmetic" explains to the reader that arithmetic is about to happen. Give the
arithmetic. A count of rare events follows a Poisson distribution, the relative
spread is `1 / sqrt(n)`, and that is 5 percent at 400 events.

The third is commentary on the analysis itself. Sentences that say a finding is
useful, expected, surprising, worth stating or worth recording are commentary. A
result that surprised the analyst can be reported as unexpected once, in the
discussion, where a paper puts it. It does not belong in the results.

## Tools

`rg` is not installed on this machine. Use `grep`. `grep -nP` is available for a
Unicode range, which is how you find an em dash.

## Work in passes on a long file

Reading the whole file and then planning the whole rewrite fails on the longer
reports. The session ends after the reads and it writes nothing.

Edit as you go. Take one section, rewrite it, call the edit tool, and move on. On
a file over 25 KB, split the work at the part headers.

## Acceptance test

Before you report that you are done, confirm all of these.

1. A `## Methods` section exists, it sits before the first result, and it holds no
   code chunk.
2. `grep -niE "this repository|earlier report|other report|PBMC report|OMIP-0?39 report|OMIP-0?43 report"` returns nothing outside a path such as `gating/` or `scripts/`.
3. No line contains an em dash or an en dash.
4. Every ```{r} chunk is byte identical to the original and sits in the same
   position.
5. The YAML header is byte identical to the original.
6. Every `>` block quote is byte identical to the original.
7. Every number in the prose appears in the original file.

Report which of the seven you checked and how.

The author runs `scripts/check_prose_rewrite.py` on your output afterwards. It
compares the file with its git revision and it fails on a changed chunk, a
changed block quote, an em dash or a number that is not in the original. It also
lists any header you added, which is expected for `## Methods` and its
subheadings and for nothing else.

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
outside the single file you were given. Do not commit.
