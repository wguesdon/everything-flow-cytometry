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

1. The YAML header. Not one character.
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

## Do not

Do not run the analysis scripts. Do not render the report. Do not touch anything
outside the single file you were given. Do not commit.
