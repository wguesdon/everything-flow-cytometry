# Rewrite one report from bullet points into prose

You are rewriting the body text of one Quarto report in `reports/`. The reports
were written as bullet points on purpose, so that the prose could be written
later. You are writing that prose.

Read this whole brief before you change a line.

## The one thing that decides whether you succeeded

A reader who knows flow cytometry should finish a section without once thinking
about the writing. If a paragraph reads like a press release, a grant abstract or
a model that has been asked to sound confident, you have failed even if every
fact is correct.

## Which mode

Use Mode 2, the prose rules in `AGENTS.md`. This overrides the rule that sends
repository documentation to Simplified Technical English. The author chose Mode 2
for these files deliberately, on 2026-08-17.

Everything under "Rules for both modes" in `AGENTS.md` still binds you. Read it
again before you start. The rules that get broken most often here are the ban on
em dashes, the ban on a bold phrase followed by a full stop, and the ban on a
label, a colon and a value inside running text.

## What you may change

Only the bullet lists that carry body text. A bullet list becomes one or more
paragraphs.

## What you must not change

1. The YAML header. Not one character.
2. Any ```{r} chunk. Not the code, not the label, not the chunk options.
3. Any `caption =` string inside a chunk, and any `fig-cap:` line. Those render
   inline and they are already written.
4. Any block quote that starts with `>`. Those are the words of the paper's
   authors and they are quoted. Never edit quoted material.
5. Any section header.
6. The `## Reproducing this` bash block.
7. The `## References` list.
8. Every number, every p value, every accession, every file name, every marker
   name, every package name.

## The numbers rule

You may not introduce a number that is not already in the file. You may not round
a number that is already there. If a bullet says the correlation is 0.972, the
prose says 0.972 and not "above 0.97" and not "almost perfect".

If you find yourself wanting a number that is not in the file, write the sentence
without it.

## How to turn bullets into paragraphs

A run of bullets is usually one argument. Find the argument first, then write it.

Do not write one sentence per bullet and join them with a full stop. That is the
most common failure and it produces a paragraph with the rhythm of a list, which
is worse than the list.

Do the opposite. Decide what the run of bullets is claiming, write that claim in
a sentence that carries it, and let the supporting facts fall into the sentences
after it. Some bullets will merge into one clause. Some will earn a paragraph of
their own. A bullet that carried no argument should disappear.

## Rhythm

Vary sentence length on purpose, and get the variance from the whole paragraph
rather than from a clipped sentence at the end. A two word payoff after a long
setup is the machine pattern, not the cure for it.

Vary the skeleton of adjacent sentences. If one starts with the subject, the next
one can start with the condition.

Vary paragraph length. A paragraph of one sentence is allowed when that sentence
carries a result. A paragraph of six sentences is allowed when the argument needs
six.

A short sentence must carry a fact. "Smallest won." is a fact. "That matters." is
rhythm, and it is banned.

## Things that mark the text as machine written

Do not do any of these.

1. Open a paragraph with "Notably", "Importantly", "Interestingly", "Crucially",
   "It is worth noting" or "It should be noted".
2. Close a section with a sentence that restates the section.
3. Write "This is not just X, it is Y" or "It is not about X, it is about Y".
4. Use a colon in place of a verb.
5. Put a bold lead-in on a paragraph.
6. Force three items when the content has two or four.
7. Use an em dash or an en dash. A comma, a full stop or a semicolon.
8. Use a parenthetical aside. Fold it into the sentence or give it its own
   sentence.
9. Hedge a result that the data settles. If the p value is 0.0001, say the effect
   is there.
10. Announce the structure of your own text, as in "Three findings follow".

## Voice

The author is a computational biologist writing for another one. They are
allowed to say that something surprised them, that a first attempt failed, and
that a question is still open. Several of these reports already say exactly that
in bullet form, and the prose should keep the admission rather than smooth it
away.

Understate a result rather than sell it. The numbers are strong enough.

Write in the third person about the analysis and the first person plural only
where the bullets already do. Do not add a narrator who was not there.

## The order of work, for one file

1. Read the whole file first, including the code chunks, so you know what each
   table and figure shows.
2. Rewrite section by section, from the top.
3. After each section, read your paragraphs aloud in your head. If you would not
   say a line in that form, write it again.
4. When the file is done, check it against the list in "What you must not
   change".
5. Check every number in your prose against the bullet text you replaced.

## Tools

`rg` is not installed on this machine. Use `grep`. `grep -nP` is available for a
Unicode range, which is how you find an em dash.

## Work in passes on a long file

Reading the whole file and then planning the whole rewrite fails on the longer
reports. The session ends after the reads and it writes nothing.

Edit as you go instead. Take one section, rewrite it, call the edit tool, and
move to the next section. On a file over 25 KB, split the work at the part
headers and treat each part as its own pass.

## One session at a time

Two of these runs in parallel will cut one of them off without an error message.
The author runs them one after another.

## Acceptance test

Before you report that you are done, confirm all of these.

1. `grep -c '^- '` on the file returns zero, or returns only bullets that are a
   genuine list a reader would want.
2. No line contains an em dash or an en dash.
3. Every ```{r} chunk is byte identical to the original.
4. The YAML header is byte identical to the original.
5. Every `>` block quote is byte identical to the original.
6. Every number in the prose appears in the original file.
7. No word from the banned vocabulary list in `AGENTS.md` appears, unless the
   field uses it as a technical term with one exact meaning. A hit inside the
   title of a cited paper does not count, because you may not change a citation.

Report which of the seven you checked and how.

The author runs `scripts/check_prose_rewrite.py` on your output afterwards. It
compares the file with its git revision and it fails on a changed chunk, a
changed header, a changed block quote, an em dash, a leftover bullet or a number
that is not in the original. Your own check and that script should agree.

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

## Do not

Do not run the analysis scripts. Do not render the report. Do not touch anything
outside the single file you were given. Do not commit.
