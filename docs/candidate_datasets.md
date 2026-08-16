# Candidate datasets for the next report

This page lists every dataset in `data/datasets/flowrepository/` that carries a published paper. It
says what each deposit can test and what it cannot. Four reports are complete. This page chooses the
fifth.

Every accession was read from the local folder on 2026-08-16. Every citation was retrieved from
Europe PMC on the same date. No citation is written from memory.

## The gap the four reports left

Each report names a check that its dataset cannot support. The same three items repeat.

| Gap | Named in | What it needs |
|---|---|---|
| Analyst to analyst variation | `automated_gating_pbmc`, `omip39_automated_gating`, `omip43_asc_analysis` | The same files gated by many people |
| A subject random effect in place of per sample tests | `yu2021_spectral_mait` | Repeated measures with a known design |
| Batch correction across acquisition dates | `yu2021_spectral_mait` | Several instruments or several dates, with a shared control |
| A second measurement technology | none of the four | Two technologies on the same donors |

## The shortlist

### 1. FR-FCM-Z282, the thirteen operator harmonisation study

Macchia I, La Sorsa V, Ruspantini I, Sanchez M, Tirelli V, Carollo M, Fedele G, Leone P, Schiavoni G,
Buccione C, Rizza P, Nisticò P, Palermo B, Morrone S, Stabile H, Rughetti A, Nuti M, Zizzari IG,
Fionda C, Maggio R, Capuano C, Quintarelli C, Sinibaldi M, Agrati C, Casetti R, Rozo Gonzalez A,
Iacobone F, Gismondi A, Belardelli F, Biffoni M, Urbani F. Multicentre Harmonisation of a Six-Colour
Flow Cytometry Panel for Naïve/Memory T Cell Immunomonitoring. J Immunol Res 2020;2020:1938704. PMID
32322591. PMC7153001. doi:10.1155/2020/1938704.

| Item | Value |
|---|---|
| Files | 248 sample FCS, plus 2 more in `attachments/`, 5.8 GB, extracted |
| Design | 13 operators, 5 centres, 7 instruments, 3 rounds, 3 vials, PBMC and whole blood |
| Panel | CD3, CD4, CD8, CD45RA, CCR7 and a viability dye |
| Instruments | FACSCanto, FACSCantoII, Gallios, LSRFortessa |
| Manual gates | `attachments/representative_flowjo_analysis.wsp` |
| Result table | `attachments/De_Identified_Proficiency_naive_mem_DataBase.xlsx`, 468 rows, 15 frequencies per row |

The file name encodes the operator, the instrument, the centre, the material and the vial. An example
is `OpA_I2_C2_IM1_PBMC1_R1.fcs`.

The result table holds 234 local rows and 234 central rows. A local row is the operator's own analysis
of the operator's own file. A central row is one analyst at the coordinating centre who analysed the
same file. The design is balanced. Each of the 13 operators contributes 36 rows. Twelve rows carry an
exclusion flag.

This deposit answers the first gap directly. It is the only dataset in the archive where the same
files were gated by many people and both answers were published.

It also carries a second problem that no other deposit has. The 248 files come from four instrument
models, the parameter count runs from 10 to 25, and the channel names differ between the makers. Some
files carry no marker name at all. One template that gates all 248 files must resolve a channel by the
fluorochrome and not by the detector name.

What the deposit cannot do: the operators stained and acquired their own tubes, so a local and a
central value differ by the analysis alone, and two operators differ by the staining, the instrument
and the analysis together.

### 2. FR-FCM-ZYQ9 with FR-FCM-ZYQB, flow cytometry against mass cytometry

Oetjen KA, Lindblad KE, Goswami M, Gui G, Dagur PK, Lai C, Dillon LW, McCoy JP, Hourigan CS. Human
bone marrow assessment by single-cell RNA sequencing, mass cytometry, and flow cytometry. JCI Insight
2018;3(23):e124928. PMID 30518681. PMC6328018. doi:10.1172/jci.insight.124928.

| Item | Value |
|---|---|
| Flow files | 132 FCS in FR-FCM-ZYQ9, 14 GB, extracted |
| Mass cytometry files | 8 FCS in FR-FCM-ZYQB, Helios, 49 named markers, extracted |
| Donors | 20 healthy adults, 10 male and 10 female, age 24 to 84, median 57 |
| Panels | B cell, T cell, NK, monocyte, dendritic cell, plus unstained controls |
| Manual gates | `workspace.xml` holds one sample only. `attachments/FR-FCM-ZYQ9_Gating_Hierarchy.pdf` holds the strategy as a figure |
| Cohort table | Table 1 of the paper gives the sex and the age of each donor |

The paper makes three claims from the flow data. The memory T cell frequency rises with age. Mass
cytometry correlates with flow cytometry on the T cell subsets of 8 donors. Flow cytometry and single
cell RNA sequencing disagree on the T cell and the NK cell frequency.

The first two claims are testable from this archive. The third needs the sequencing data, which sits
in GEO and not in `data/`.

The eight mass cytometry files match the eight donor subset that the paper compares. This is the only
pair of technologies on the same donors in the archive, and no report has read a mass cytometry file
yet.

What the deposit cannot do: there is no full manual gating workspace, so the manual against automated
comparison of the OMIP-039 report cannot repeat here.

### 3. FR-FCM-Z2KP, a 265 population workspace on COVID-19 samples

Vanderbeke L, Van Mol P, Van Herck Y, De Smet F, Humblet-Baron S, Martinod K, Antoranz A, Arijs I,
Boeckx B, Bosisio FM, Casaer M, Dauwe D, De Wever W, Dooms C, Dreesen E, Emmaneel A, Filtjens J, et
al. Monocyte-driven atypical cytokine storm and aberrant neutrophil activation as key mediators of
COVID-19 disease severity. Nat Commun 2021;12(1):4117. PMID 34226537. PMC8257697.
doi:10.1038/s41467-021-24360-w.

| Item | Value |
|---|---|
| Files | 49 FCS, 872 MB, extracted |
| Groups | 6 healthy, 23 mild to moderate, 20 severe |
| Panel | Intracellular cytokine staining. FOXP3, GATA3, Tbet, IL-2, IL-17a, IFNg, CD3, CD4, CD45RA and more |
| State | Compensated and pre-gated to live cells |
| Manual gates | `attachments/01-May-2020_Human_COVID_analysis_template.wsp`, 265 populations across 52 samples |

The accession is named in the data availability statement of the paper, which is how the paper was
identified. FlowRepository lists no manuscript for this record.

This is the richest workspace in the archive after OMIP-043. It is the natural dataset for a report on
transcription factor gating and cytokine gating, which no report has covered.

What the deposit cannot do: the files are already pre-gated to live cells, so the early quality
control steps cannot be shown.

### 4. FR-FCM-ZZZU and FR-FCM-ZZZV, the FlowCAP-II benchmarks

Aghaeepour N, Finak G, FlowCAP Consortium, DREAM Consortium, Hoos H, Mosmann TR, Brinkman R, Gottardo
R, Scheuermann RH. Critical assessment of automated flow cytometry data analysis techniques. Nat
Methods 2013;10(3):228-238. PMID 23396282. PMC3906045. doi:10.1038/nmeth.2365.

| Accession | Content |
|---|---|
| FR-FCM-ZZZU | 308 FCS. 44 infants, 7 stimulation conditions, 140 HIV exposed uninfected against 168 unexposed |
| FR-FCM-ZZZV | 240 FCS. Intracellular cytokine staining from trial HVTN 049, 120 training and 120 testing files |

Both deposits carry a metadata CSV with the sample class. `Challenge3Metadata.csv` also carries the
training and testing split that the challenge used.

These are the only datasets in the archive where the paper reports the performance of an algorithm and
not a biological frequency. A report here compares our pipeline against a published ranking of
methods, which is a different question from the four reports so far.

The HVTN 049 trial paper is Spearman P, Lally MA, Elizaga M, Montefiori D, Tomaras GD, McElrath MJ,
Hural J, De Rosa SC, et al. A trimeric, V2-deleted HIV-1 envelope glycoprotein vaccine elicits potent
neutralizing antibodies but limited breadth of neutralization in human volunteers. J Infect Dis
2011;203(8):1165-1173. PMID 21451004. PMC3068023. doi:10.1093/infdis/jiq175.

### 5. FR-FCM-Z244, a mass cytometry clinical trial panel

Hartmann FJ, Babdor J, Gherardini PF, Amir ED, Jones K, Sahaf B, Marquez DM, Krutzik P, O'Donnell E,
Sigal N, Chang HY, Rosenberg-Hasson Y, Mach KE, Liao YJ, Rieger K, Miklos DB, Maecker HT, Bendall SC.
Comprehensive Immune Monitoring of Clinical Trials to Advance Human Immunotherapy. Cell Rep
2019;28(3):819-831.e4. PMID 31315057. PMC6656694. doi:10.1016/j.celrep.2019.06.049.

| Item | Value |
|---|---|
| Files | 28 FCS, 1.3 GB, extracted, day 30 and day 90 after transplantation |
| Claim | Patients with graft versus host disease carry fewer CD27 negative B cells and fewer naive CD4 T cells |
| Attachments | None. There is no workspace and no clinical table in the deposit |

The claim needs the disease status of each patient, and the deposit does not carry it. The status must
be read from the paper before this dataset can be used.

## The recommendation

Take FR-FCM-Z282.

It closes the gap that three of the four reports name, it holds a published local result and a
published central result for every file, and its four instrument models make the channel resolution
problem real rather than hypothetical. The report writes itself as three arms.

1. The local analysis, which is 13 operators on their own files. The values are in the xlsx.
2. The central analysis, which is one analyst on every file. The values are in the same xlsx.
3. An openCyto template, which is this repository on every file. The values do not exist yet.

The paper claims that the central analysis reduces the cross-centre variability, and that the
reduction is larger in whole blood. Arm 1 and arm 2 test that claim from the deposit alone. Arm 3 asks
whether a template moves the coefficient of variation below the central analyst, which is the question
that the Maecker study raised and that no report in this repository has been able to ask.

Take FR-FCM-ZYQ9 with FR-FCM-ZYQB second. It brings mass cytometry into the repository and it tests an
age effect across 20 donors.

## Every accession in the archive

The `Used` column marks the four datasets that a report already covers.

| Folder | Accession | FCS | Paper | State |
|---|---|---|---|---|
| `FR-FCM-Z282` | FR-FCM-Z282 | 250 | PMID 32322591 | Extracted. Workspace and result table |
| `FR-FCM-ZYQ9` | FR-FCM-ZYQ9 | 132 | PMID 30518681 | Extracted. Panel and hierarchy PDFs |
| `FR-FCM-ZYQB` | FR-FCM-ZYQB | 8 | PMID 30518681 | Extracted. Helios mass cytometry |
| `FR-FCM-ZZZU` | FR-FCM-ZZZU | 308 | PMID 23396282 | Extracted. Class labels in `HEUvsUE.csv` |
| `FlowRepository_FR-FCM-ZZZV_files` | FR-FCM-ZZZV | 240 | PMID 23396282, PMID 21451004 | Extracted. Labels in `Challenge3Metadata.csv` |
| `FR-FCM-ZZZV` | FR-FCM-ZZZV | 60 | as above | A partial second copy of the same accession |
| `FR-FCM-Z2KP` | FR-FCM-Z2KP | 49 | PMID 34226537 | Extracted. Workspace of 265 populations |
| `FlowRepository_FR-FCM-Z244_files` | FR-FCM-Z244 | 28 | PMID 31315057 | Extracted. No attachment |
| `FlowRepository_FR-FCM-Z3WR_files` | FR-FCM-Z3WR | 83 | PMID 33870241 | Used by `yu2021_spectral_mait` |
| `FlowRepository_FR-FCM-Z4KT_files` | FR-FCM-Z4KT | 16 | PMID 34868024 | Extracted. The published R workflow is in `attachments/script.Rmd` |
| `OMIP-39` | FR-FCM-ZYY6 | 13 | PMID 28715616 | Used by `omip39_automated_gating` |
| `OMIP-43` | FR-FCM-ZYBP | 233 | PMID 29286577 | Used by `omip43_asc_analysis` |
| `OMIP-018` | FR-FCM-ZZ36 | 17 | PMID 23504907 | Extracted. The workspace is a legacy binary FlowJo file and it does not parse |
| `OMIP-40` | FR-FCM-ZY6D | 9 | PMID 28834328 | Extracted. The gates are a Cytobank XML with 1877 gate elements |
| `OMIP-44` | FR-FCM-ZYC2 | 36 | PMID 29356334 | Extracted. No workspace |
| `OMIP-47` | FR-FCM-ZYFB | 7 | PMID 29782066 | Extracted. No workspace |
| `OMIP-51` | FR-FCM-ZYN4 | 60 | PMID 30549419 | Extracted. No workspace |
| `OMIP-60` | FR-FCM-ZYRX | 33 | PMID 31334913 | Extracted. No workspace |
| `FR-FCM-ZYRN` | FR-FCM-ZYRN | 61 | PMID 31334918 | Extracted. This is the OMIP-058 deposit and it holds compensation beads |
| `FR-FCM-Z6UG` | FR-FCM-Z6UG | 8 | none listed | Extracted. Mouse depletion check. A workspace is present |
| `FR-FCM-ZZCA` | FR-FCM-ZZCA | 5 | see the PDF in the folder | Extracted. Five files only |
| `FlowRepository_FR-FCM-ZZLV_files` | FR-FCM-ZZLV | 3 | none listed | Extracted. Quality control scripts in `attachments/` |
| `Pytometry` | FR-FCM-ZYQ9, FR-FCM-ZYQB | 140 | PMID 30518681 | A second copy of the two bone marrow accessions |
| `OMIP-16` | FR-FCM-ZZ2T, ZZ2V, ZZ3Y, ZZ3Z | 0 | PMID 23184609 | Still in four zip files |
| `OMIP-23` | FR-FCM-ZZ74 | 0 | PMID 25132115 | Still in a zip file |
| `OMIP-24` | FR-FCM-ZZEB | 0 | PMID 25352070 | Still in a zip file |
| `OMIP-030` | FR-FCM-ZZWU | 0 | PMID 26506224 | Still in a zip file |
| `OMIP-058` | FR-FCM-ZYRN | 0 | PMID 31334918 | Three zip files. The folder `FR-FCM-ZYRN` holds the extracted copy |
| `OMIP-80` | FR-FCM-Z3JB | 0 | PMID 34693626 | Seven zip files, 4.8 GB |
| `FR-FCM-Z32U` | FR-FCM-Z32U | 0 | PMID 33200508 | Still in a zip file, 1.3 GB |
| `Spectral_Flow_Workflow-main` | none | 0 | PMID 34868024 | The R workflow of the FR-FCM-Z4KT paper, with no data |

## What to unzip first

Seven folders still hold their data in a zip file. Unzip a folder only when a report needs it, because
`data/` is already about 103 GB.

| Folder | Reason to unzip |
|---|---|
| `FR-FCM-Z32U` | Mair and Liechti published a full phenotyping study of dendritic cells and monocytes. It pairs with OMIP-044 by the same first author |
| `OMIP-80` | A 29 colour panel for NK and T cells. It is the largest unopened deposit |
| `OMIP-030` | A small T cell subset panel. It is the cheapest one to open |
