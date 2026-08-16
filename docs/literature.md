# Literature

Every citation on this page was checked against a retrieved record on 2026-08-13. The sources were
Europe PMC, PubMed, Crossref, Semantic Scholar and the PMC full text. A claim that failed verification
is not on this page.

## Start here

Read these four first. They give the field, the method and the current state.

| Reference | Why you read it |
|---|---|
| Saeys Y, Van Gassen S, Lambrecht BN. Computational flow cytometry: helping to make sense of high-dimensional immunology data. Nat Rev Immunol 2016;16:449-62. PMID 27320317 | The best entry point to computational flow cytometry. |
| Liechti T, Weber LM, Ashhurst TM, Stanley N, Prlic M, Van Gassen S, Mair F. An updated guide for the perplexed: cytometry in the high-dimensional era. Nat Immunol 2021;22:1190-7. PMID 34489590 | A practical guide to high dimensional panels. |
| Mair F, et al. The end of gating? An introduction to automated analysis of high dimensional cytometry data. Eur J Immunol 2016;46:34-43. PMID 26548301 | An introduction to automated analysis. |
| O'Neill K, Aghaeepour N, Spidlen J, Brinkman R. Flow cytometry bioinformatics. PLoS Comput Biol 2013;9:e1003365. PMID 24363631 | The founding overview. It is open access. |

## The scaling problem

Robinson JP, Gmyrek GB, Rajwa B. Flow Cytometry: Advances, Challenges and Trends. BioEssays
2026;48(1):e70091. PMID 41311164. PMC12706142. The paper is open access under CC BY.

The paper gives the arithmetic for the number of bi-axial plots. The count is `(N x (N - 1)) / 2`,
where `N` is the number of channels.

| Panel | Possible bi-axial plots |
|---|---|
| 3 colours | 3 |
| 10 colours | 45 |
| 45 markers | 990 |

The same paper reports that spectral analysers with five lasers carry 51 to 144 detectors. Those
instruments measure 15 to 50 signals at one time.

## Gating causes more variation than the laboratory does

The studies below share one design. The same FCS files go to many analysts. Donor, reagent, staining
and instrument variation are then zero, and the remaining variation comes from the gates.

| Reference | Design | Result |
|---|---|---|
| Maecker HT, McCoy JP Jr, for the FOCIS Human Immunophenotyping Consortium. A model for harmonizing flow cytometry in clinical trials. Nat Immunol 2010;11(11):975-978. PMID 20959798 | Pre-stained cells to 15 experienced laboratories | Local gating gave a mean CV of 20.5%. Central gating of the same raw files gave 4%. |
| Maecker HT, et al. Standardization of cytokine flow cytometry assays. BMC Immunol 2005;6:13. PMID 15978127 | The primary source behind the 2010 paper | Central analysis with a dynamic gating template reduced the CVs to 3 to 7%. |
| Finak G, et al. Standardizing Flow Cytometry Immunophenotyping Analysis from the Human ImmunoPhenotyping Consortium. Sci Rep 2016;6:20686. PMID 26861911 | Nine sites, five standardised eight colour panels, one SOP | Site specific gating had a larger effect on assay sensitivity than centre to centre technical variability. |
| McNeil LK, et al. A harmonized approach to intracellular cytokine staining gating. Cytometry A 2013;83A:728-38. PMID 23788464 | 110 laboratories, one shared set of four colour files | The 110 laboratories used 110 different gating approaches. False positive calls fell from 23% to 9% when only the gating changed. |
| Westera L, et al. Clin Transl Gastroenterol 2017;8(11):e126. PMID 29095427 | Prospective multicentre study, three gating arms | Local gating gave a mean CV of 4.4 to 102.1%. Central gating gave 1.8 to 20.9%. |
| Gouttefangeas C, et al. Cancer Immunol Immunother 2015;64:585-98. PMID 25854580 | 17 laboratories re-gated shared multimer files | Data analysis is a source of variation in the multimer assay outcome. |
| Mandruzzato S, et al. Cancer Immunol Immunother 2016;65:161-9. PMID 26728481 | 23 laboratories, MDSC phenotyping | The gating strategy was the main parameter associated with variation. |
| Gratama JW, et al. Cytometry 1997;30:10-22. PMID 9056737 | The 1997 original | A uniform instrument setup reduced variability by 13%. Standard list mode analysis reduced it by 43%. |

Liu P, et al. Comprehensive evaluation and practical guideline of gating methods for high-dimensional
cytometry data. Brief Bioinform 2024;26:bbae633. PMID 39656848. Five expert raters gated one PBMC
dataset. The pairwise kappa index ran from 0.44 to 0.86.

## Automation compared with manual gating

Delgado AH, et al. Front Immunol 2023;14:1268686. PMID 37915569. The study compared EuroFlow automated
gating with manual gating on a 14 colour, 18 antibody panel.

| Measure | Automated | Manual |
|---|---|---|
| Median analysis time per sample | 6 minutes | 40 minutes |
| Inter-expert median CV | 3.9% | 17.3% |
| Intra-sample median CV | 1.7% | 10.4% |

The panel resolved 117 distinct B-cell and plasma-cell subsets.

Chen J, et al. Automated cytometric gating with human-level performance using bivariate segmentation.
Nat Commun 2025;16:1576. PMID 39939580. The UNITO method deviates from the human consensus by no more
than any single human does.

## How much automation the field actually uses

Cheung M, Campbell JJ, Whitby L, Thomas RJ, Braybrook J, Petzing J. Current trends in flow cytometry
automated data analysis software. Cytometry A 2021;99(10):1007-1021. PMID 33606354. The survey covered
51 software tools and 49 respondents from clinical laboratories.

| Measure | Value |
|---|---|
| Respondents who never use automated software | 26 of 49, which is 53% |
| Respondents who mainly use automated software | 1 of 49 |
| Tools available in R | 59% |
| Tools available in Matlab | 29% |
| Tools available in Python | 18% |
| Tools with a graphical interface | 41% |

Liu P, Liu S, Fang Y, Xue X, Zou J, Tseng G, Konnikova L. Recent Advances in Computer-Assisted
Algorithms for Cell Subtype Identification of Cytometry Data. Front Cell Dev Biol 2020;8:234.
PMID 32411698. The authors could run 21 of the 32 published methods. The other 11 were no longer
available or could not be installed.

## Data sharing

Leipold MD, Olsen LR. A literature study and public survey on mass cytometry dataset release and
reuse. Cytometry A 2022;101(2):109-13. PMID 34757690. The journal version is behind a paywall. The
authors posted a CC BY preprint with the full text at
<https://doi.org/10.6084/m9.figshare.15087294.v1>.

The study covered 692 CyTOF papers from 1 January 2018 to 17 June 2020.

| Measure | Value |
|---|---|
| Papers that generated new data | 563 |
| Papers with an accessible dataset | 24.5% |
| Papers with no link to the dataset | 72.1% |
| Public datasets released | 89 |
| Times those datasets were reused | 167 |
| Share of all reuse taken by two datasets | 19% |

The two datasets are Levine 2015 and Samusik 2016. Both were used many times to develop and to
benchmark analysis algorithms. In papers that produced both data types, sequencing data was deposited
at more than three times the rate of the mass cytometry data.

Other papers that make the same argument:

- Lin D, Gururaj A, Lin-Gibson S, Wang L. AI and flow cytometry. J Immunol 2026;215(2):vkaf292.
  PMID 41212078. A group that includes NIST reports that millions of flow cytometry datasets are
  siloed and cannot be used for AI applications.
- Hu Z, Bhattacharya S, Butte AJ. Application of Machine Learning for Cytometry Data. Front Immunol
  2021;12:787574. PMID 35046945. FlowRepository held 1,375 cytometry datasets in September 2021. GEO
  held 160,010 transcriptomics datasets.
- Spidlen J, Brinkman RR. Use FlowRepository to share your clinical data upon study publication.
  Cytometry B 2018. PMID 27342384.
- Lucas F, et al. MiSet RFC Standards. Cytometry A 2020;97(2):148-155. PMID 31769204.

## Standards and reporting

- Lee JA, et al. MIFlowCyt: the minimum information about a flow cytometry experiment. Cytometry A
  2008. PMID 18752282.
- Spidlen J, et al. FlowRepository: a resource of annotated flow cytometry datasets. Cytometry A 2012.
  PMID 22887982.
- Cossarizza A, et al. Guidelines for the use of flow cytometry and cell sorting in immunological
  studies, third edition. Eur J Immunol 2021. PMID 34910301.

Read the requirement carefully. Cytometry A makes MIFlowCyt reporting mandatory. It requests the
public release of the original list mode data. The two levels are not the same, and a cytometrist
knows the difference.

## Benchmarks

- Aghaeepour N, et al. Critical assessment of automated flow cytometry data analysis techniques
  (FlowCAP). Nat Methods 2013.
- Weber LM, Robinson MD. Comparison of clustering methods for high-dimensional single-cell flow and
  mass cytometry data. Cytometry A 2016;89:1084-96. PMID 27992111.
- Van Gassen S, et al. FlowSOM: Using self-organizing maps for visualization and interpretation of
  cytometry data. Cytometry A 2015. PMID 25573116.
