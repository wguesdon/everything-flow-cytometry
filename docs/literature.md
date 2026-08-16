# Literature

Every citation on this page was checked against a retrieved record. The sources were Europe PMC,
PubMed, Crossref, Semantic Scholar and the PMC full text. A claim that failed verification is not on
this page.

| Section | Date of the check |
|---|---|
| Every section down to "Benchmarks" | 2026-08-13 |
| "Papers behind the packages" | 2026-08-16 |

`docs/packages.md` lists the software. The last section of this page gives the paper for each package
in that list. Read `scripts/verify_package_papers.sh` to repeat the check.

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
  cytometry data. Cytometry A 2015. PMID 25573116. The full record is in the next section, because
  FlowSOM is also a package in `docs/packages.md`.

## Papers behind the packages

Cite the paper below when you use the package in `docs/packages.md`. Every record here was retrieved
from Europe PMC on 2026-08-16. `scripts/verify_package_papers.sh` runs the same queries again.

The `Open access` column states what the Europe PMC record reports. A paper marked `No` can still
have a free author copy, so check the publisher page before you pay for it.

### R packages

| Package | Reference | Open access |
|---|---|---|
| flowCore | Hahne F, LeMeur N, Brinkman RR, Ellis B, Haaland P, Sarkar D, Spidlen J, Strain E, Gentleman R. flowCore: a Bioconductor package for high throughput flow cytometry. BMC Bioinformatics 2009;10:106. PMID 19358741. PMC2684747. doi 10.1186/1471-2105-10-106 | Yes |
| CytoML | Finak G, Jiang W, Gottardo R. CytoML for cross-platform cytometry data sharing. Cytometry A 2018;93(12):1189-1196. PMID 30551257. PMC6443375 | Yes |
| ggcyto | Van P, Jiang W, Gottardo R, Finak G. ggCyto: next generation open-source visualization software for cytometry. Bioinformatics 2018;34(22):3951-3953. PMID 29868771. PMC6223365 | Yes |
| PeacoQC | Emmaneel A, Quintelier K, Sichien D, Rybakowska P, Marañón C, Alarcón-Riquelme ME, Van Isterdael G, Van Gassen S, Saeys Y. PeacoQC: Peak-based selection of high quality cytometry data. Cytometry A 2022;101(4):325-338. PMID 34549881. PMC9293479 | Yes |
| flowAI | Monaco G, Chen H, Poidinger M, Chen J, de Magalhães JP, Larbi A. flowAI: automatic and interactive anomaly discerning tools for flow cytometry data. Bioinformatics 2016;32(16):2473-2480. PMID 27153628 | No |
| CytoNorm | Van Gassen S, Gaudilliere B, Angst MS, Saeys Y, Aghaeepour N. CytoNorm: A Normalization Algorithm for Cytometry Data. Cytometry A 2020;97(3):268-278. PMID 31633883. PMC7078957 | Yes |
| CytoNorm 2.0 | Quintelier KLA, Willemsen M, Bosteels V, Aerts JGJV, Saeys Y, Van Gassen S. CytoNorm 2.0: A flexible normalization framework for cytometry data without requiring dedicated controls. Cytometry A 2025;107(2):69-87. PMID 39871681 | No |
| openCyto | Finak G, Frelinger J, Jiang W, Newell EW, Ramey J, Davis MM, Kalams SA, De Rosa SC, Gottardo R. OpenCyto: an open source infrastructure for scalable, robust, reproducible, and automated, end-to-end flow cytometry data analysis. PLoS Comput Biol 2014;10(8):e1003806. PMID 25167361. PMC4148203 | Yes |
| flowDensity | Malek M, Taghiyar MJ, Chong L, Finak G, Gottardo R, Brinkman RR. flowDensity: reproducing manual gating of flow cytometry data by automated density-based cell population identification. Bioinformatics 2015;31(4):606-607. PMID 25378466. PMC4325545 | No, but PMC holds the full text |
| flowClust | Lo K, Hahne F, Brinkman RR, Gottardo R. flowClust: a Bioconductor package for automated gating of flow cytometry data. BMC Bioinformatics 2009;10:145. PMID 19442304. PMC2701419 | Yes |
| FlowSOM | Van Gassen S, Callebaut B, Van Helden MJ, Lambrecht BN, Demeester P, Dhaene T, Saeys Y. FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data. Cytometry A 2015;87(7):636-645. PMID 25573116. doi 10.1002/cyto.a.22625 | No |
| diffcyt | Weber LM, Nowicka M, Soneson C, Robinson MD. diffcyt: Differential discovery in high-dimensional cytometry via high-resolution clustering. Commun Biol 2019;2:183. PMID 31098416. PMC6517415 | Yes |
| CATALYST | Chevrier S, Crowell HL, Zanotelli VRT, Engler S, Robinson MD, Bodenmiller B. Compensation of Signal Spillover in Suspension and Imaging Mass Cytometry. Cell Syst 2018;6(5):612-620.e5. PMID 29605184. PMC5981006 | Yes |
| CATALYST | Crowell HL, Chevrier S, Jacobs A, Sivapatham S, Tumor Profiler Consortium, Bodenmiller B, Robinson MD. An R-based reproducible and user-friendly preprocessing pipeline for CyTOF data. F1000Res 2020;9:1263. PMID 36072920. PMC9411975. doi 10.12688/f1000research.26073.2 | Yes |
| cydar | Lun ATL, Richard AC, Marioni JC. Testing for differential abundance in mass cytometry data. Nat Methods 2017;14(7):707-709. PMID 28504682. PMC6155493 | Yes |
| HDCytoData | Weber LM, Soneson C. HDCytoData: Collection of high-dimensional cytometry benchmark datasets in Bioconductor object formats. F1000Res 2019;8:1459. PMID 31857895. PMC6904983 | Yes |

CATALYST has two papers. Cite Chevrier 2018 for the spillover compensation method. Cite Crowell 2020
for the preprocessing pipeline and the differential analysis workflow.

CytoNorm also has two papers. Cite Van Gassen 2020 when your batches carry a dedicated control sample.
Cite Quintelier 2025 when they do not.

### Python packages

| Package | Reference | Open access |
|---|---|---|
| FlowKit | White S, Quinn J, Enzor J, Staats J, Mosier SM, Almarode J, Denny TN, Weinhold KJ, Ferrari G, Chan C. FlowKit: A Python Toolkit for Integrated Manual and Automated Cytometry Analysis Workflows. Front Immunol 2021;12:768541. PMID 34804056. PMC8602902 | Yes |
| FlowSOM (Python) | Couckuyt A, Rombaut B, Saeys Y, Van Gassen S. Efficient cytometry analysis with FlowSOM in Python boosts interoperability with other single-cell tools. Bioinformatics 2024;40(4):btae179. PMID 38632080. PMC11052654 | Yes |
| Pytometry | Büttner M, Hempel F, Ryborz T, Theis FJ, Schultze JL. Pytometry: Flow and mass cytometry analytics in Python. bioRxiv 2022. doi 10.1101/2022.10.10.511546. Europe PMC id PPR557097 | Yes, as a preprint |
| cytoflow | Teague L. Cytoflow: User-Friendly Python Software for Computational Flow Cytometry. Cytometry A 2026. PMID 42596036. doi 10.1002/cyto.a.70058 | No |
| scanpy | Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol 2018;19(1):15. PMID 29409532. PMC5802054 | Yes |

Two entries need care. Pytometry has a bioRxiv preprint and no journal article, and a search of Europe
PMC and Crossref on 2026-08-16 found no published version. The cytoflow record carries no volume and
no page numbers yet, because the journal published it in 2026. An earlier cytoflow preprint exists at
doi 10.1101/2022.07.22.501078, and its author string differs from the journal record.

### Packages with no paper

A search of Europe PMC on 2026-08-16 returned no article that describes these packages. Cite the
package and its version, or cite the paper in the third column.

| Package | Language | Cite this instead |
|---|---|---|
| flowWorkspace | R | Finak 2018 for CytoML, or Finak 2014 for openCyto. Both use the GatingSet object that flowWorkspace defines. |
| flowStats | R | Nothing. Cite the Bioconductor package and its version. |
| readfcs | Python | Nothing. Cite the PyPI package and its version. |
| FlowRepositoryR | R | Spidlen J, et al. FlowRepository: a resource of annotated flow cytometry datasets. Cytometry A 2012. PMID 22887982. The package is removed and the repository paper stays valid. |
