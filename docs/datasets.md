# Public datasets

This page lists the places that hold public FCS files. Every count was checked on 2026-08-13.

Start with HDCytoData if you want a dataset today. Start with ImmPort if you want a durable link for a
paper.

## Where to get data now

| Repository | Funding | Content | State |
|---|---|---|---|
| HDCytoData | Bioconductor | Ten standard benchmark datasets on ExperimentHub | The fastest route. The ten datasets came from FlowRepository. |
| ImmPort | NIAID contract | 1,502 studies at release DR66. 330 are tagged Flow Cytometry and 86 are CyTOF. | The durable option. The API is open and needs no authentication. |
| Zenodo | CERN | 304 records hold at least one FCS file | It has no cytometry metadata model. Treat it as a store, not as an archive. |
| Dryad | Non-profit | 174 datasets match the exact phrase "flow cytometry" | Small but stable. |
| NanoFlow Repository | Academic | 279 datasets, for extracellular vesicles | It is still marked BETA. Its own API states that the IRIs are not permanent links. Do not cite it as an archival home. |
| HuBMAP | NIH | 23 CyTOF and 16 imaging mass cytometry datasets | Small. |

The Zenodo curve is steep. The record counts by year are 26 in 2022, 31 in 2023, 59 in 2024, 98 in
2025 and 62 by 13 August 2026.

## FlowRepository

FlowRepository is the repository that the community guidelines point to. It has three problems.

1. It stopped accepting new experiments. A site-wide banner states that the creation of new
   experiments is disabled because the service runs out of space. Wayback snapshots place the banner
   between 15 April and 10 May 2025. It is in every snapshot since that date.
2. The R client is gone. Bioconductor removed `FlowRepositoryR` at release 3.14 in October 2021. No
   supported programmatic client exists.
3. The TLS certificate expired on 18 March 2023. Browsers and default HTTP clients refuse the
   connection. The site is reachable over plain HTTP or with certificate verification off.

The funding page names the Wallace H. Coulter Foundation and the International Society for Advancement
of Cytometry. It names no government source. GEO runs on NCBI. ImmPort runs on a NIAID contract. That
difference explains why a disk shortfall closes submissions at one service and not at the others.

## Scale

| Repository | Datasets | Date |
|---|---|---|
| GEO | 293,250 series | 2026-08-13 |
| FlowRepository | 1,375 datasets | September 2021, the last published count |

## Where the data does not go

The Allen Institute Human Immune Health Atlas, published in Nature in 2025, profiled more than 300
adults with flow cytometry. Its data availability statement names GEO for the scRNA-seq data, dbGaP for
the raw fastq files, and HISE and Zenodo for the code. It names no repository for the flow cytometry
data.

## To do

1. Choose one public dataset for the worked example in `examples/`.
2. Record the accession, the panel and the licence for that dataset.
3. Write a download script so the example starts from the accession and not from a local file.
