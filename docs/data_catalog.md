# Data catalogue

The `data/` folder is gitignored and lives in S3. This page is committed, so
you can choose what to pull before you transfer anything.

Generated on 2026-08-16 by `scripts/make_data_catalog.sh`.

## How to pull one folder

```bash
./sync.sh catalog                                  # the same table, read from S3
./sync.sh pull datasets/flowrepository/FR-FCM-ZZZU  # one accession
./sync.sh pull literature repositories             # two folders at once
./sync.sh pull                                     # everything, about 100 GB
```

A full pull transfers about 100 GB. Name the folders you need instead.

## Top level

| Folder | Size | Files | FCS files |
|---|---|---|---|
| `datasets` | 100G | 1967 | 1811 |
| `literature` | 1.9M | 1 | 0 |
| `repositories` | 2.2G | 1786 | 92 |

## datasets/

| Folder | Size | Files | FCS files |
|---|---|---|---|
| `datasets/FlowCal_data` | 4.0K | 0 | 0 |
| `datasets/flowcore_data` | 5.5M | 1 | 0 |
| `datasets/flowjo` | 207M | 10 | 8 |
| `datasets/FlowKit_data` | 4.0K | 0 | 0 |
| `datasets/flowrepository` | 100G | 1951 | 1799 |
| `datasets/NIH_ImmPort` | 4.0K | 0 | 0 |
| `datasets/readfcs_data` | 81M | 4 | 4 |

## datasets/flowrepository/

These are FlowRepository downloads. The source site stopped accepting new
experiments in 2025 and its TLS certificate expired on 18 March 2023, so
several of these accessions are hard to download again. Treat this copy as
the working copy.

| Folder | Size | Files | FCS files |
|---|---|---|---|
| `datasets/flowrepository/FlowRepository_FR-FCM-Z244_files` | 1.3G | 28 | 28 |
| `datasets/flowrepository/FlowRepository_FR-FCM-Z3WR_files` | 3.8G | 83 | 83 |
| `datasets/flowrepository/FlowRepository_FR-FCM-Z4KT_files` | 518M | 18 | 16 |
| `datasets/flowrepository/FlowRepository_FR-FCM-ZZLV_files` | 50M | 11 | 3 |
| `datasets/flowrepository/FlowRepository_FR-FCM-ZZZV_files` | 2.1G | 241 | 240 |
| `datasets/flowrepository/FR-FCM-Z282` | 5.8G | 254 | 250 |
| `datasets/flowrepository/FR-FCM-Z2KP` | 872M | 50 | 49 |
| `datasets/flowrepository/FR-FCM-Z32U` | 1.3G | 2 | 0 |
| `datasets/flowrepository/FR-FCM-Z6UG` | 38M | 9 | 8 |
| `datasets/flowrepository/FR-FCM-ZYQ9` | 14G | 141 | 132 |
| `datasets/flowrepository/FR-FCM-ZYQB` | 3.1G | 8 | 8 |
| `datasets/flowrepository/FR-FCM-ZYRN` | 1.4G | 61 | 61 |
| `datasets/flowrepository/FR-FCM-ZZCA` | 5.2M | 8 | 5 |
| `datasets/flowrepository/FR-FCM-ZZZU` | 5.6G | 309 | 308 |
| `datasets/flowrepository/FR-FCM-ZZZV` | 557M | 61 | 60 |
| `datasets/flowrepository/OMIP-018` | 292M | 20 | 17 |
| `datasets/flowrepository/OMIP-030` | 128M | 2 | 0 |
| `datasets/flowrepository/OMIP-058` | 1.2G | 4 | 0 |
| `datasets/flowrepository/OMIP-16` | 407M | 5 | 0 |
| `datasets/flowrepository/OMIP-23` | 45M | 2 | 0 |
| `datasets/flowrepository/OMIP-24` | 276M | 2 | 0 |
| `datasets/flowrepository/OMIP-39` | 343M | 16 | 13 |
| `datasets/flowrepository/OMIP-40` | 103M | 13 | 9 |
| `datasets/flowrepository/OMIP-43` | 4.2G | 243 | 233 |
| `datasets/flowrepository/OMIP-44` | 5.4G | 39 | 36 |
| `datasets/flowrepository/OMIP-47` | 1.3G | 9 | 7 |
| `datasets/flowrepository/OMIP-51` | 2.0G | 62 | 60 |
| `datasets/flowrepository/OMIP-60` | 275M | 35 | 33 |
| `datasets/flowrepository/OMIP-80` | 4.8G | 9 | 0 |
| `datasets/flowrepository/Pytometry` | 12G | 145 | 140 |
| `datasets/flowrepository/Spectral_Flow_Workflow-main` | 2.1M | 4 | 0 |

## Total

| Measure | Value |
|---|---|
| Size of `data/` | 103G |
| Files | 3755 |
| FCS files | 1903 |
