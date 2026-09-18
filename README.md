<div align="center">

# 🧬 NoCoSMiCC
### Non-Coding Somatic Mutations in Colorectal Cancer

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)](link-to-publication)
[![Genomics](https://img.shields.io/badge/Field-Cancer%20Genomics-purple.svg)]()

*Identifying non-coding driver mutations in colorectal cancer using 5,000+ whole genomes and colon epithelial-specific regulatory elements*

[Publication](#) • [Documentation](#) • [Data Access](#data-access)

</div>

---

## 🔬 Overview

NoCoSMiCC aims to identify non-coding somatic driver mutations in colorectal cancer (CRC) by integrating multi-modal regulatory genomics data with large-scale whole genome sequencing [web:1][web:26]. This pipeline generates candidate cis-regulatory elements (cCREs) in colon epithelial cells and facilitates the integration of somatic mutations 
**Key Features:**
- 🎯 Tissue-specific regulatory element identification from scATAC-seq, and histone/TF ChIP-seq
- 🧪 Integration of data from ENCODE and HuBMAP
- Generation of a candidate cis-regulatory element map in colon epithelial cells
- 📊 Evolutionary conservation and synteny analysis across 241 mammals
- 🏗️ TF motif disruption and chromatin contacts
- 🔍 Analysis of somatic mutations in >5,000 CRC whole genomes across five independent cohorts

## 📊 Citation

If you use NoCoSMiCC in your research, please cite:

@article{nocosmcc2026,
title={Non-Coding Somatic Mutations in Colorectal Cancer},
author={Hagan et al.},
journal={Journal Name},
year={2025},
doi={pending}
}

## 🚀 Quick Start
Clone the repository
git clone https://github.com/yourusername/NoCoSMiCC.git
cd NoCoSMiCC

Install dependencies
conda env create -f environment.yml
conda activate nocosmicc

Run the cCRE generation pipeline
./scripts/run_pipeline.sh

### Dependencies

## Data access

A key component of this project is the construction of colon epithelium-specific 
regulatory element maps, adapted from the ENCODE cCRE generation pipeline 
(https://www.encodeproject.org/), using single-cell chromatin accessibility data 
in addition to histone and transcription factor ChIP-seq.

**Genome build:** [FILL IN: e.g., GRCh38/hg38]

**Reference files required:**
- Genome FASTA: [FILL IN: source/version]
- TSS annotation (GTF/BED): [FILL IN: source/version, e.g. GENCODE vX]
- Blacklist/low-mappability regions: [FILL IN, if used]

The following data types are used to generate colon cCREs:
- scATAC-seq
- H3K4me3 ChIP-seq
- H3K27ac ChIP-seq
- CTCF ChIP-seq

## Running the cCRE generation pipeline

### Step 1: Generate peaks from scATAC-seq data
Fragment files are obtained for each sample and used to create arrow files for an 
ArchR project. Quality control removes low-quality cells (e.g., doublets, cells 
with low TSS enrichment and/or low fragment counts). Accessibility around marker 
genes is used to select epithelial cell populations for peak calling; a peak file 
and bigWig file are exported per sample.

```bash
./1_Process-sc-fragments.R
```

### Step 2: Process single-cell chromatin accessibility peaks
Peaks from the ArchR analysis are filtered and saved as chromatin accessible 
regions (CARs) for each sample.

```bash
./2_Process-sc-CARs.sh
```

### Step 3: Download and process bulk chromatin accessibility peaks
Bulk DNase and ATAC peaks are downloaded from ENCODE, filtered, and saved as CARs 
for each sample.

```bash
./3_Obtain-Bulk-Peaks.sh
```

### Step 4: Process CARs
The bigWig signal file for each bulk sample is downloaded from ENCODE and used to 
generate an "output.signal" file via `bigWigAverageOverBed` on the CARs. Processed 
CARs are saved to a separate directory.

```bash
./4_Process-Bulk-CARs.sh
```

### Step 5: Create representative CARs (rCARs)
The 10th percentile of average signal over each region is computed separately for 
DNase, ATAC, and scATAC CARs. Each assay is processed individually:
1. Combined into one file
2. Filtered by [FILL IN: exact cutoff value/rationale] and excluded chromosomes (e.g., chrM)
3. Sorted, with the best peak selected per region using `pick.best.peak.py` 
   (maximum signal across all samples)

DNase, ATAC, and scATAC peaks are then merged, and a unique identifier is assigned 
to each region.

```bash
./5_Create-rCARs.sh
```

### Step 6: Retrieve signals
Average signal over CARs is retrieved for each individual assay:

```bash
./6_Retrieve-Dnase-Signal.sh
./6_Retrieve-ATAC-Signal.sh
./6_Retrieve-scATAC-Signal.sh
./6_Retrieve-H3K27ac-Signal.sh
./6_Retrieve-H3K4me3-Signal.sh
./6_Retrieve-CTCF-Signal.sh
```

### Step 7: Determine maximum z-scores
Maximum z-scores are computed for each CAR across all assay types.

```bash
./7_Determine-Max-Zscores.sh
```

### Step 8: Classify and resize cCREs
Elements are classified based on assay signal thresholds and position relative to 
transcription start sites (TSSs) [FILL IN: which TSS annotation, thresholds used 
for each class — CA, CA-CTCF, CA-H3K4me3, dELS, pELS, PLS]. Each element is then 
centered on its peak summit and resized to a fixed [FILL IN: confirm 501bp] window.

```bash
./8_Classify-cCREs.sh
```

![cCRE_example](https://github.com/user-attachments/assets/706e6311-5911-46e1-9806-06c08f92b7e6)

### Output

The final cCRE map (`cCRE_resource_master_table.tsv`) contains [FILL IN: total 
count, e.g. 80,294] candidate cis-regulatory elements, each 501bp, with the 
following columns:

| Column | Description |
|---|---|
| `chrom`, `start`, `end` | Genomic coordinates ([FILL IN: 0-based BED / confirm]) |
| `cCRE_ID` | Unique identifier (`Colon-cCRE-#`) |
| `cCRE_class` | ENCODE-style classification (CA, CA-CTCF, CA-H3K4me3, dELS, pELS, PLS) |
| `synteny_status` | Syntenic / Not Syntenic (human, mouse, dog; see Synteny) |
| `zoonomia_group` | Conservation group G1–G3, or NA (see Evolutionary sequence conservation) |

## Additional analyses

### Synteny

[FILL IN: tool(s) used — e.g., halLiftover, Cactus alignment — and a brief 
description of how syntenic regions across human, mouse, and dog were defined]

```bash
[FILL IN: script name]
```

### Evolutionary sequence conservation

[FILL IN: tool/data source — e.g., Zoonomia phyloP/phastCons scores across 240 
mammals — and how the G1/G2/G3 groups were derived]

```bash
[FILL IN: script name]
```
