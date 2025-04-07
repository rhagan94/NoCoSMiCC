# NoCoSMiCC - Non-Coding Somatic Mutations in Colorectal Cancer

This repository contains the files and scripts required to enable reproduction of the analysis described in:

(link to publication)

### Abstract
Colorectal cancer (CRC) is a significant malignancy, representing the third most common cancer type and second most fatal globally. The identification of driver mutations in CRC genomes is a major research area and has a key role in facilitating the understanding of colon tumour development. Although the landscape of coding drivers has been relatively well studied, significantly less is known about the role of mutations in non-coding regions. However, recent large-scale projects have produced data that can address bottlenecks associated with the identification of non-coding driver mutations. For example, Whole genome sequencing (WGS) data availability has increased substantially, and thousands of datasets are now accessible across many cancer types, including CRC. In addition, knowledge of the non-coding genome has improved through initiatives such as the ENCODE project, which have enabled the characterisation of regulatory elements across diverse cell and tissue types.

In this study, we utilised the ENCODE database to collate data from human colon tissue and developed a pipeline to generate a catalogue of candidate cis-regulatory elements (cCREs). Additionally, we explored single cell chromatin accessibility data from the human biomolecular atlas program (HubMAP) to identify thousands of cCREs found specifically in colon epithelial cells and categorised these elements (e.g. enhancers) through integration with histone and transcription factor ChIP-seq data. Further analysis of our colon epithelial cCREs revealed differences in synteny and conservation scores (phyloP) across 241 mammals, providing insight into the potential functional significance of these elements.

### Dependencies

## 1. Data access
A key component of the initial phase of this project is the construction of colorectal regulatory element maps using data from the Encyclopedia of DNA Elements (ENCODE) project (https://www.encodeproject.org/). The ENCODE project aims to deliniate functional elements in the genome and contains integrated data across tissue types in both human and mouse, forming a catalogue of candidate cis-regulatory elements (cCREs). To create a regulatory element map for the analysis of CRC somatic mutaions, we adapt the ENCODE cCRE generation pipeline and utilise single cell chromatin accessibility data in addition to histone and TF ChIP-seq.

The following data types are used to generate colon cCREs:

- scATAC-seq
- H3K4me3 ChIP-seq
- H3K27ac ChIP-seq
- CTCF ChIP-seq

Single cell chromatin accessibility data was obtained from different sources including ENCODE, the Human Biomolcular Atlas Program (HuBMAP) and the following publications:

https://www.nature.com/articles/s41586-023-05915-x 

## 2. Running the cCRE generation pipeline

### Step 1: Generate peaks from scATAC-seq data
Fragment files are first obtained for each individual sample and used to create arrow files for an ArchR project. Quality control is performed to remove any low-quality cells from the analysis such as doublets and cells with low TSS enrichment and/or number of fragments. Accessibility around marker genes is used to select epithelial cell populations for peak calling - a peak file and bigWig file are created and exported for further analysis.

```

./0_Process-sc-fragments.R

```
### Step 2: Process single cell chromatin accessibility peaks 
Peaks from the ArchR analysis are filtered and saved as chromatin accessible regions (CARs) for each sample.

```

./1_Process-sc-CARs.sh

```

### Step 3: Download and process bulk chromatin accessibility peaks 
Bulk DNase and ATAC peaks are downloaded from ENCODE, filtered and saved as chromatin accessible regions (CARs) for each sample.

```

./2_Obtain-Bulk-Peaks.sh

```

### Step 4: Process CARs
The BigWig signal file for each bulk sample is downloaded from ENCODE and used to generate an "output.signal" file using "bigwigaverageoverbed" on the CARs. These are then saved as processed CARs in a separate directory.

```

./3_Process-Bulk-CARs

```

### Step 5: Create representative CARs (rCARs)
The 10th percentile of average signal over each region is computed seperately for DNase, ATAC and scATAC CARs and the three assays are processed individually:

- Combined into one file
- Filtered by specific cutoff and Chr names (e.g. ChrM)
- Sorted and the best peak is picked using the "pick.best.peak.py" script - this is the max signal across all samples

DNase, ATAC and scATAC peaks are then merged and a unique identifier label is added for each region.

```

./4_Create_rCARs

```

### Step 6: Retrieve signals
Retrieve the average signal over CARs for each of the individual assays

```
# Dnase
./5_Retrieve_Dnase_Signal.sh

# Bulk ATAC
./5_Retrieve_ATAC_Signal.sh

# scATAC
./5_Retrieve_scATAC_Signal.sh

# H3K27ac
./5_Retrieve_H3K27ac_Signal.sh

# H3K4me3
./5_Retrieve_H3K4me3_Signal.sh

# CTCF
./5_Retrieve_CTCF_Signal.sh

```

### Step 7: Determine max z-scores
Maximum z-scores for each CAR in all assay types are determined 

```

./6_Determine-Max-Zscores.sh

```

### Step 8: Classify cCREs
Based on assay signals and position in the genome relative to transcription start sites (TSSs)

```

./7_Classify-cCREs.sh

```
![cCRE_example](https://github.com/user-attachments/assets/706e6311-5911-46e1-9806-06c08f92b7e6)

# Additonal analyses

## Construction of a syntenic regulatory map
For the next stage of the analysis, a synteny map is constructed using the human, mouse and dog genomes. 

## Construction of TAD maps
Topologically associated domains (TADs) are evolutionarily conserved and self-interacting regions of the genome with important functions in gene regulation. Using the syntenic colorectal regulatory element map, a TAD map is created to identify elements that are conserved between TAD boundaries in the different species.


## Analysing mutations

A relatively large number of CRC WGS databases are now available and can be accessed by following the relevant data access procedures:

- Genomics England (https://www.genomicsengland.co.uk/research/research-environment)
- Hartwig Medical Foundation (https://www.hartwigmedicalfoundation.nl/en/data/data-access-request/)
- Pan Cancer Analysis of Whole Genomes (https://docs.icgc.org/pcawg/data/)

We have also included our own datasets in the analysis, accessible here:

(link to our datasets)

**Total CRC WGS datasets used:** > 3,000

### Download WGS data 
```
# Genomics England
wget https:\\

# Hartwig
wget https:\\

# PCAWG
wget https:\\

# Our data
wget https:\\
```


## TO REMOVE

ENCODE colon datasets used in this pipeline are indicated below:

**Assay** = The assay used to generate the data

**Biosample** = The biological material used as input for the experiment.

**Accession** = The specific ENCODE experiment from which the data is derived.

**Associated files** = ENCODE Ids corresponding to relevant files in each sample

|   Assay     |   Biosample  | Accession   | Associated files             
| :--------------:| :---------:|:-----------:|:-----------:
| Dnase-seq |   Transverse    | ENCSR276ITP | ENCFF274XSP.bed.gz & ENCFF299OOV.bigWig
| Dnase-seq |   Transverse    | ENCSR340MRJ | ENCFF690KZJ.bed.gz & ENCFF405NTZ.bigWig
| Dnase-seq |   Transverse    | ENCSR923JYH | ENCFF503AZD.bed.gz & ENCFF753MXL.bigWig
| Dnase-seq |   Transverse    | ENCSR763AKE | ENCFF903UTH.bed.gz & ENCFF291LZE.bigWig
| Dnase-seq |   Left colon   | ENCSR279SXQ | ENCFF042YMQ.bed.gz & ENCFF452PQB.bigWig
| Dnase-seq |   Left colon   | ENCSR867XFA | ENCFF341AEQ.bed.gz & ENCFF613BDA.bigWig
| ----------------| ---------- |-----------  | -------------------------------------
| ATAC-seq | Transverse | ENCSR386HAZ | ENCFF355GMG.bed.gz & ENCFF668GUI.bigWig
| ATAC-seq | Transverse | ENCSR404LLJ | ENCFF563UYM.bed.gz & ENCFF033RPN.bigWig
| ATAC-seq | Transverse | ENCSR761TKU | ENCFF591HMI.bed.gz & ENCFF811ERB.bigWig
| ATAC-seq | Transverse | ENCSR668VCT | ENCFF054QTC.bed.gz & ENCFF509NRA.bigWig
| ATAC-seq | Left colon | ENCSR600ZHS | ENCFF009YES.bed.gz & ENCFF057BIJ.bigWig
| ATAC-seq | Sigmoid | ENCSR355SGJ | ENCFF846PJS.bed.gz & ENCFF961XDO.bigWig
| ATAC-seq | Sigmoid | ENCSR846VLJ | ENCFF452WBJ.bed.gz & ENCFF784HME.bigWig
| ATAC-seq | Sigmoid | ENCSR086OGH | ENCFF270JNZ.bed.gz & ENCFF049WJI.bigWig
| ATAC-seq | Sigmoid | ENCSR548QCP | ENCFF018EMP.bed.gz & ENCFF796DRU.bigWig
| ----------------| ---------- |-----------  | -------------------------------------
| H3K27ac | Transverse | ENCSR792VLP | ENCFF741NZM.bigWig
| H3K27ac | Transverse | ENCSR640XRV | ENCFF532ZGB.bigWig
| H3K27ac | Transverse | ENCSR069EGE | ENCFF427MZX.bigWig
| H3K27ac | Transverse | ENCSR208QRN | ENCFF318ECM.bigWig
| H3K27ac | Sigmoid | ENCSR807XUB | ENCFF608MDC.bigWig
| H3K27ac | Sigmoid | ENCSR641SDI | ENCFF468UEP.bigWig
| H3K27ac | Sigmoid | ENCSR937EVN | ENCFF322NLT.bigWig
| H3K27ac | Sigmoid | ENCSR268ZCF | ENCFF111DLN.bigWig
| H3K27ac | Sigmoid | ENCSR561YSH | ENCFF124AKT.bigWig
| ----------------| ---------- |-----------  | -------------------------------------
| H3K4me3 | Transverse | ENCSR315EZG | ENCFF487CTD.bigWig
| H3K4me3 | Transverse | ENCSR557OWY | ENCFF568IBR.bigWig
| H3K4me3 | Transverse | ENCSR813ZEY | ENCFF252OBP.bigWig
| H3K4me3 | Transverse | ENCSR933BVL | ENCFF339CRV.bigWig
| H3K4me3 | Sigmoid | ENCSR960AAL | ENCFF122GSI.bigWig
| H3K4me3 | Sigmoid | ENCSR793IKH | ENCFF800VAP.bigWig
| H3K4me3 | Sigmoid | ENCSR172LVU | ENCFF886LUE.bigWig
| H3K4me3 | Sigmoid | ENCSR900UIP | ENCFF237VMY.bigWig
| H3K4me3 | Sigmoid | ENCSR421HUB | ENCFF423YBA.bigWig
| ----------------| ---------- |-----------  | -------------------------------------
| CTCF | Transverse | ENCSR449SEF | ENCFF626YRZ.bigWig
| CTCF | Transverse | ENCSR907BES | ENCFF435CDF.bigWig
| CTCF | Transverse | ENCSR236YGF | ENCFF686TFX.bigWig
| CTCF | Transverse | ENCSR608WPS | ENCFF646EZE.bigWig
| CTCF | Transverse | ENCSR558HTE | ENCFF349LWW.bigWig
| CTCF | Transverse | ENCSR769WKR | ENCFF493XMW.bigWig
| CTCF | Transverse | ENCSR833FWC | ENCFF517GQE.bigWig
| CTCF | Transverse | ENCSR102CSD | ENCFF170KVV.bigWig
| CTCF | Sigmoid | ENCSR857RJQ | ENCFF634JUC.bigWig
| CTCF | Sigmoid | ENCSR925GDS | ENCFF154FOF.bigWig
| CTCF | Sigmoid | ENCSR721AHD | ENCFF460ECX.bigWig
| CTCF | Sigmoid | ENCSR222SQE | ENCFF985EXZ.bigWig
| ----------------| ---------- |-----------  | -------------------------------------
| scATAC-seq | Transverse | ENCSR997YNO | ENCFF490PBS.fragments.tsv.gz
| scATAC-seq | Transverse | ENCSR506YMX | ENCFF389FMC.fragments.tsv.gz
| scATAC-seq | Transverse | ENCSR434SXE | ENCFF888LMJ.fragments.tsv.gz 
| scATAC-seq | Transverse | ENCSR349XKD | ENCFF551UBJ.fragments.tsv.gz 
| scATAC-seq | Sigmoid | ENCSR388NCA | ENCFF994HOW.fragments.tsv.gz 
| scATAC-seq | Sigmoid | ENCSR367GKP | ENCFF010JFQ.fragments.tsv.gz 
| scATAC-seq | Left colon | ENCSR916RYB | ENCFF810NQZ.fragments.tsv.gz 
| scATAC-seq | Left colon | ENCSR830FPR | ENCFF040VYY.fragments.tsv.gz 

Additional information on the data organisation structure within the ENCODE project can be found here: https://www.encodeproject.org/help/data-organization/ 

### Note
Files containing all of the relevant IDs for each dataset are located in the ```Files``` directory of this repository. These files are used to locate the relevant data for downloading and processing during the pipeline run.



