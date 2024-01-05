# NoCoSMiCC - Non-Coding Somatic Mutations in Colorectal Cancer

This repository contains the files and scripts required to enable reproduction of the analysis described in:

(link to publication)

## Data access

The analysis of non-coding somatic mutations relevant to colorectal cancer (CRC) first requires collation of Whole Genome Sequencing (WGS) data. A large number of CRC WGS databases are now available and can be accessed by following the relevant data access procedures:

- Genomics England (https://www.genomicsengland.co.uk/research/research-environment)
- Hartwig Medical Foundation (https://www.hartwigmedicalfoundation.nl/en/data/data-access-request/)
- Pan Cancer Analysis of Whole Genomes (https://docs.icgc.org/pcawg/data/)

We have also included our own datasets in the analysis, accessible here:

(link to our datasets)

**Total CRC WGS datasets used:** > 3,000

## Download WGS data 
```
# Genomics England
wget https:\\

# Hartwig
wget https:\\

# PCAWG
wget https:\\
```

## Mapping of conserved regulatory elements 
A key component of the initial phase of this project is the construction of colorectal regulatory element maps using data from the Encyclopedia of DNA Elements (ENCODE) project (https://www.encodeproject.org/). The ENCODE project aims to deliniate functional elements in the genome and contains integrated data across tissue types in human and mouse, forming a catalogue of cis-regulatory elements (CREs). To create specific regulatory element maps for the analysis of CRC somatic mutaions in the two species, we utilise data originating only from normal colonic tissue.

Relevant ENCODE information used in the study is outlined in the table below:

|   Biosample     |   Species  | Accession                
| ----------------| ---------- |-----------
| Large intestine |   Mouse    | ENCSR557MSF
| Large intestine |   Human    | ENCSR073QUR
| Large intestine |   Human    | ENCSR751PAL
| Large intestine |   Human    | ENCSR229WYJ
| Large intestine |   Human    | ENCSR512QPR
| Large intestine |   Human    | ENCSR206POU
| Large intestine |   Human    | ENCSR834NCV
| Large intestine |   Human    | ENCSR389UFY
| Large intestine |   Human    | ENCSR481AZU
| Large intestine |   Human    | ENCSR207LTA
| Large intestine |   Human    | ENCSR211XAO
| Large intestine |   Human    | ENCSR735GNP
| Large intestine |   Human    | ENCSR373KEK
| Large intestine |   Human    | ENCSR033DHK
| Large intestine |   Human    | ENCSR368VZW
| ----------------| ---------- |-----------
|   Intestine     |   Mouse    | ENCSR089NNM
|   Intestine     |   Mouse    | ENCSR496ZLW
|   Intestine     |   Mouse    | ENCSR884KZR
|   Intestine     |   Mouse    | ENCSR371LBI
| ----------------| ---------- |-----------
| Sigmoid colon   |   Human    | ENCSR102HGJ
| Sigmoid colon   |   Human    | ENCSR855THX
| Sigmoid colon   |   Human    | ENCSR191LCU
| Sigmoid colon   |   Human    | ENCSR328EZS
| Sigmoid colon   |   Human    | ENCSR752MGN
| Sigmoid colon   |   Human    | ENCSR300XSY
| ----------------| ---------- |-----------
| Transverse colon|   Human    | ENCSR364HRH
| Transverse colon|   Human    | ENCSR884ABV
| Transverse colon|   Human    | ENCSR811EFU
| Transverse colon|   Human    | ENCSR019GHX
| ----------------| ---------- |-----------
| Colonic mucosa  |   Human    | ENCSR269SMA
| Colonic mucosa  |   Human    | ENCSR830QVV
| Colonic mucosa  |   Human    | ENCSR834AVK
| ----------------| ---------- |-----------
| Descending colon mucosa|   Human    | ENCSR116DCM
| Descending colon mucosa|   Human    | ENCSR417SFM
| Descending colon mucosa|   Human    | ENCSR267EIV
| ----------------| ---------- |-----------
| Left colon|   Human    | ENCSR696WQN
| Left colon|   Human    | ENCSR103FRM
| ----------------| ---------- |-----------
| Muscle layer of colon|   Human    | ENCSR196DHC
| Muscle layer of colon|   Human    | ENCSR412JVQ

**Total human CRE files** = 34

**Total mouse CRE files** = 5

### Dowload CRE data from ENCODE
A file containing all of the URLs necessary to download human and mouse CRE data is located in the ```Files``` directory of this repository. To download CRE data run the following code:

```
xargs -L 1 curl -O -J -L < ./Files/encode_cre_files.txt
```
