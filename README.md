# NoCoSMiCC - Non-Coding Somatic Mutations in Colorectal Cancer

This repository contains the files and scripts required to enable reproduction of the analysis described in:

(link to publication)

### Abstract

### Dependencies

## 1. Data access

The analysis of non-coding somatic mutations implictaed in colorectal cancer (CRC) development first requires collation of Whole Genome Sequencing (WGS) data. A large number of CRC WGS databases are now available and can be accessed by following the relevant data access procedures:

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

## 2. Mapping of conserved regulatory elements 
A key component of the initial phase of this project is the construction of colorectal regulatory element maps using data from the Encyclopedia of DNA Elements (ENCODE) project (https://www.encodeproject.org/). The ENCODE project aims to deliniate functional elements in the genome and contains integrated data across tissue types in both human and mouse, forming a catalogue of cis-regulatory elements (CREs). To create specific regulatory element maps for the analysis of CRC somatic mutaions in the two species, we utilise data originating only from *normal colonic tissue*.

Relevant ENCODE information used in the study is outlined in the table below:

|   Biosample     |   Species  | Accession   | Associated files             
| :--------------:| :---------:|:-----------:|:-----------:
| Large intestine |   Mouse    | ENCSR557MSF | ENCFF152AFP.bed.gz & ENCFF954ZED.bigBed
| Large intestine |   Human    | ENCSR073QUR | ENCFF024PQS.bed.gz & ENCFF778PTM.bigBed
| Large intestine |   Human    | ENCSR751PAL | ENCFF034BHY.bed.gz & ENCFF524SID.bigBed
| Large intestine |   Human    | ENCSR229WYJ | ENCFF572DWG.bed.gz & ENCFF281JFE.bigBed
| Large intestine |   Human    | ENCSR512QPR | ENCFF205FKB.bed.gz & ENCFF970REE.bigbed
| Large intestine |   Human    | ENCSR206POU | ENCFF661PLY.bed.gz & ENCFF048AGW.bigBed
| Large intestine |   Human    | ENCSR834NCV | ENCFF666XFC.bed.gz & ENCFF898SJN.bigBed
| Large intestine |   Human    | ENCSR389UFY | ENCFF025KEO.bed.gz & ENCFF850JRQ.bigBed
| Large intestine |   Human    | ENCSR481AZU | ENCFF900LHN.bed.gz & ENCFF169NIT.bigBed
| Large intestine |   Human    | ENCSR207LTA | ENCFF803TAS.bed.gz & ENCFF362GPG.bigBed
| Large intestine |   Human    | ENCSR211XAO | ENCFF284MXI.bed.gz & ENCFF063IWT.bigBed
| Large intestine |   Human    | ENCSR735GNP | ENCFF715SJC.bed.gz & ENCFF118IPF.bigBed
| Large intestine |   Human    | ENCSR373KEK | ENCFF257QEB.bed.gz & ENCFF312SCY.bigBed
| Large intestine |   Human    | ENCSR033DHK | ENCFF633GHY.bed.gz & ENCFF661FXP.bigBed
| Large intestine |   Human    | ENCSR368VZW | ENCFF841KWY.bed.gz & ENCFF669TTI.bigBed
| ----------------| ---------- |-----------  | -------------------------------------
|   Intestine     |   Mouse    | ENCSR089NNM | ENCFF974UWG.bed.gz & ENCFF824OZP.bigBed
|   Intestine     |   Mouse    | ENCSR496ZLW | ENCFF751RMP.bed.gz & ENCFF604EPU.bigBed
|   Intestine     |   Mouse    | ENCSR884KZR | ENCFF556NXC.bed.gz & ENCFF278JOB.bigBed
|   Intestine     |   Mouse    | ENCSR371LBI | ENCFF931ZNP.bed.gz & ENCFF479RIY.bigBed
| ----------------| ---------- |-----------  | -------------------------------------
| Sigmoid colon   |   Human    | ENCSR102HGJ | ENCFF516ZJM.bed.gz & ENCFF244IEB.bigBed
| Sigmoid colon   |   Human    | ENCSR855THX | ENCFF507HVH.bed.gz & ENCFF386TZM.bigBed
| Sigmoid colon   |   Human    | ENCSR191LCU | ENCFF277SLP.bed.gz & ENCFF844JZF.bigBed
| Sigmoid colon   |   Human    | ENCSR328EZS | ENCFF466PWQ.bed.gz & ENCFF057KVE.bigBed
| Sigmoid colon   |   Human    | ENCSR752MGN | ENCFF322DZO.bed.gz & ENCFF838QAW.bigBed
| Sigmoid colon   |   Human    | ENCSR300XSY | ENCFF137THP.bed.gz & ENCFF571FOS.bigbed
| ----------------| ---------- |-----------  | -------------------------------------
| Transverse colon|   Human    | ENCSR364HRH | ENCFF164RLV.bed.gz & ENCFF224CJA.bigBed
| Transverse colon|   Human    | ENCSR884ABV | ENCFF424YXV.bed.gz & ENCFF186HRJ.bigBed
| Transverse colon|   Human    | ENCSR811EFU | ENCFF937AZG.bed.gz & ENCFF488QJI.bigbed
| Transverse colon|   Human    | ENCSR019GHX | ENCFF606DOD.bed.gz & ENCFF790PDO.bigBed
| ----------------| ---------- |-----------  | -------------------------------------
| Colonic mucosa  |   Human    | ENCSR269SMA | ENCFF898UAX.bed.gz & ENCFF389JSK.bigBed
| Colonic mucosa  |   Human    | ENCSR830QVV | ENCFF316MRH.bed.gz & ENCFF163GCV.bigBed
| Colonic mucosa  |   Human    | ENCSR834AVK | ENCFF475ITB.bed.gz & ENCFF946AMW.bigBed
| ----------------| ---------- |-----------  | -------------------------------------
| Descending colon mucosa|   Human    | ENCSR116DCM | ENCFF632LJU.bed.gz & ENCFF582GTV.bigBed
| Descending colon mucosa|   Human    | ENCSR417SFM | ENCFF350NWA.bed.gz & ENCFF045ZQC.bigBed
| Descending colon mucosa|   Human    | ENCSR267EIV | ENCFF095SMO.bed.gz & ENCFF028NAN.bigbed
| ----------------| ---------- |----------- | -------------------------------------
| Left colon|   Human    | ENCSR696WQN | ENCFF318JTT.bed.gz & ENCFF187RLW.bigBed
| Left colon|   Human    | ENCSR103FRM | ENCFF954AFR.bed.gz & ENCFF046TXU.bigBed
| ----------------| ---------- |----------- | -------------------------------------
| Muscle layer of colon|   Human    | ENCSR196DHC | ENCFF750WMI.bed.gz & ENCFF375JWG.bigBed
| Muscle layer of colon|   Human    | ENCSR412JVQ | ENCFF360RXI.bed.gz & ENCFF935DDY.bigBed

**Total human CRE files** = 34

**Total mouse CRE files** = 5

### Dowload CRE data from ENCODE
A file containing all of the URLs necessary to download human and mouse CRE data shown in the table above is located in the ```Files``` directory of this repository. To download CRE data run the following code:

```
# Create a new directory to store the CRE data files
mkdir ./raw_data/encode_cres

# Switch to the new directory and download
cd ./raw_data/encode_cres
xargs -L 1 curl -O -J -L < ./Files/encode_cre_files.txt
```

## 3) Construction of a syntenic regulatory map
For the next stage of the analysis, a synteny map is constructed using the human, mouse and dog genomes. 

## 4) Construction of TAD maps
Topologically associated domains (TADs) are evolutionarily conserved and self-interacting regions of the genome with important functions in gene regulation. Using the syntenic colorectal regulatory element map, a TAD map is created to identify elements that are conserved between TAD boundaries in the different species.






