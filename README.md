# NoCoSMiCC - Non-Coding Somatic Mutations in Colorectal Cancer

This repository contains the files and scripts required to enable reproduction of the analysis described in:

(link to publication)

### Abstract

### Dependencies

## 1. Data access

The analysis of non-coding somatic mutations involved in colorectal cancer (CRC) development first requires the collation of Whole Genome Sequencing (WGS) data from CRC patients. A relatively large number of CRC WGS databases are now available and can be accessed by following the relevant data access procedures:

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

**Biosample** = The biological material used as input for the experiment.

**Species** = The organism used in the experiment - human or mouse.

**Accession** = The specific ENCODE assay (experiment) from which the data is derived.

**Associated files** = Files available for download for each accession. 


Sample number |   Biosample     |   Species  | Accession   | Associated files             
| :--------------:| :--------------:| :---------:|:-----------:|:-----------:
| 1| Large intestine |   Mouse    | ENCSR557MSF | ENCFF152AFP.bed.gz & ENCFF954ZED.bigBed
| 2| Large intestine |   Human    | ENCSR073QUR | ENCFF024PQS.bed.gz & ENCFF778PTM.bigBed
| 3| Large intestine |   Human    | ENCSR751PAL | ENCFF034BHY.bed.gz & ENCFF524SID.bigBed
| 4| Large intestine |   Human    | ENCSR229WYJ | ENCFF572DWG.bed.gz & ENCFF281JFE.bigBed
| 5| Large intestine |   Human    | ENCSR512QPR | ENCFF205FKB.bed.gz & ENCFF970REE.bigbed
| 6| Large intestine |   Human    | ENCSR206POU | ENCFF661PLY.bed.gz & ENCFF048AGW.bigBed
| 7| Large intestine |   Human    | ENCSR834NCV | ENCFF666XFC.bed.gz & ENCFF898SJN.bigBed
| 8| Large intestine |   Human    | ENCSR389UFY | ENCFF025KEO.bed.gz & ENCFF850JRQ.bigBed
| 9| Large intestine |   Human    | ENCSR481AZU | ENCFF900LHN.bed.gz & ENCFF169NIT.bigBed
| 10| Large intestine |   Human    | ENCSR207LTA | ENCFF803TAS.bed.gz & ENCFF362GPG.bigBed
| 11| Large intestine |   Human    | ENCSR211XAO | ENCFF284MXI.bed.gz & ENCFF063IWT.bigBed
| 12| Large intestine |   Human    | ENCSR735GNP | ENCFF715SJC.bed.gz & ENCFF118IPF.bigBed
| 13| Large intestine |   Human    | ENCSR373KEK | ENCFF257QEB.bed.gz & ENCFF312SCY.bigBed
| 14| Large intestine |   Human    | ENCSR033DHK | ENCFF633GHY.bed.gz & ENCFF661FXP.bigBed
| 15| Large intestine |   Human    | ENCSR368VZW | ENCFF841KWY.bed.gz & ENCFF669TTI.bigBed
| ----------------| ---------- |-----------  | -------------------------------------
| 16|  Intestine     |   Mouse    | ENCSR089NNM | ENCFF974UWG.bed.gz & ENCFF824OZP.bigBed
| 17|  Intestine     |   Mouse    | ENCSR496ZLW | ENCFF751RMP.bed.gz & ENCFF604EPU.bigBed
| 18|  Intestine     |   Mouse    | ENCSR884KZR | ENCFF556NXC.bed.gz & ENCFF278JOB.bigBed
| 19|  Intestine     |   Mouse    | ENCSR371LBI | ENCFF931ZNP.bed.gz & ENCFF479RIY.bigBed
| ----------------| ---------- |-----------  | -------------------------------------
| 20| Sigmoid colon   |   Human    | ENCSR102HGJ | ENCFF516ZJM.bed.gz & ENCFF244IEB.bigBed
| 21| Sigmoid colon   |   Human    | ENCSR855THX | ENCFF507HVH.bed.gz & ENCFF386TZM.bigBed
| 22| Sigmoid colon   |   Human    | ENCSR191LCU | ENCFF277SLP.bed.gz & ENCFF844JZF.bigBed
| 23| Sigmoid colon   |   Human    | ENCSR328EZS | ENCFF466PWQ.bed.gz & ENCFF057KVE.bigBed
| 24| Sigmoid colon   |   Human    | ENCSR752MGN | ENCFF322DZO.bed.gz & ENCFF838QAW.bigBed
| 25| Sigmoid colon   |   Human    | ENCSR300XSY | ENCFF137THP.bed.gz & ENCFF571FOS.bigbed
| ----------------| ---------- |-----------  | -------------------------------------
| 26| Transverse colon|   Human    | ENCSR364HRH | ENCFF164RLV.bed.gz & ENCFF224CJA.bigBed
| 27| Transverse colon|   Human    | ENCSR884ABV | ENCFF424YXV.bed.gz & ENCFF186HRJ.bigBed
| 28| Transverse colon|   Human    | ENCSR811EFU | ENCFF937AZG.bed.gz & ENCFF488QJI.bigbed
| 29| Transverse colon|   Human    | ENCSR019GHX | ENCFF606DOD.bed.gz & ENCFF790PDO.bigBed
| ----------------| ---------- |-----------  | -------------------------------------
| 30| Colonic mucosa  |   Human    | ENCSR269SMA | ENCFF898UAX.bed.gz & ENCFF389JSK.bigBed
| 31| Colonic mucosa  |   Human    | ENCSR830QVV | ENCFF316MRH.bed.gz & ENCFF163GCV.bigBed
| 32| Colonic mucosa  |   Human    | ENCSR834AVK | ENCFF475ITB.bed.gz & ENCFF946AMW.bigBed
| ----------------| ---------- |-----------  | -------------------------------------
| 33| Descending colon mucosa|   Human    | ENCSR116DCM | ENCFF632LJU.bed.gz & ENCFF582GTV.bigBed
| 34| Descending colon mucosa|   Human    | ENCSR417SFM | ENCFF350NWA.bed.gz & ENCFF045ZQC.bigBed
| 35| Descending colon mucosa|   Human    | ENCSR267EIV | ENCFF095SMO.bed.gz & ENCFF028NAN.bigbed
| ----------------| ---------- |----------- | -------------------------------------
| 36| Left colon|   Human    | ENCSR696WQN | ENCFF318JTT.bed.gz & ENCFF187RLW.bigBed
| 37| Left colon|   Human    | ENCSR103FRM | ENCFF954AFR.bed.gz & ENCFF046TXU.bigBed
| ----------------| ---------- |----------- | -------------------------------------
| 38| Muscle layer of colon|   Human    | ENCSR196DHC | ENCFF750WMI.bed.gz & ENCFF375JWG.bigBed
| 39| Muscle layer of colon|   Human    | ENCSR412JVQ | ENCFF360RXI.bed.gz & ENCFF935DDY.bigBed
| ----------------| ---------- |----------- | -------------------------------------
| 40| Rectal mucosa|   Human    | ENCSR242DUB | ENCFF842GND.bed.gz & ENCFF141OKN.bigBed
| 41| Rectal mucosa|   Human    | ENCSR556FAX | ENCFF478BCY.bed.gz & ENCFF738NDA.bigBed
| ----------------| ---------- |----------- | -------------------------------------
| 42| Rectal smooth muscle|   Human    | ENCSR843VDL | ENCFF653ZCA.bed.gz & ENCFF620THB.bigBed

**Total human cCRE files** = 37

**Total mouse cCRE files** = 5

Additional information on the data organisation structure within the ENCODE project can be found here: https://www.encodeproject.org/help/data-organization/ 

### Dowload cCRE data from ENCODE
A file containing all of the URLs necessary to download human and mouse CRE data shown in the table above is located in the ```Files``` directory of this repository. To download CRE data run the following code:

```
# Create a new directory to store the CRE data files
mkdir -p ./raw_data/encode_cres

# Switch to the new directory and download the data from ENCODE
cd ./raw_data/encode_cres
xargs -L 1 curl -O -J -L < ./Files/encode_cre_files.txt
```

## 3) Construction of a syntenic regulatory map
For the next stage of the analysis, a synteny map is constructed using the human, mouse and dog genomes. 

### Download the 'chain files' from UCSC for human, mouse and dog comparisons
```
# Make a directory to store the chain files
mkdir reference_genome/ucsc_chain

# Mouse to human
wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/liftOver/mm39ToHg38.over.chain.gz" -P reference_genome/ucsc_chain

# Dog to human
wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/liftOver/canFam6ToHg38.over.chain.gz" -P reference_genome/ucsc_chain

```

## 4) Construction of TAD maps
Topologically associated domains (TADs) are evolutionarily conserved and self-interacting regions of the genome with important functions in gene regulation. Using the syntenic colorectal regulatory element map, a TAD map is created to identify elements that are conserved between TAD boundaries in the different species.






