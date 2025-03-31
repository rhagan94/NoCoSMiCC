#!/bin/bash

# Adapted from:

# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for aggregating phyloP scores on human cCREs

# Get the PhyloP BigWig file
system('wget "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/cactus241way/cactus241way.phyloP.bw" -P /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files')

# Reformat labels (if necessary)
# awk '{gsub(/_/,"",$4)}1' /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3-reformat.bed

cCRE_file=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs-2025.bed
phyloP_signal=/Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/cactus241way.phyloP.bw
output_matrix=/Users/ryanhagan/NoCoSMiCC/Zoonomia/Outputs/phyloP241_matrix.mtx

scriptDir=/Users/ryanhagan/NoCoSMiCC/Zoonomia/Scripts

n=500 # total number of bins

## cut cCRE-center±250bp into 1bp bins
awk -v n="$n" '{FS=OFS="\t"}{
    center=int(($2+$3)/2);
    start=center-250;
    for(i=1;i<=n;i++){
        print $1,int(start+(i-1)),int(start+i),$4"_"i;1;"."
    }
}' <(cut -f 1-4 ${cCRE_file}) > tmp.ccREs_split.bed 

## calcualte signal for each bin
bigWigAverageOverBed ${phyloP_signal} tmp.ccREs_split.bed tmp.ccREs_split.tab
cut -f 1,5 tmp.ccREs_split.tab | awk '{FS=OFS="\t"}{
    split($1,a,"_");print $1,a[1],a[2],$2
}' | sort -k2,2 -k3,3n > tmp.ccREs_signal.txt

## make matrix
awk -v n="$n" '{FS=OFS="\t"}{if(NR%n==1){split($1,a,"_");printf a[1]"\t"$4"\t"}else if(NR%n!=0){printf $4"\t"}else{printf $4"\n"}}'  tmp.ccREs_signal.txt > ${output_matrix}

## remove the tmp files
rm tmp.ccREs_split.bed tmp.ccREs_signal.txt

## END OF SCRIPT ##
