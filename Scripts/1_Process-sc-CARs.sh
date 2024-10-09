#!/bin/bash

# Define the genome and directories
genome=hg38

# Define the data directory, script directory, and hotspots file
dataDir=/Users/ryanhagan/NoCoSMiCC/ArchR_analysis
scriptDir=/Users/ryanhagan/NoCoSMiCC/Scripts/Custom_ENCODE_Pipeline_V1
samples=/Users/ryanhagan/NoCoSMiCC/Files/singlecell-to-run.txt
signalDir=/Users/ryanhagan/NoCoSMiCC/ArchR_analysis/GroupBigWigs
outdir=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs

# Loop through each line in the hotspots file
while IFS=$'\t' read -r dset; do

    # Create a temporary directory and change to it
    tmp_dir=$(mktemp -d -t rhagan-XXXXXX)
    cd "$tmp_dir" || exit 1

    signal=$signalDir/$dset-ArchR.bw

    # Prepare the CARs file
    CAR=$dataDir/MACS2_peaks/$dset.epithelial.bed 
    cp "$CAR" bed
    
    awk '{print $1 "\t" $2 "\t" $3 "\t" "'$dset'-"NR "\t" $4}' bed | sort -k4,4 > new
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' new > new.bed
    /Users/ryanhagan/miniconda3/envs/ryan_env/bin/bigwigaverageoverbed "$signal" new.bed out.tab
    sort -k1,1 out.tab > tmp.1
    paste new tmp.1 | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $10}' > output.sc.$dset

    # Move the output to the data directory
    mkdir -p "$outdir/Processed-CARs/"
    mv output.sc.$dset "$outdir/Processed-CARs/"

    # Clean up the temporary directory
    rm -r "$tmp_dir"

done < "$samples"
