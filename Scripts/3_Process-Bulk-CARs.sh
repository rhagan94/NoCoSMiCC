#!/bin/bash

# Define the genome and directories
genome=hg38

# Define the data directory, script directory, and hotspots file
dataDir=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/$genome-rCARs
scriptDir=/Users/ryanhagan/NoCoSMiCC/Scripts/Custom_ENCODE_Pipeline_V1
peaks=/Users/ryanhagan/NoCoSMiCC/Files/peaks-to-run.txt

# Loop through each line in the hotspots file
while IFS=$'\t' read -r dset dpeak dsig _; do

    # Create a temporary directory and change to it
    tmp_dir=$(mktemp -d -t rhagan-XXXXXX)
    cd "$tmp_dir" || exit 1

    # Check if the signal file exists locally, otherwise download it
    if [ -f /data/projects/encode/data/$dset/$dsig.bigWig ]; then
        signal=/data/projects/encode/data/$dset/$dsig.bigWig
    else
        wget https://www.encodeproject.org/files/$dsig/@@download/$dsig.bigWig
        signal=$dsig.bigWig
    fi

    # Prepare the DHS file
    dhs=$dataDir/CARs/$dpeak.CARs.bed 
    cp "$dhs" bed
    awk '{print $1 "\t" $2 "\t" $3 "\t" "'$dpeak'-"NR "\t" $4}' bed | sort -k4,4 > new
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' new > new.bed
    /Users/ryanhagan/miniconda3/envs/ryan_env/bin/bigwigaverageoverbed "$signal" new.bed out.tab
    sort -k1,1 out.tab > tmp.1
    paste new tmp.1 | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $10}' > output.$dsig

    # Move the output to the data directory
    mkdir -p "$dataDir/Processed-CARs/"
    mv output.$dsig "$dataDir/Processed-CARs/"

    # Clean up the temporary directory
    rm -r "$tmp_dir"

done < "$peaks"
