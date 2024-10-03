#!/bin/bash

# Define the genome
genome=hg38

# Define the txt file with sample info
peaks=/Users/ryanhagan/NoCoSMiCC/Files/ENCODE-bulk-peaks.txt

# Set directory path(s)
dataDir=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/$genome-rCARs

cd $dataDir

# Iterate over each line in the sample info file
while IFS=$'\t' read -r exp peaks _; do

    # Check if the output file already exists
    if [ ! -f "$dataDir/CARs/$peaks.CARs.bed" ]; then

        tmp_dir=$(mktemp -d -t rhagan-XXXXXX)
        cd "$tmp_dir" || exit 1

        wget "https://www.encodeproject.org/files/$peaks/@@download/$peaks.bed.gz"

        gunzip "$peaks.bed.gz"

        awk '{if ($3-$2 >= 20 && $3-$2 <= 1000) print $0}' "$peaks.bed" > peaks

        mkdir -p "$dataDir/CARs"

        mv peaks "$dataDir/CARs/$peaks.CARs.bed"

        rm -r "$tmp_dir"
    fi

done < "$peaks"

