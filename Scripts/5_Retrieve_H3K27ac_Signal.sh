#!/bin/bash

# Download the necessary bigWig files for each sample and store locally
cd /Users/ryanhagan/NoCoSMiCC/raw_data/Histone\ ChIPseq/Human/bigWigs/H3K27ac
# Transverse colon
wget https://www.encodeproject.org/files/ENCFF741NZM/@@download/ENCFF741NZM.bigWig &
wget https://www.encodeproject.org/files/ENCFF532ZGB/@@download/ENCFF532ZGB.bigWig &
wget https://www.encodeproject.org/files/ENCFF427MZX/@@download/ENCFF427MZX.bigWig &
wget https://www.encodeproject.org/files/ENCFF318ECM/@@download/ENCFF318ECM.bigWig &

# Sigmoid colon
wget https://www.encodeproject.org/files/ENCFF608MDC/@@download/ENCFF608MDC.bigWig &
wget https://www.encodeproject.org/files/ENCFF468UEP/@@download/ENCFF468UEP.bigWig &
wget https://www.encodeproject.org/files/ENCFF322NLT/@@download/ENCFF322NLT.bigWig &
wget https://www.encodeproject.org/files/ENCFF111DLN/@@download/ENCFF111DLN.bigWig &
wget https://www.encodeproject.org/files/ENCFF124AKT/@@download/ENCFF124AKT.bigWig

cd /Users/ryanhagan/NoCoSMiCC

# Set variables
genome=hg38
mode=H3K27ac

if [ $mode == "DNase" ] || [ $mode == "CTCF" ] || [ $mode == "POL2" ] || [ $mode == "ATAC" ]
then
    width=0
else
    width=500
fi
echo $width

peaks=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/$genome-rCARs/hg38-rCARs-v3-labelled.bed

files=/Users/ryanhagan/NoCoSMiCC/Files/$mode-List_colon.txt
output=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/H3K27ac

dataDir=/Users/ryanhagan/NoCoSMiCC/raw_data/Histone_ChIPseq/Human/bigWigs/H3K27ac
scriptDir=/Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication

while IFS=$'\t' read -r dset dsig _; do

    # Create a temporary directory and change to it
    tmp_dir=$(mktemp -d -t rhagan-XXXXXX)
    cd "$tmp_dir" || exit 1

awk -F "\t" '{printf "%s\t%.0f\t%.0f\t%s\n", $1,$2-'$width',$3+'$width',$4}' \
$peaks | awk '{if ($2 < 0) print $1 "\t" 0 "\t" $3 "\t" $4 ; else print $0}' \
| sort -u > little

/Users/ryanhagan/miniconda3/envs/ryan_env/bin/bigwigaverageoverbed -bedOut=out2.bed $dataDir/$dsig.bigWig little out2

python $scriptDir/log.zscore.normalization.py out2 > l

sort -k2,2rg l | awk 'BEGIN {rank=0; before=0; running=1}{if ($2 != before) \
    rank = running; print $1 "\t" $2 "\t" $3 "\t" rank; before=$2; \
    running += 1}' | sort -k1,1 > $output/$dset"-"$dsig".txt"

rm -r "$tmp_dir"

done < "$files"





