#!/bin/bash

# Change directory and download the necessary bigWig files for each sample and store locally
cd /Users/ryanhagan/NoCoSMiCC/raw_data/atac_data/bigWigs

wget https://www.encodeproject.org/files/ENCFF668GUI/@@download/ENCFF668GUI.bigWig &
wget https://www.encodeproject.org/files/ENCFF033RPN/@@download/ENCFF033RPN.bigWig &
wget https://www.encodeproject.org/files/ENCFF811ERB/@@download/ENCFF811ERB.bigWig &
wget https://www.encodeproject.org/files/ENCFF509NRA/@@download/ENCFF509NRA.bigWig &
wget https://www.encodeproject.org/files/ENCFF057BIJ/@@download/ENCFF057BIJ.bigWig &
wget https://www.encodeproject.org/files/ENCFF961XDO/@@download/ENCFF961XDO.bigWig &
wget https://www.encodeproject.org/files/ENCFF784HME/@@download/ENCFF784HME.bigWig &
wget https://www.encodeproject.org/files/ENCFF049WJI/@@download/ENCFF049WJI.bigWig &
wget https://www.encodeproject.org/files/ENCFF796DRU/@@download/ENCFF796DRU.bigWig

cd /Users/ryanhagan/NoCoSMiCC

# Set variables
genome=hg38
mode=ATAC

if [ $mode == "DNase" ] || [ $mode == "CTCF" ] || [ $mode == "POL2" ] || [ $mode == "ATAC" ]
then
    width=0
else
    width=500
fi
echo $width

peaks=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/$genome-rCARs/hg38-rCARs-v3-labelled.bed

files=/Users/ryanhagan/NoCoSMiCC/Files/$mode-List_colon.txt
output=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/signal_output/ATAC

dataDir=/Users/ryanhagan/NoCoSMiCC/raw_data/atac_data/bigWigs
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





