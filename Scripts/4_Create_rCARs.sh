#!/bin/bash

genome=hg38
dir=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/$genome-rCARs
scriptDir=/Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication

cd $dir/Processed-CARs

# Calculate the 10th percentile of average signal in DNase regions
paste output.ENCFF291LZE output.ENCFF299OOV output.ENCFF405NTZ output.ENCFF452PQB output.ENCFF613BDA output.ENCFF753MXL | awk '{ print $6}' > merged_DNase_signal.txt
sort -n merged_DNase_signal.txt  | awk '{all[NR] = $0} END{print all[int(NR*0.1 - 0.5)]}' # 0.110121

# Calculate the 10th percentile of average signal in ATAC regions
paste output.ENCFF668GUI output.ENCFF033RPN output.ENCFF811ERB output.ENCFF509NRA output.ENCFF057BIJ output.ENCFF961XDO \
output.ENCFF784HME output.ENCFF049WJI output.ENCFF796DRU | awk '{ print $6}' > merged_ATAC_signal.txt
sort -n merged_ATAC_signal.txt  | awk '{all[NR] = $0} END{print all[int(NR*0.1 - 0.5)]}' # 2.2947

# Calculate the 10th percentile of average signal in scATAC regions
paste output.sc.Transverse_ENCSR349XKD output.sc.Transverse_ENCSR434SXE output.sc.Transverse_ENCSR506YMX \
output.sc.Transverse_ENCSR997YNO output.sc.Left_ENCSR830FPR | awk '{ print $6}' > merged_scATAC_signal.txt
sort -n merged_scATAC_signal.txt  | awk '{all[NR] = $0} END{print all[int(NR*0.1 - 0.5)]}' # 0.00688278

## Process the DNase data
echo -e "Combining DNase peaks..."
cat output.ENCFF291LZE output.ENCFF299OOV output.ENCFF405NTZ output.ENCFF452PQB output.ENCFF613BDA output.ENCFF753MXL > tmp
mv tmp $dir/$genome-DNase-All.bed

if [ $genome == "mm10" ]
then
#cutoff=0.1454
dnase_cutoff=0.158899 #2021 cutoff
elif [ $genome == "hg38" ]
then
dnase_cutoff=0.110121 #custom cutoff

fi

echo -e "Filtering DNase peaks..."
cd $dir
awk '{if ($1 !~ /_/ && $3-$2 >= 20 && $6 >= '$dnase_cutoff') print $0}' \
    $genome-DNase-All.bed | grep -v "chrEBV" | grep -v "chrM" > $genome-DNase-Filtered.bed

mkdir scratch
cp $genome-DNase-Filtered.bed scratch/tmp.bed
cd scratch

echo -e "Sorting DNase peaks..."
sort -k1,1 -k2,2n tmp.bed > sorted
rm -f rPeaks
num=$(wc -l sorted | awk '{print $1}')

echo -e "Merging DNase peaks..."
while [ $num -gt 0 ]
do
    echo -e "\t" $num
    bedtools merge -i sorted -c 4,6 -o collapse,collapse > merge
    python $scriptDir/pick.best.peak.py merge > peak-list
    awk 'FNR==NR {x[$1];next} ($4 in x)' peak-list sorted >> rPeaks
    bedtools intersect -v -a sorted -b rPeaks > remaining
    mv remaining sorted
    num=$(wc -l sorted | awk '{print $1}')
done

mv rPeaks ../tmp.bed
cd ../
rm -r scratch

sort -k1,1 -k2,2n tmp.bed > sorted.dnase.bed

awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' sorted.dnase.bed> $genome-rDNase-Filtered.bed

## Process the ATAC data
echo -e "Combining ATAC peaks..."
cd $dir/Processed-CARs
cat output.ENCFF668GUI output.ENCFF033RPN output.ENCFF811ERB output.ENCFF509NRA output.ENCFF057BIJ output.ENCFF961XDO \
output.ENCFF784HME output.ENCFF049WJI output.ENCFF796DRU > tmp
mv tmp $dir/$genome-ATAC-All.bed

if [ $genome == "mm10" ]
then
#cutoff=0.1454
atac_cutoff=0.158899 #2021 cutoff
elif [ $genome == "hg38" ]
then
atac_cutoff=2.2947 #custom cutoff

fi

echo -e "Filtering ATAC peaks..."
cd $dir
awk '{if ($1 !~ /_/ && $6 >= '$atac_cutoff') print $0}' \
    $genome-ATAC-All.bed | grep -v "chrEBV" | grep -v "chrM" > $genome-ATAC-Filtered.bed

mkdir scratch
cp $genome-ATAC-Filtered.bed scratch/tmp.bed
cd scratch

echo -e "Sorting ATAC peaks..."
sort -k1,1 -k2,2n tmp.bed > sorted
rm -f rPeaks
num=$(wc -l sorted | awk '{print $1}')

echo -e "Merging ATAC peaks..."
while [ $num -gt 0 ]
do
    echo -e "\t" $num
    bedtools merge -i sorted -c 4,6 -o collapse,collapse > merge
    python $scriptDir/pick.best.peak.py merge > peak-list
    awk 'FNR==NR {x[$1];next} ($4 in x)' peak-list sorted >> rPeaks
    bedtools intersect -v -a sorted -b rPeaks > remaining
    mv remaining sorted
    num=$(wc -l sorted | awk '{print $1}')
done

mv rPeaks ../tmp.bed
cd ../
rm -r scratch

sort -k1,1 -k2,2n tmp.bed > sorted.atac.bed

awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' sorted.atac.bed > $genome-rATAC-Filtered.bed

## Process the single cell ATAC data
echo -e "Combining single cell ATAC peaks..."
cat output.sc.Transverse_ENCSR349XKD output.sc.Transverse_ENCSR434SXE output.sc.Transverse_ENCSR506YMX \
output.sc.Transverse_ENCSR997YNO output.sc.Left_ENCSR830FPR  > tmp
mv tmp $dir/$genome-scATAC-All.bed

if [ $genome == "mm10" ]
then
scatac_cutoff=0
elif [ $genome == "hg38" ]
then
scatac_cutoff=0.00688278
fi

echo -e "Filtering single cell ATAC peaks..."
cd $dir
awk '{if ($1 !~ /_/ && $6 >= '$scatac_cutoff') print $0}' \
    $genome-scATAC-All.bed | grep -v "chrEBV" | grep -v "chrM" > $genome-scATAC-Filtered.bed

mkdir scratch
cp $genome-scATAC-Filtered.bed scratch/tmp.bed
cd scratch

echo -e "Sorting single cell ATAC peaks..."
sort -k1,1 -k2,2n tmp.bed > sorted
rm -f rPeaks
num=$(wc -l sorted | awk '{print $1}')

echo -e "Merging single cell ATAC peaks..."
while [ $num -gt 0 ]
do
    echo -e "\t" $num
    bedtools merge -i sorted -c 4,6 -o collapse,collapse > merge
    python $scriptDir/pick.best.peak.py merge > peak-list
    awk 'FNR==NR {x[$1];next} ($4 in x)' peak-list sorted >> rPeaks
    bedtools intersect -v -a sorted -b rPeaks > remaining
    mv remaining sorted
    num=$(wc -l sorted | awk '{print $1}')
done

mv rPeaks ../tmp.bed
cd ../
rm -r scratch

sort -k1,1 -k2,2n tmp.bed > sorted.scatac.bed

awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' sorted.scatac.bed > $genome-rscATAC-Filtered.bed

####### Combine DNase and ATAC peaks ############

cat $genome-rDNase-Filtered.bed $genome-rATAC-Filtered.bed $genome-rscATAC-Filtered.bed | sort -k1,1 -k2,2n | mergeBed -i stdin > $genome-rCARs-v3.bed

# Add identifier label to the cCREs
awk '{printf "%s\t%s\t%s\tColoncCRE%d\n", $1, $2, $3, NR}' "$genome-rCARs-v3.bed" > "$genome-rCARs-v3-labelled.bed"

