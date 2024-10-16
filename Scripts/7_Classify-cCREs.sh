genome=hg38

dataDir=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs
rdhs=$dataDir/sc-only-rCARs-labelled.bed
scriptDir=/Users/ryanhagan/NoCoSMiCC/Scripts/ENCODE_cCRE_pipeline_replication

cd $dataDir/MaxZ

# Sort the TSS bed file
bedtools sort -i /Users/ryanhagan/NoCoSMiCC/Files/hg38_tss.bed > /Users/ryanhagan/NoCoSMiCC/Files/hg38_tss_sorted.bed

# Create the necessary TSS files for annotation
bedtools slop -i /Users/ryanhagan/NoCoSMiCC/Files/hg38_tss_sorted.bed -g /Users/ryanhagan/NoCoSMiCC/Files/hg38.chrom.sizes.txt -b 2000 > /Users/ryanhagan/NoCoSMiCC/Files/Hg38_TSS.Basic.4k.bed

# remove chrM
grep -v "chrM" /Users/ryanhagan/NoCoSMiCC/Files/Hg38_TSS.Basic.4k.bed > /Users/ryanhagan/NoCoSMiCC/Files/Hg38_TSS.Basic.4k_no_chrM.bed

if [[ $genome == "mm10" ]]
then
    prox=/Users/ryanhagan/NoCoSMiCC/Files/Hg38_TSS.Basic.4k_no_chrM.bed
    tss=/Users/ryanhagan/NoCoSMiCC/Files/hg38_tss_sorted.bed
    ChromInfo=~/Lab/Reference/Mouse/ChromInfo.txt
elif [[ $genome == "hg38" ]]
then
    prox=/Users/ryanhagan/NoCoSMiCC/Files/Hg38_TSS.Basic.4k_no_chrM.bed
    tss=/Users/ryanhagan/NoCoSMiCC/Files/hg38_tss_sorted.bed
    #ChromInfo=~/Lab/Reference/Human/hg38/chromInfo.txt
fi

echo "Splitting cCREs into groups..."
sed 's/\t$//' $rdhs > rdhs-cleaned.bed
awk '{if ($2 > 1.28) print $0}' $genome-scATAC-maxZ.txt > list
awk 'FNR==NR {x[$1];next} ($4 in x)' list rdhs-cleaned.bed > bed
bedtools intersect -u -a bed -b $tss > tss
bedtools intersect -v -a bed -b $tss > a1
bedtools intersect -u -a a1 -b $prox | sort -k1,1 -k2,2n > prox
bedtools intersect -v -a bed -b $prox > distal

bedtools closest -d -a prox -b $tss > tmp
python $scriptDir/calculate-center-distance.py tmp agnostic > new
awk '{if ($2 >= -200 && $2 <= 200) print $0}' new > center-distance
awk '{if ($2 < -2000 || $2 > 2000) print $0}' new > far
awk 'FNR==NR {x[$1];next} ($4 in x)' center-distance prox >> tss
awk 'FNR==NR {x[$1];next} ($4 in x)' far prox >> distal
cat center-distance far > new
awk 'FNR==NR {x[$1];next} !($4 in x)' new prox > tmp
mv tmp prox
rm new


###TSS elements###
awk 'FNR==NR {x[$4];next} ($1 in x)' tss $genome-H3K4me3-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > PLS

awk 'FNR==NR {x[$4];next} ($1 in x)' tss $genome-H3K4me3-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no1

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 $genome-H3K27ac-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > pELS

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 $genome-H3K27ac-maxZ.txt | \
    awk '{if ($2 <= 1.28) xprint $0}' > no2

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 $genome-CTCF-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > CTCFonly

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 $genome-CTCF-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > CAonly


###Proximal elements###
awk 'FNR==NR {x[$4];next} ($1 in x)' prox $genome-H3K27ac-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' >> pELS

awk 'FNR==NR {x[$4];next} ($1 in x)' prox $genome-H3K27ac-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no1

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 $genome-H3K4me3-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > DNaseK4

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 $genome-H3K4me3-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no2

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 $genome-CTCF-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' >> CTCFonly

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 $genome-CTCF-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' >> CAonly


###Distal elements###
awk 'FNR==NR {x[$4];next} ($1 in x)' distal $genome-H3K27ac-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > dELS

awk 'FNR==NR {x[$4];next} ($1 in x)' distal $genome-H3K27ac-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no1

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 $genome-H3K4me3-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' >> DNaseK4

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 $genome-H3K4me3-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no2

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 $genome-CTCF-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' >> CTCFonly

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 $genome-CTCF-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' >> CAonly


echo "Accessioning ccREs..."
awk 'FNR==NR {x[$1];next} ($4 in x)' PLS $rdhs | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "PLS" }' > l.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' pELS $rdhs | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "pELS"}' >> l.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' dELS $rdhs | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "dELS"}' >> l.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' DNaseK4 $rdhs | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "DNase-H3K4me3"}' >> l.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' CTCFonly $rdhs | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "CTCF-only"}' >> l.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' CAonly $rdhs | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" "CA-only"}' >> l.bed


awk 'FNR==NR {x[$4];next} ($1 in x)' l.bed $genome-CTCF-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > CTCFall

awk 'FNR==NR {x[$1];next} ($4 in x)' CTCFall l.bed | \
    awk '{print $0",CTCF-bound"}' > m.bed
awk 'FNR==NR {x[$1];next} !($4 in x)' CTCFall l.bed | \
    awk '{print $0}' >> m.bed

sort -k1,1 -k2,2n m.bed > l.bed

mv l.bed sc-only-cCREs.bed

### END OF PIPELINE
