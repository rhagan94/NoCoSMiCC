## This script is used to identify cCREs lifted over to the mouse and dog genomes ##

# Define paths
filesDir=/project/home/p201120/ryan/cCRE_pipeline/files
liftOverDir=/project/home/p201120/ryan/cCRE_pipeline/outputs/LiftOver
cCREs=/project/home/p201120/ryan/cCRE_pipeline/ENCODE_outputs/colon-epithelial-cCREs-Z1.28.bed

# Download liftOver files
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/vsHg38/reciprocalBest/hg38.mm10.rbest.chain.gz -P $filesDir
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsCanFam6/reciprocalBest/hg38.canFam6.rbest.chain.gz -P $filesDir

# Unzip
gunzip $filesDir/*rbest.chain.gz

# liftover cCREs to mouse
# 0.5
liftOver $cCREs $filesDir/hg38.mm10.rbest.chain $liftOverDir/sc_cCRE_mouse_lifted.bed $liftOverDir/sc_cCRE_mouse_unlifted.bed -minMatch=0.5
# 0.9
liftOver $cCREs $filesDir/hg38.mm10.rbest.chain $liftOverDir/sc_cCRE_mouse_lifted_0.9.bed $liftOverDir/sc_cCRE_mouse_unlifted_0.9.bed -minMatch=0.9

# liftover cCREs to dog
# 0.5
liftOver $cCREs $filesDir/hg38.canFam6.rbest.chain $liftOverDir/sc_cCRE_dog_lifted.bed $liftOverDir/sc_cCRE_dog_unlifted.bed -minMatch=0.5
# 0.9
liftOver $cCREs $filesDir/hg38.canFam6.rbest.chain $liftOverDir/sc_cCRE_dog_lifted_0.9.bed $liftOverDir/sc_cCRE_dog_unlifted_0.9.bed -minMatch=0.9

# Filter to keep only those lifted over to BOTH mouse and dog
# Cut
cut -f4 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_mouse_lifted.bed \
  > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids.txt

cut -f4 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_dog_lifted.bed \
  > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids.txt

# Sort 
sort /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids.txt \
  -o /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids_sorted.txt

sort /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids.txt \
  -o /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids_sorted.txt

# Get the common IDs
comm -12 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids_sorted.txt \
         /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids_sorted.txt \
  > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/overlapping_cCRE_ids.txt

# Define the files
lifted_cCREs=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/overlapping_cCRE_ids.txt
raw_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs-2025.bed
filtered_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs-2025.bed

# Filter
awk 'NR==FNR {gsub(/"/, "", $1); lifted_cCREs[$1]; next} $4 in lifted_cCREs' $lifted_cCREs $raw_cCRE_bed > $filtered_cCRE_bed

## END OF SCRIPT ##
