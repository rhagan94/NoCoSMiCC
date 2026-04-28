## This script is used to identify cCREs lifted over to the mouse and dog genomes ##

# Define paths
filesDir=/project/home/p201120/ryan/cCRE_pipeline/files
liftOverDir=/project/home/p201120/ryan/cCRE_pipeline/outputs/LiftOver
cCREs=/project/home/p201120/ryan/cCRE_pipeline/ENCODE_outputs/colon-epithelial-cCREs.bed
lifted_output_file=$liftOverDir/colon-epithelial-lifted-cCREs.bed
non_lifted_output_file=$liftOverDir/colon-epithelial-nonlifted-cCREs.bed

# Download liftOver files and unzip
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm39/vsHg38/reciprocalBest/hg38.mm39.rbest.chain.gz -P $filesDir
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsCanFam6/reciprocalBest/hg38.canFam6.rbest.chain.gz -P $filesDir

gunzip $filesDir/*rbest.chain.gz

# liftover cCREs to mouse
# 0.5
liftOver $cCREs $filesDir/hg38.mm39.rbest.chain $liftOverDir/sc_cCRE_mouse_lifted_0.5.bed $liftOverDir/sc_cCRE_mouse_unlifted_0.5.bed -minMatch=0.5
# 0.9
liftOver $cCREs $filesDir/hg38.mm39.rbest.chain $liftOverDir/sc_cCRE_mouse_lifted_0.9.bed $liftOverDir/sc_cCRE_mouse_unlifted_0.9.bed -minMatch=0.9

# liftover cCREs to dog
# 0.5
liftOver $cCREs $filesDir/hg38.canFam6.rbest.chain $liftOverDir/sc_cCRE_dog_lifted_0.5.bed $liftOverDir/sc_cCRE_dog_unlifted_0.5.bed -minMatch=0.5
# 0.9
liftOver $cCREs $filesDir/hg38.canFam6.rbest.chain $liftOverDir/sc_cCRE_dog_lifted_0.9.bed $liftOverDir/sc_cCRE_dog_unlifted_0.9.bed -minMatch=0.9

# Filter to keep cCREs lifted to mouse and dog
cut -f4 $liftOverDir/sc_cCRE_mouse_lifted_0.5.bed \
  | sort -u -o $liftOverDir/cCRE_mouse_ids_sorted.txt

cut -f4 $liftOverDir/sc_cCRE_dog_lifted_0.5.bed \
  | sort -u -o $liftOverDir/cCRE_dog_ids_sorted.txt

comm -12 $liftOverDir/cCRE_mouse_ids_sorted.txt $liftOverDir/cCRE_dog_ids_sorted.txt > $liftOverDir/overlapping_cCRE_ids.txt

# Lifted output
awk 'NR==FNR {gsub(/"/, "", $1); ids[$1]; next} $4 in ids' "$liftOverDir/overlapping_cCRE_ids.txt" "$cCREs" > "$lifted_output_file"

# Non-lifted output
cut -f4 "$cCREs" | sort -u > $liftOverDir/cCRE_ids_sorted.txt

comm -23 $liftOverDir/cCRE_ids_sorted.txt $liftOverDir/overlapping_cCRE_ids.txt > $liftOverDir/non_overlapping_cCRE_ids.txt

awk 'NR==FNR {ids[$1]; next} $4 in ids' "$liftOverDir/non_overlapping_cCRE_ids.txt" "$cCREs" > "$non_lifted_output_file"
