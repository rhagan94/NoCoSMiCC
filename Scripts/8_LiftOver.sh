## This script is used to identify cCREs lifted over to the mouse and dog genomes ##

# First, obtain the necessary liftOver files
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/vsHg38/reciprocalBest/hg38.mm10.rbest.chain.gz -P /Users/ryanhagan/NoCoSMiCC/Files
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsCanFam6/reciprocalBest/canFam6.hg38.rbest.chain.gz -P /Users/ryanhagan/NoCoSMiCC/Files

# Unzip
gunzip /Users/ryanhagan/NoCoSMiCC/Files/*rbest.chain.gz

# liftover human cCREs to mouse - these are mouse co-ordinates
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs-2025.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm10.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_mouse_lifted.bed \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_mouse_unlifted.bed -minMatch=0.5

# liftover human cCREs to dog - these are dog co-ordinates
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs-2025.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.canFam6.rbest.chain  /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_dog_lifted.bed \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_dog_unlifted.bed -minMatch=0.5

# Now filter the human cCREs and keep only those lifted over to BOTH mouse and dog
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
