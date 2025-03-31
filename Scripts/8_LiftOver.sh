# liftover human cCREs to mouse
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v4.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm10.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_lifted \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_UNlifted -minMatch=0.5

###############################################################
# liftover human cCREs to mouse for single cell only cCREs
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm10.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_mouse_lifted \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_mouse_UNlifted -minMatch=0.9

# Extract fourth column from the lifted cCREs in mouse
cut -f 4 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_mouse_lifted > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids.txt

# Sort the ID files
sort /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids.txt -o /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids_sorted.txt

lifted_cCREs=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids_sorted.txt
raw_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs.bed
filtered_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs.bed
#filtered_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs-0_9.bed

awk 'NR==FNR {gsub(/"/, "", $1); lifted_cCREs[$1]; next} $4 in lifted_cCREs' $lifted_cCREs $raw_cCRE_bed > $filtered_cCRE_bed

# Find the lifted cCREs located within syntenic blocks
bedtools intersect -u -f 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/sc-lifted-syntenic-cCREs.bed

# Find the lifted and syntenic cCREs located within compartments A and I
bedtools intersect -wa -u -f 1 -a /Users/ryanhagan/NoCoSMiCC/synteny/sc-lifted-syntenic-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/compartment_A.bed \
> /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-cCREs.bed 
wc -l /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-cCREs.bed

bedtools intersect -wa -u -f 1 -a /Users/ryanhagan/NoCoSMiCC/synteny/sc-lifted-syntenic-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/compartment_I.bed \
> /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-cCREs.bed 
wc -l /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-cCREs.bed

# Fing the lifted, syntenic and compartment-overlapping cCREs that also overlap Zoonomia elements
bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-G1-cCREs.bed 
wc -l /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-G1-cCREs.bed

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G2_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-G2-cCREs.bed 
wc -l /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-G2-cCREs.bed

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G3_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-G3-cCREs.bed 
wc -l /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compA-G3-cCREs.bed 

##

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-G1-cCREs.bed 
wc -l /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-G1-cCREs.bed

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G2_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-G2-cCREs.bed 
wc -l /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-G2-cCREs.bed

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G3_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-G3-cCREs.bed 
wc -l /Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-G3-cCREs.bed




### Alternative filtering

## lifted and Zoonomia elements
bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/sc-lifted-G1-cCREs.bed

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G2_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/sc-lifted-G2-cCREs.bed

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G3_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/sc-lifted-G3-cCREs.bed

## Zoonomia elements only
bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/sc-G1-cCREs.bed

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G2_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/sc-G2-cCREs.bed

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G3_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/sc-G3-cCREs.bed

## Compartments only
bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/compartment_A.bed \
> /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/sc-CompA-cCREs.bed

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/compartment_I.bed \
> /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/sc-CompI-cCREs.bed

bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/compartment_B.bed \
> /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/sc-CompB-cCREs.bed

## Synteny only
bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/sc-synteny-cCREs.bed

bedtools intersect -wa -v -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/sc-no-synteny-cCREs.bed

## Liftover only
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs.bed
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs-0_9.bed


#### Bulk analysis

# LiftOver for mouse

liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/bulk-only-rCARs/MaxZ/hg38-cCREs-bulk-only.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm10.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/bulk_cCREs_lifted.bed \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/bulk_cCREs_UNlifted.bed -minMatch=0.5

# Extract fourth column from the lifted cCREs in mouse
cut -f 4 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/bulk_cCREs_lifted.bed > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids.txt

# Sort the ID files
sort /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids.txt -o /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids_sorted.txt

lifted_cCREs=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids_sorted.txt
raw_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/bulk-only-rCARs/MaxZ/hg38-cCREs-bulk-only.bed
filtered_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/bulk-only-lifted-cCREs.bed

awk 'NR==FNR {gsub(/"/, "", $1); lifted_cCREs[$1]; next} $4 in lifted_cCREs' $lifted_cCREs $raw_cCRE_bed > $filtered_cCRE_bed

## lifted and Zoonomia elements
bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/bulk-only-lifted-cCREs.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/bulk-lifted-G1-cCREs.bed
wc -l /Users/ryanhagan/NoCoSMiCC/synteny/bulk-lifted-G1-cCREs.bed

##############################################################







# liftover human cCREs to dog
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v4.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.canFam6.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_lifted \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_UNlifted -minMatch=0.5

# liftover human cCREs to chimpanzee
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v4.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.panTro6.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_chimp_lifted \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_chimp_UNlifted -minMatch=0.5


bedtools intersect -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_lifted -b /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_lifted \
> /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_common_lifted.bed


### Adjust the minMatch to 0.9

# liftover human cCREs to mouse
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm10.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_lifted_strict \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_UNlifted -minMatch=0.9

# liftover human cCREs to dog
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.canFam6.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_lifted \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_UNlifted -minMatch=0.9

# liftover human cCREs to chimpanzee
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.panTro6.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_chimp_lifted \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_chimp_UNlifted -minMatch=0.9


# using mm39 instead of mm10
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm39.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dolphin_lifted \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dolphin_UNlifted -minMatch=0.9



# Overlap the syntenic regions from GENESPACE with the cCREs

# Find common syntenic regions between human_mouse and human_dog
bedtools intersect -a /Users/ryanhagan/NoCoSMiCC/synteny/human_mouse_syntenic_blocks.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/human_dog_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed

# Assess the cCREs that have complete overlap in a syntenic region - syntenic cCREs
bedtools intersect -u -f 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/syntenic_cCREs.bed

# Assess the cCREs that have do not have complete overlap in a syntenic region - nonsyntenic cCREs
bedtools intersect -v -f 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/nonsyntenic_cCREs.bed




### Filter the cCRE bed file by retaining only cCREs that have a score > 0.5 in BOTH mouse and dog 

# Extract fourth column from the lifted cCREs in mouse
cut -f 4 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_lifted > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids.txt

# Extract fourth column from the lifted cCREs in dog
cut -f 4 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_lifted > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids.txt

# Sort the ID files
sort /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids.txt -o /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids_sorted.txt
sort /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids.txt -o /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids_sorted.txt

# Find common elements
comm -12 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids_sorted.txt /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids_sorted.txt > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_common_ids.txt


common_lifted_cCREs=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_common_ids.txt
raw_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v4.bed
filtered_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/common-liftover-cCREs.bed

awk 'NR==FNR {gsub(/"/, "", $1); common_lifted_cCREs[$1]; next} $4 in common_lifted_cCREs' $common_lifted_cCREs $raw_cCRE_bed > $filtered_cCRE_bed











# intersect human_mouse and human_dog syntenic regions 
bedtools intersect -e -f 1 -F 1 -a /Users/ryanhagan/NoCoSMiCC/synteny/human_mouse_syntenic_blocks.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/human_dog_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed




bedtools intersect -f 1 -a /Users/ryanhagan/NoCoSMiCC/synteny/human_mouse_syntenic_blocks.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/human_dog_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/test1_common_syntenic_blocks.bed

bedtools intersect -F 1 -a /Users/ryanhagan/NoCoSMiCC/synteny/human_mouse_syntenic_blocks.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/human_dog_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/test2_common_syntenic_blocks.bed




# Overlap the common syntenic regions with the cCREs
bedtools intersect -u -f 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/syntenic_cCREs.bed

bedtools intersect -v -f 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/nonsyntenic_cCREs.bed





### test 2

bedtools intersect -a /Users/ryanhagan/NoCoSMiCC/synteny/human_mouse_syntenic_blocks.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed \
 > /Users/ryanhagan/NoCoSMiCC/synteny/test3_common_syntenic_blocks.bed

bedtools intersect -a /Users/ryanhagan/NoCoSMiCC/synteny/common_syntenic_blocks.bed -b /Users/ryanhagan/NoCoSMiCC/synteny/human_mouse_syntenic_blocks.bed \
 > /Users/ryanhagan/NoCoSMiCC/synteny/test4_common_syntenic_blocks.bed






bedtools intersect -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/zooUCEs.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/UltraConserved_colon_cCREs.bed



bedtools intersect -a /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/zooUCEs.bed -b /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed \
> /Users/ryanhagan/NoCoSMiCC/synteny/UltraConserved_ENCODE_regions_in_colon.bed


bedtools intersect -F 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/Zoonomia/Conserved_colon_G1_cCREs.bed

bedtools intersect -F 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G2_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/Zoonomia/Evolving_colon_G2_cCREs.bed

bedtools intersect -F 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G3_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/Zoonomia/Primate_colon_G3_cCREs.bed


bedtools intersect -F 0.9 -f 0.9 -e -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed \
> /Users/ryanhagan/NoCoSMiCC/Zoonomia/Conserved_colon_G1_cCREs.bed

## TAD boundary overlaps

bedtools intersect -f 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/Files/Q4of4_40kbBoundaries_removeGaps.bed \
> /Users/ryanhagan/NoCoSMiCC/conserved_TAD_overlaps.bed

bedtools intersect -f 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/Files/Q4of4_100kbBookendBoundaries_removeGaps.bed \
> /Users/ryanhagan/NoCoSMiCC/conserved_TAD_overlaps.bed

wc -l /Users/ryanhagan/NoCoSMiCC/conserved_TAD_overlaps.bed


bedtools intersect -v -f 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed -b /Users/ryanhagan/NoCoSMiCC/Files/Q4of4_100kbBookendBoundaries_removeGaps.bed \
> /Users/ryanhagan/NoCoSMiCC/conserved_TAD_no_overlaps.bed

wc -l /Users/ryanhagan/NoCoSMiCC/conserved_TAD_no_overlaps.bed






####### Constrained final list of cCREs for plotting

# liftover human cCREs to mouse
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/human_nc_ccre_G1.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm10.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/nc_cCRE_mouse_G1_lifted.bed \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/nc_cCRE_mouse__G1_UNlifted.bed -minMatch=0.5

liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/human_nc_ccre_nonG1.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm10.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/nc_cCRE_mouse_nonG1_lifted.bed \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/nc_cCRE_mouse__G1_nonG1_UNlifted.bed -minMatch=0.5



##### Single cell only cCREs - found in compartment A
# liftover human cCREs to mouse
liftOver /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/single_cell_compartment_A_cCREs.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm10.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_compA_lifted.bed \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_compA_UNlifted.bed -minMatch=0.5

## How many are located within syntenic blocks?
bedtools intersect -wa -u -f 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_compA_lifted.bed -b /Users/ryanhagan/NoCoSMiCC/Synteny/common_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/Synteny/sc_compA_lift_synt.bed
wc -l /Users/ryanhagan/NoCoSMiCC/Synteny/sc_compA_lift_synt.bed

##### Bulk only cCREs - found in compartment A
liftOver /Users/ryanhagan/NoCoSMiCC/HiC_Analysis/bulk_compartment_A_cCREs.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm10.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/bulk_compA_lifted.bed \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/bulk_compA_UNlifted.bed -minMatch=0.5

## How many are located within syntenic blocks?
bedtools intersect -wa -u -f 1 -a /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/bulk_compA_lifted.bed -b /Users/ryanhagan/NoCoSMiCC/Synteny/common_syntenic_blocks.bed \
> /Users/ryanhagan/NoCoSMiCC/Synteny/bulk_compA_lift_synt.bed
wc -l /Users/ryanhagan/NoCoSMiCC/Synteny/bulk_compA_lift_synt.bed


## Additional filtering - overlap with G1 cCREs - single cell
# Overlapping G1 (highly-conserved) cCREs
bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/Synteny/sc_compA_lift_synt.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed > /Users/ryanhagan/NoCoSMiCC/Zoonomia/sc_liftover_compA_synt_G1_cCREs.bed 

# count cCREs 
wc -l /Users/ryanhagan/NoCoSMiCC/Zoonomia/sc_liftover_compA_synt_G1_cCREs.bed

# Not overlapping G1 (highly-conserved) cCREs
bedtools intersect -wa -v -a /Users/ryanhagan/NoCoSMiCC/Synteny/sc_compA_lift_synt.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed > /Users/ryanhagan/NoCoSMiCC/Zoonomia/sc_liftover_compA_synt_nonG1_cCREs.bed 

# count cCREs 
wc -l /Users/ryanhagan/NoCoSMiCC/Zoonomia/sc_liftover_compA_synt_nonG1_cCREs.bed


## Additional filtering - overlap with G1 cCREs - bulk
# Overlapping G1 (highly-conserved) cCREs
bedtools intersect -wa -u -a /Users/ryanhagan/NoCoSMiCC/Synteny/bulk_compA_lift_synt.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed > /Users/ryanhagan/NoCoSMiCC/Zoonomia/bulk_liftover_compA_synt_G1_cCREs.bed 

# count cCREs 
wc -l /Users/ryanhagan/NoCoSMiCC/Zoonomia/bulk_liftover_compA_synt_G1_cCREs.bed

# Not overlapping G1 (highly-conserved) cCREs
bedtools intersect -wa -v -a /Users/ryanhagan/NoCoSMiCC/Synteny/bulk_compA_lift_synt.bed -b /Users/ryanhagan/NoCoSMiCC/Zoonomia/Files/G1_cCREs.bed > /Users/ryanhagan/NoCoSMiCC/Zoonomia/bulk_liftover_compA_synt_nonG1_cCREs.bed 

# count cCREs 
wc -l /Users/ryanhagan/NoCoSMiCC/Zoonomia/bulk_liftover_compA_synt_nonG1_cCREs.bed




### Single cell analysis in 2025

# liftover human cCREs to mouse 
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs-2025.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.mm10.rbest.chain /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_mouse_lifted.bed \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_mouse_unlifted.bed -minMatch=0.5

# liftover human cCREs to dog 
liftOver /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs-2025.bed \
/Users/ryanhagan/NoCoSMiCC/Files/hg38.canFam6.rbest.chain  /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_dog_lifted.bed \
/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_dog_unlifted.bed -minMatch=0.5


# Extract fourth column (ID) from the mouse lifted file
cut -f4 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_mouse_lifted.bed \
  > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids.txt

# Extract fourth column (ID) from the dog lifted file
cut -f4 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_dog_lifted.bed \
  > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids.txt

# Sort both ID files
sort /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids.txt \
  -o /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids_sorted.txt

sort /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids.txt \
  -o /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids_sorted.txt

# Extract the overlapping IDs between mouse and dog
comm -12 /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_ids_sorted.txt \
         /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_ids_sorted.txt \
  > /Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/overlapping_cCRE_ids.txt

lifted_cCREs=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/overlapping_cCRE_ids.txt
raw_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/sc-only-cCREs-2025.bed
filtered_cCRE_bed=/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs-2025.bed

awk 'NR==FNR {gsub(/"/, "", $1); lifted_cCREs[$1]; next} $4 in lifted_cCREs' $lifted_cCREs $raw_cCRE_bed > $filtered_cCRE_bed




