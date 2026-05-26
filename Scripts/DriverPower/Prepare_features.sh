# Define elements
elements=1kb_elements_primary.sorted.bed6

# make coverage table and nucleotide content table
./generate_cg.sh ./hg38.primary.fa $elements k100.callable_final.faidxsorted.bed.gz ./tmp/example.cg ./tmp/example.totcg

# make nucleotide content features
$PYTHON generate_nuc_covar.py ./tmp/example.cg ./output/colon_nuc_features.tsv

# make bigwig features
./bws2cv.sh ./bigWigs/ $elements ./output/bw_features.tsv 4

## STEP 4: combine all features
./combine_cv.sh ./output ./1kb_gwide_features.tsv

# Run prepare.py to prepare response tables
$PYTHON prepare.py MSS_nHM_filt_mutations.tsv.gz ./BMR_features/colon-epithelial-cCREs.bed4 ./BMR_features/k100.callable_final.faidxsorted.bed.gz cCREs_y.tsv
$PYTHON prepare.py MSS_nHM_filt_mutations.tsv.gz ./BMR_features/1kb_elements_primary.sorted.bed4 ./BMR_features/k100.callable_final.faidxsorted.bed.gz 1kb_y.tsv
$PYTHON prepare.py MSS_nHM_filt_mutations.tsv.gz ./BMR_features/5kb_elements_primary.sorted.bed4 ./BMR_features/k100.callable_final.faidxsorted.bed.gz 5kb_y.tsv

# split by chr arm
$PYTHON split_cCRE_tsv_by_arm.py
$PYTHON split_1kb_tsv_by_arm.py
$PYTHON split_5kb_tsv_by_arm.py

# convert TSV to HDF5 format
$PYTHON convert_tsv_hdf5.py
