#!/bin/bash

ConsensusPeaks=/mnt/tier2/project/p201120/ryan/cCRE_pipeline/ArchR_Epithelial_Combined/PeakCalls/BED/consensus_peaks.bed
signalDir=/mnt/tier2/project/p201120/ryan/cCRE_pipeline/ArchR_Epithelial_Combined/GroupBigWigs/Sample
outdir=/mnt/tier2/project/p201120/ryan/cCRE_pipeline/outouts/colon-epithelial-rCARs/Processed-CARs
BWAvOverBed=/mnt/tier2/project/p201120/ryan/envs/get/bin/bigWigAverageOverBed

# Make the out dir 
mkdir -p "$outdir"

# Prepare peaks
tmp_peaks=$(mktemp -t peaks_XXXXXX.bed)
awk 'BEGIN{OFS="\t"} {
    print $1, $2, $3, "peak"NR
}' "$ConsensusPeaks" > "$tmp_peaks"

for bw in "$signalDir"/*.bw; do
    dset=$(basename "$bw" -TileSize-100-normMethod-ReadsInTSS-ArchR.bw)
    tmp_out=$(mktemp -t out_XXXXXX.tab)

    "$BWAvOverBed" "$bw" "$tmp_peaks" "$tmp_out"

    # Join
    join -1 4 -2 1 \
        <(sort -k4,4 "$tmp_peaks") \
        <(sort -k1,1 "$tmp_out") \
        | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$1,$5,$9}' \
        > "$outdir/output.sc.$dset"

    rm -f "$tmp_out"
done
