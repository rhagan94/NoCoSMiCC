#!/bin/bash

scriptDir=/project/home/p201120/ryan/cCRE_pipeline/scripts
signalDir=/project/home/p201120/ryan/cCRE_pipeline/ENCODE_outputs/signal-output
outDir=/project/home/p201120/ryan/cCRE_pipeline/ENCODE_outputs

for mode in scATAC H3K4me3 H3K27ac CTCF; do
    echo "Determining maxZ for $mode..."
    cd "$signalDir/$mode"
    "$PYTHON" "$scriptDir/select-max-zscore.py" *.txt \
        > "$outDir/$mode-maxZ.txt"
    echo "  Done: $outDir/$mode-maxZ.txt"
done

echo "All complete."

# Verify
wc -l "$outDir"/*-maxZ.txt
