#!/usr/bin/env python3

import pyBigWig
import subprocess

umap_bw = "/Users/ryanhagan/swedish_paper/k36.Umap.MultiTrackMappability.bw"
encode_blacklist = "/Users/ryanhagan/swedish_paper/DriverPower/ENCFF356LFX.bed.gz"
output_bed = "/Users/ryanhagan/swedish_paper/DriverPower/callable.bed"

# Define Chrs and create tmp bed
standard_chroms = [f"chr{i}" for i in range(1,22)] + ["chrX", "chrY"]
temp_bed = output_bed + ".tmp.bed"

# Open BigWig
bw = pyBigWig.open(umap_bw)

# Stream uniquely mappable intervals to BED file
with open(temp_bed, 'w') as out:
    for chrom in bw.chroms():
        if chrom not in standard_chroms:
            continue
        for start, end, score in bw.intervals(chrom) or []:
            if score == 1:
                out.write(f"{chrom}\t{start}\t{end}\n")

bw.close()

# Sort and merge intervals using bedtools
merged_bed = temp_bed + ".merged.bed"
subprocess.run([
    "bedtools", "sort", "-i", temp_bed
], stdout=open(temp_bed + ".sorted", "w"))
subprocess.run([
    "bedtools", "merge", "-i", temp_bed + ".sorted"
], stdout=open(merged_bed, "w"))

# Subtract ENCODE blacklist
subprocess.run([
    "bedtools", "subtract", "-A",
    "-a", merged_bed,
    "-b", encode_blacklist
], stdout=open(output_bed, "w"))

# Cleanup temp files
for f in [temp_bed, temp_bed + ".sorted", merged_bed]:
    subprocess.run(["rm", f])

print(f"Callable BED written to: {output_bed}")
