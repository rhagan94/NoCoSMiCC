BASE="/project/home/p201120/ryan"
VCF_DIR="$BASE/VCFs_filtered/Hartwig"
VEP_CACHE_DIR="$BASE/VEP/vep_cache"
FASTA="$BASE/refs/hg38.analysisSet.fa"
VEP_OUT_DIR="$BASE/VEP/vep_hartwig"
BCFTOOLS="/project/home/p201120/ryan/envs/bcftools/bin/bcftools"
GTF="$BASE/refs/gtf/Homo_sapiens.GRCh38.113.gtf.gz"

mkdir -p "$VEP_OUT_DIR"

# Extract full gene body coordinates from GTF
get_region() {
    local gene="$1"
    zcat "$GTF" | awk -v g="$gene" '
        $3=="gene" {
            for(i=9;i<=NF;i++) {
                if($i=="gene_name") {
                    name=$(i+1); gsub(/[";]/,"",name)
                    if(name==g) print "chr"$1":"$4"-"$5
                }
            }
	}' | head -1
}

echo "Extracting gene coordinates from GTF..."
KRAS_REGION=$(get_region "KRAS")
BRAF_REGION=$(get_region "BRAF")
POLE_REGION=$(get_region "POLE")

echo "  KRAS: $KRAS_REGION"
echo "  BRAF: $BRAF_REGION"
echo "  POLE: $POLE_REGION"

REGIONS="${BRAF_REGION},${KRAS_REGION},${POLE_REGION}"
echo "  Full region string: $REGIONS"

run_vep() {
    local VCF="$1"
    local sample out filtered final

    sample=$(basename "$VCF" .vcf.gz)
    out="$VEP_OUT_DIR/${sample}.vep.vcf"
    filtered="$VEP_OUT_DIR/${sample}.regions.vcf.gz"
    final="$VEP_OUT_DIR/${sample}.canonical_coding.vcf"

    [[ -s "$final" ]] && { echo "Skipping $sample (already done)"; return; }

    echo "Filtering to target regions: $sample"
    "$BCFTOOLS" view "$VCF" -r "$REGIONS" -f PASS -O z -o "$filtered"
    "$BCFTOOLS" index "$filtered"

    if [[ "$("$BCFTOOLS" view -H "$filtered" | wc -l)" -eq 0 ]]; then
        echo "No variants in target regions for $sample — skipping"
        rm -f "$filtered" "${filtered}.csi"
        return
    fi

    echo "Running VEP: $sample"
    vep \
      -i "$filtered" \
      -o "$out" \
      --vcf \
      --symbol \
      --hgvs \
      --numbers \
      --canonical \
      --mane_select \
      --biotype \
      --cache --offline \
      --dir_cache "$VEP_CACHE_DIR" \
      --fasta "$FASTA" \
      --assembly GRCh38 \
      --fork 16 \
      --no_stats \
      --force_overwrite \
      --pick_allele \
      --pick_order canonical,mane_select,appris,biotype,tsl,rank,ccds,length

    echo "Filtering with filter_vep: $sample"
    filter_vep \
      -i "$out" \
      -o "$final" \
      --format vcf \
      --filter "CANONICAL is YES \
        and BIOTYPE is protein_coding \
        and (IMPACT is HIGH or IMPACT is MODERATE)" \
      --force_overwrite

    rm -f "$filtered" "${filtered}.csi" "$out"
    echo "Done: $sample → $final"
}

export -f run_vep
export VEP_CACHE_DIR FASTA VEP_OUT_DIR REGIONS BCFTOOLS

parallel -j 8 run_vep ::: "$VCF_DIR"/*.vcf.gz
