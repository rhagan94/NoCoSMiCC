#!/bin/bash
ONCOKB_DIR="/mnt/tier2/project/p201120/ryan/oncokb-annotator"
REF="refs/hg38.analysisSet.fa"
VEP_PATH="/mnt/tier2/project/p201120/ryan/envs/get/bin"
VEP_CACHE="/project/home/p201120/ryan/VEP/vep_cache"
ONCOKB_TOKEN="ONCOKB_API_TOKEN"
OUTDIR="oncokb_output"

mkdir -p $OUTDIR

for COHORT in Hartwig ICGC_DACO Mutographs Nunes; do
    echo "=== Processing $COHORT ==="
    COHORT_DIR="$OUTDIR/$COHORT"
    mkdir -p $COHORT_DIR
    rm -f ${COHORT_DIR}/*.maf ${COHORT_DIR}/*.vcf

    for vcf in VCFs_filtered/${COHORT}/*.vcf.gz; do
        sample=$(basename $vcf .filtered.vcf.gz)
        echo "  $sample..."
        bcftools view -R target_genes.bed $vcf 2>/dev/null > ${COHORT_DIR}/${sample}.filtered.vcf
        perl ${ONCOKB_DIR}/vcf2maf.pl \
            --input-vcf  ${COHORT_DIR}/${sample}.filtered.vcf \
            --output-maf ${COHORT_DIR}/${sample}.maf \
            --tumor-id   $sample \
            --ncbi-build GRCh38 \
            --ref-fasta  $REF \
            --vep-path   $VEP_PATH \
            --vep-data   $VEP_CACHE \
            --custom-enst oncokb_isoforms.txt
    done

    # Combine MAFs
    rm -f ${COHORT_DIR}/${COHORT}.maf
    sed -n '2p' $(ls ${COHORT_DIR}/*.maf | head -1) > ${COHORT_DIR}/${COHORT}.maf
    for maf in ${COHORT_DIR}/*.maf; do
        [[ $maf == *"${COHORT}.maf" ]] && continue
        tail -n +3 $maf
    done >> ${COHORT_DIR}/${COHORT}.maf

    echo "  $(wc -l < ${COHORT_DIR}/${COHORT}.maf) lines in combined MAF"

    # Run OncoKB
    python ${ONCOKB_DIR}/MafAnnotator.py \
        -i ${COHORT_DIR}/${COHORT}.maf \
        -o ${COHORT_DIR}/${COHORT}.oncokb.maf \
        -b $ONCOKB_TOKEN \
        -q HGVSp_Short

    echo "=== Done: $COHORT ==="
done
