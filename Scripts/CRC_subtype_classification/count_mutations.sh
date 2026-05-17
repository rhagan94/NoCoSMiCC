for COHORT in Hartwig ICGC_DACO Mutographs Nunes; do
    for vcf in VCFs_filtered/${COHORT}/*.vcf.gz; do
        sample=$(basename $vcf .filtered.vcf.gz)
        n=$(bcftools view -H $vcf 2>/dev/null | wc -l)
        echo -e "$sample\t$COHORT\t$n"
    done
done > nMutations.tsv
