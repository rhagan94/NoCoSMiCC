## Script to download latest GnomAD variants and filter VCFs

export OUTDIR="/project/home/p201120/ryan/gnomAD"

parallel -j 24 --eta \
    'echo "Processing chr{}..." && \
     bcftools view \
         --threads 4 \
         -i "AF > 0.01" \
         "gs://gcp-public-data--gnomad/release/4.1.1/vcf/genomes/gnomad.genomes.v4.1.1.sites.chr{}.vcf.bgz" \
         -O z -o ${OUTDIR}/gnomad_AF0.01.chr{}.vcf.gz && \
     bcftools index --tbi ${OUTDIR}/gnomad_AF0.01.chr{}.vcf.gz' \
    ::: {1..22} X Y

# merge 
bcftools concat \
    ${OUTDIR}/gnomad_AF0.01.chr{1..22}.vcf.gz \
    ${OUTDIR}/gnomad_AF0.01.chrX.vcf.gz \
    ${OUTDIR}/gnomad_AF0.01.chrY.vcf.gz \
    -O z -o ${OUTDIR}/gnomad_AF0.01.all.vcf.gz

# index
bcftools index --tbi ${OUTDIR}/gnomad_AF0.01.all.vcf.gz

mkdir -p VCFs_filtered/{Hartwig,ICGC_DACO,Nunes,Mutographs}

GNOMAD="/tmp/gnomad_AF0.01.sitesonly.vcf.gz"

mkdir -p VCFs_filtered/{Hartwig,ICGC_DACO,Nunes,Mutographs}

# Build jobs list
> jobs.txt
for f in VCFs/Hartwig_PASS/*.vcf.gz; do echo "$f|Hartwig|.purple.somatic.pass.vcf.gz"  >> jobs.txt; done
for f in VCFs/ICGC_DACO/*.vcf.gz;    do echo "$f|ICGC_DACO|.vcf.gz"                    >> jobs.txt; done
for f in VCFs/Nunes/*.vcf.gz;        do echo "$f|Nunes|.vcf.gz"                         >> jobs.txt; done
for f in VCFs/Mutographs/*.vcf.gz;   do echo "$f|Mutographs|.vcf.gz"                    >> jobs.txt; done

echo "Filtering $(wc -l < jobs.txt) VCFs..."

# Filter: remove any variant with exact match in gnomAD AF > 1% blacklist
parallel -j 64 --colsep '\|' --bar \
    'sample=$(basename {1} {3}); \
     bcftools isec -C {1} '"$GNOMAD"' -w1 --threads 2 -O z \
         -o VCFs_filtered/{2}/${sample}.filtered.vcf.gz && \
     bcftools index --tbi VCFs_filtered/{2}/${sample}.filtered.vcf.gz' \
    :::: jobs.txt

echo "Done."
