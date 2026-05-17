# Clone the OncoKB annotator
git clone https://github.com/oncokb/oncokb-annotator.git
cd oncokb-annotator
pip install -r requirements/common.txt

# Download vcf2maf
wget https://raw.githubusercontent.com/mskcc/vcf2maf/v1.6.22/vcf2maf.pl

# Download OncoKB canonical transcripts
curl -s "https://www.oncokb.org/api/v1/utils/allCuratedGenes" \
    -H "Authorization: Bearer YOUR_ONCOKB_API_TOKEN" \
    > oncokb_cancer_genes.json

python3 -c "
import json
with open('oncokb_cancer_genes.json') as f:
    genes = json.load(f)
for g in genes:
    enst = g.get('grch38Isoform', '')
    symbol = g.get('hugoSymbol', '')
    if enst and symbol:
        print(f'{symbol}\t{enst}')
" > oncokb_isoforms.txt
