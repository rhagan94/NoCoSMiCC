import pandas as pd
import numpy as np
import glob
import os

oncogenic = {"Oncogenic", "Likely Oncogenic"}
all_results = []

for cohort in ["Hartwig", "ICGC_DACO", "Mutographs", "Nunes"]:
    # Get all sample IDs from VCFs
    samples = [os.path.basename(f).replace(".filtered.vcf.gz", "")
               for f in glob.glob(f"VCFs_filtered/{cohort}/*.vcf.gz")]
    samples_df = pd.DataFrame({"sample_id": samples, "cohort": cohort})

    # Load OncoKB annotated MAF
    maf = pd.read_csv(f"oncokb_output/{cohort}/{cohort}.oncokb.maf",
                      sep="\t", low_memory=False)
    maf_keep = maf[maf["ONCOGENIC"].isin(oncogenic)]

    pivot = maf_keep.groupby(["Tumor_Sample_Barcode", "Hugo_Symbol"])["HGVSp_Short"] \
                    .apply(",".join).unstack()

    for gene in ["BRAF", "KRAS", "POLE"]:
        if gene not in pivot.columns:
            pivot[gene] = None

    pivot = pivot[["BRAF", "KRAS", "POLE"]].fillna("WT").reset_index()
    pivot.rename(columns={"Tumor_Sample_Barcode": "sample_id"}, inplace=True)

    result = samples_df.merge(pivot, on="sample_id", how="left").fillna("WT")
    all_results.append(result)

final = pd.concat(all_results, ignore_index=True)

# Add nMutations
nmut = pd.read_csv("nMutations.tsv", sep="\t", header=None,
                   names=["sample_id", "cohort", "nMutations"])
final = final.merge(nmut[["sample_id", "nMutations"]], on="sample_id", how="left")

# HM classification per cohort (upper quartile + 1.5 * IQR)
def classify_hm(group):
    q75 = group["nMutations"].quantile(0.75)
    q25 = group["nMutations"].quantile(0.25)
    iqr = q75 - q25
    threshold = q75 + 1.5 * iqr
    print(f"  {group['cohort'].iloc[0]}: Q75={q75:.0f}, IQR={iqr:.0f}, threshold={threshold:.0f}")
    group["HM_status"] = np.where(group["nMutations"] > threshold, "HM", "nHM")
    return group

print("HM thresholds per cohort:")
final = final.groupby("cohort", group_keys=False).apply(classify_hm)

# Summary
print(f"\nTotal samples: {len(final)}")
for cohort in ["Hartwig", "ICGC_DACO", "Mutographs", "Nunes"]:
    print(f"\n{cohort}:")
    sub = final[final["cohort"] == cohort]
    for gene in ["BRAF", "KRAS", "POLE"]:
        mut = (sub[gene] != "WT").sum()
        wt  = (sub[gene] == "WT").sum()
        print(f"  {gene}: {mut} mutant, {wt} WT")
    hm  = (sub["HM_status"] == "HM").sum()
    nhm = (sub["HM_status"] == "nHM").sum()
    print(f"  HM: {hm}, nHM: {nhm}")
  
final.to_csv("sample_mutation_status_all_cohorts.tsv", sep="\t", index=False)
