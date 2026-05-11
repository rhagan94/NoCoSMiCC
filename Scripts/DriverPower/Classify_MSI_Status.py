import pandas as pd

# Load Nunes cohort activities
# Update this path if your SigProfiler output directory differs
DF_PATH = 'SPA_nunes/output_SBS96/Assignment_Solution/Activities/Assignment_Solution_Activities.txt'
COUNTS_PATH = 'VEP/mutation_counts.tsv'
OUT_PATH = 'SPA_nunes/MSI_POLE_MSS_classification_with_counts.txt'

# Load activities matrix
df = pd.read_csv(
    DF_PATH,
    sep='\t',
    index_col=0
)

print(f"Nunes samples loaded: {len(df)}")

# --------------------------------------------------
# Define signature groups
# --------------------------------------------------
msi_sigs = [
    'SBS6', 'SBS14', 'SBS15',
    'SBS20', 'SBS21', 'SBS26', 'SBS44'
]

pole_sigs = [
    'SBS10a', 'SBS10b'
]

# --------------------------------------------------
# Calculate proportions
# --------------------------------------------------
totals = df.sum(axis=1)

msi_cols = [s for s in msi_sigs if s in df.columns]
pole_cols = [s for s in pole_sigs if s in df.columns]

if len(msi_cols) == 0:
    raise ValueError('No MSI signatures found in activities table.')

if len(pole_cols) == 0:
    print('Warning: No POLE signatures found in activities table.')

msi_prop = df[msi_cols].sum(axis=1) / totals

if len(pole_cols) > 0:
    pole_prop = df[pole_cols].sum(axis=1) / totals
else:
    pole_prop = pd.Series(0, index=df.index)

# --------------------------------------------------
# Classification logic
# --------------------------------------------------
def classify(idx):
    if pole_prop[idx] > 0.20:
        return 'POLE'
    elif msi_prop[idx] > 0.30:
        return 'MSI'
    else:
	return 'MSS'


df['Status'] = [classify(i) for i in df.index]
df['MSI_proportion'] = msi_prop.round(4)
df['POLE_proportion'] = pole_prop.round(4)

print('\nClassification counts:')
print(df['Status'].value_counts())

# --------------------------------------------------
# Load mutation counts
# --------------------------------------------------
counts = pd.read_csv(COUNTS_PATH, sep='\t')

# Update cohort filter if needed
# Replace 'NUNES' with your actual cohort label in mutation_counts.tsv
nunes_counts = counts[counts['Cohort'] == 'NUNES'].copy()

print(f"\nMutation count rows for NUNES cohort: {len(nunes_counts)}")

# --------------------------------------------------
# Match sample IDs
# --------------------------------------------------
df_reset = df.reset_index()

# Rename first column to Samples if necessary
if df_reset.columns[0] != 'Samples':
    df_reset = df_reset.rename(columns={df_reset.columns[0]: 'Samples'})

# Try direct matching first
merged = df_reset.merge(
    nunes_counts[['sample_id', 'nMutations']],
    left_on='Samples',
    right_on='sample_id',
    how='left'
)

# Count successfully matched samples
matched = merged['nMutations'].notna().sum()

print(f"\\nDirect sample matches: {matched} / {len(merged)}")

print(borderline_pole.to_string())
