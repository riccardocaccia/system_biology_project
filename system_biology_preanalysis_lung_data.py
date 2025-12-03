import pandas as pd
import numpy as np
import os

# --- Configuration ---
OUTDIR = os.path.join(os.getcwd(), "DE_Results_N2")
os.makedirs(OUTDIR, exist_ok=True)
# ---------------------

exp_data = pd.read_csv("/Users/giuse/Downloads/GSE133206_raw_counts_GRCh38.p13_NCBI.tsv", sep="\t")
## 1. Load Data (Use your actual loading code here)

exp_data = exp_data.astype(float)

print(f"Data shape: {exp_data.shape}, Columns: {list(exp_data.columns)}")

# Define the column names for clarity
HEALTHY_COL = exp_data.columns[0]  # e.g., 'Healthy_S1'
TUMORAL_COL = exp_data.columns[1]  # e.g., 'Tumoral_S1'

# Check for zero values, which break log2 calculations
if (exp_data[[HEALTHY_COL, TUMORAL_COL]] == 0).any(axis=None):
    print("WARNING: Zero values found. Adding a pseudo-count (1.0) to avoid log(0).")
    # Add a small pseudo-count (e.g., 1.0 or the minimum non-zero value)
    exp_data += 1.0

## 2. Calculate Fold Change and log2FC
print("\nCalculating log2(Fold Change) directly...")

# Create the results DataFrame
res_df = exp_data.copy()

# A. Calculate Fold Change (FC)
# FC = Expression(Tumoral) / Expression(Healthy)
res_df['fold_change'] = res_df[TUMORAL_COL] / res_df[HEALTHY_COL]

# B. Calculate log2(FC)
# log2(FC) = log2(Tumoral / Healthy)
res_df['log2fc'] = np.log2(res_df['fold_change'])

# C. Classify the result based on log2fc
# We use arbitrary thresholds common in single-cell or low-replicate analysis:
LOGFC_THRESHOLD = 1.0 # Fold change of 2 (2^1 = 2)

res_df['direction'] = 'No Change'
res_df.loc[res_df['log2fc'] > LOGFC_THRESHOLD, 'direction'] = 'Up-regulated (Tumoral)'
res_df.loc[res_df['log2fc'] < -LOGFC_THRESHOLD, 'direction'] = 'Down-regulated (Tumoral)'

# Prepare final report, including the raw expression values
final_res_df = res_df[[HEALTHY_COL, TUMORAL_COL, 'fold_change', 'log2fc', 'direction']].copy()

# Sort by the magnitude of the log2fc (largest difference first)
final_res_df['abs_log2fc'] = np.abs(final_res_df['log2fc'])
final_res_df = final_res_df.sort_values(by='abs_log2fc', ascending=False)
final_res_df = final_res_df.drop(columns=['abs_log2fc'])

## 3. Save Results
output_path = os.path.join(OUTDIR, "fold_change_results_N2.csv")
final_res_df.to_csv(output_path)

print("\n--- Summary ---")
print(f"Total genes analyzed: {len(final_res_df)}")
print(f"Genes with |log2fc| > {LOGFC_THRESHOLD}: {len(final_res_df[final_res_df['direction'] != 'No Change'])}")
print(f"Results saved to: {output_path}")

print("\nTop 5 Differentially Expressed Genes (by magnitude of log2FC):")
print(final_res_df[['log2fc', 'fold_change', 'direction']].head())

