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

