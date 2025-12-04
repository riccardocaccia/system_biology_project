import pandas as pd
import os
import numpy as np

# --- Configuration ---
OUTDIR = os.path.join(os.getcwd(), "DE_Results_N2")
N_TOP_UP = 15 
UP_DIRECTION = 'Up-regulated (Tumoral)'

# 1. Load data
results_path = "C:\\Users\\ricca\\OneDrive\\Desktop\\DE_Results_N2\\fold_change_results_N2.csv"
df = pd.read_csv(results_path)

# 2. Filter and Select Top Up-regulated Genes
up_regulated_df = df[df['direction'] == UP_DIRECTION].sort_values(by='log2fc', ascending=False)
top_up_df = up_regulated_df.head(N_TOP_UP)

# 3. Extract the list of Gene IDs (Entrez IDs)
gene_id_list = top_up_df.index.astype(int).tolist()

print("\n--- LISTA DI ENTRÉZ ID READY FOR ORA (Top 15 Up-regulated) ---")
# Stampa gli ID separati da spazi/nuova riga per il copia-incolla facile
print("\n".join(map(str, gene_id_list)))
