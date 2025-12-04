import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# --- Configuration ---
OUTDIR = os.path.join(os.getcwd(), "DE_Results_N2")

# Define filtering constants
UP_DIRECTION = 'Up-regulated (Tumoral)'
DOWN_DIRECTION = 'Down-regulated (Tumoral)'
N_TOP_UP = 15 # Extract the Top 15 over-expressed genes

# 1. Load and Prepare Data
results_path = "C:\\Users\\ricca\\OneDrive\\Desktop\\DE_Results_N2\\fold_change_results_N2.csv"
try:
    df = pd.read_csv(results_path)
    
    # Cleaning and indexing
    if 'GeneID' in df.columns:
        df = df.set_index('GeneID')
    if 'Unnamed: 0' in df.columns:
        df = df.drop(columns=['Unnamed: 0'])
    # Remove infinite or NaN log2fc values
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['log2fc', 'direction'])
    
except FileNotFoundError:
    print(f"ERROR: File not found at path: {results_path}")
    exit()

# Plot 1: Bar Plot of Differentially Expressed Gene Counts (Up vs Down) 
# Count genes by direction and filter only DEGs
deg_counts = df['direction'].value_counts()
deg_counts = deg_counts.filter(items=[UP_DIRECTION, DOWN_DIRECTION])

plt.figure(figsize=(7, 6))
# Bar plot of counts (Red for Up, Blue for Down)
sns.barplot(x=deg_counts.index, y=deg_counts.values, palette=['red', 'blue'])

plt.title('Count of Differentially Expressed Genes (Tumoral vs Normal)', fontsize=14)
plt.ylabel('Number of Genes', fontsize=12)
plt.xlabel('Direction of Change', fontsize=12)
plt.xticks(ticks=[0, 1], labels=['Up-regulated', 'Down-regulated'], rotation=0)

# Add counts on top of the bars
for i, count in enumerate(deg_counts.values):
    plt.text(i, count, str(count), ha='center', va='bottom', fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "bar_plot_deg_counts.png"), dpi=300)
plt.close() 

# Plot 2: Horizontal Bar Plot of Top Over-expressed Genes (log2FC magnitude)
# Filter Up-regulated genes and sort them by log2FC
up_regulated_df = df[df['direction'] == UP_DIRECTION].sort_values(by='log2fc', ascending=False)
top_up_df = up_regulated_df.head(N_TOP_UP)

plt.figure(figsize=(10, 8))
sns.barplot(
    x='log2fc',
    y=top_up_df.index.astype(str),
    data=top_up_df,
    color='red'
)

plt.title(f'Top {N_TOP_UP} Over-expressed Genes (Up-regulated) in Tumoral', fontsize=14)
plt.xlabel('log₂(Fold Change)', fontsize=12) 
plt.ylabel('Gene ID (Entrez ID)', fontsize=12)
plt.gca().invert_yaxis() # To have the highest value at the top
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, f"barplot_top_{N_TOP_UP}_up_genes.png"), dpi=300)
plt.close()