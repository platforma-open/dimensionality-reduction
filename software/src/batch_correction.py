import pandas as pd
import polars as pl
import scanpy as sc
import harmonypy as hm
import argparse
import os
import numpy as np
import time
from scipy.sparse import csr_matrix

np.random.seed(0)

def log_message(message, status="INFO"):
    """Logs messages in a structured format."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{status}] {message}")

# Argument parsing
parser = argparse.ArgumentParser(description="Batch correction for scRNA-seq using Harmony (embeddings)")
parser.add_argument("--counts", help="Path to raw counts CSV file", required=True)
parser.add_argument("--metadata", help="Path to metadata CSV file", required=True)
parser.add_argument("--output", help="Path to output directory", required=True)
args = parser.parse_args()

# Create output directory if it doesn't exist
os.makedirs(args.output, exist_ok=True)

# 1. LOAD DATA WITH POLARS & CATEGORICAL OPTIMIZATION
log_message("Loading raw counts with Polars and Categorical optimization", "STEP")

# Peek schema to handle flexible headers and identify columns for Categorical casting
temp_scan = pl.scan_csv(args.counts)
file_schema = temp_scan.collect_schema()
column_names = set(file_schema.keys())

## Then, make sure we have the required columns
base_required = {"Sample", "Ensembl Id", "Raw gene expression"}
cell_headers = {"Cell Barcode", "Cell ID"}

missing_base = base_required - set(column_names)
has_cell_header = any(h in column_names for h in cell_headers)
if missing_base or not has_cell_header:
    expected_desc = f"{sorted(base_required)} and one of {sorted(cell_headers)}"
    raise KeyError(f"Counts CSV must contain columns: {expected_desc}. Found: {list(column_names)}")

## Second, build schema_overrides only for columns that actually exist
schema_overrides = {}
categorical_candidates = ["Sample", "Ensembl Id", "Cell Barcode", "Cell ID"]
for col in categorical_candidates:
    if col in column_names:
        schema_overrides[col] = pl.Categorical

counts_pl = pl.read_csv(args.counts, schema_overrides=schema_overrides)

# Normalize to legacy internal name 'Cell Barcode'
if "Cell ID" in counts_pl.columns and "Cell Barcode" not in counts_pl.columns:
    counts_pl = counts_pl.rename({"Cell ID": "Cell Barcode"})

# Create a unique identifier for each cell
SEPARATOR = '|||'
counts_pl = counts_pl.with_columns(
    (pl.col('Sample').cast(str) + pl.lit(SEPARATOR) + pl.col('Cell Barcode').cast(str))
    .cast(pl.Categorical)
    .alias('UniqueCellId')
)

log_message("Creating sparse matrix from long format data", "STEP")

# Extract integer codes directly from categorical columns (much faster than join)
row_codes_raw = counts_pl['UniqueCellId'].to_physical().to_numpy()
col_codes_raw = counts_pl['Ensembl Id'].to_physical().to_numpy()
# Use float32 for expression values (Scanpy standard)
expression_values = counts_pl['Raw gene expression'].cast(pl.Float32).to_numpy()

# Remap codes to 0-indexed contiguous using np.unique (efficient integer-based mapping)
u_row_phys, row_idx = np.unique(row_codes_raw, return_inverse=True)
u_col_phys, col_idx = np.unique(col_codes_raw, return_inverse=True)

# Map labels to sorted ranks and get sorted unique IDs (efficiently processing only unique labels)
unique_cell_ids, row_map = np.unique(counts_pl['UniqueCellId'].cat.get_categories().gather(u_row_phys).to_numpy(), return_inverse=True)
unique_gene_ids, col_map = np.unique(counts_pl['Ensembl Id'].cat.get_categories().gather(u_col_phys).to_numpy(), return_inverse=True)

# Final row and column codes are the mapped indices
row_codes = row_map[row_idx].astype(np.int32)
col_codes = col_map[col_idx].astype(np.int32)

# Delete Polars objects to free memory - all needed data has been extracted
del counts_pl, row_codes_raw, col_codes_raw, u_row_phys, u_col_phys, row_idx, col_idx, row_map, col_map

n_cells = len(unique_cell_ids)
n_genes = len(unique_gene_ids)

# Vectorized extraction of Sample and Cell Barcode for obs (more efficient than doing it later)
obs_df = pd.DataFrame(index=unique_cell_ids)
split_ids = pd.Series(unique_cell_ids).str.split(SEPARATOR, n=1, expand=True, regex=False)
obs_df['Sample'] = split_ids[0].values
obs_df['Cell Barcode'] = split_ids[1].values

# Create the sparse matrix
adata = sc.AnnData(
    X=csr_matrix((expression_values, (row_codes, col_codes)), shape=(n_cells, n_genes), dtype=np.float32),
    obs=obs_df,
    var=pd.DataFrame(index=unique_gene_ids)
)
log_message(f"AnnData object created: {n_cells} cells x {n_genes} genes", "DONE")

# 3. LOAD & MERGE METADATA
log_message("Loading and merging metadata", "STEP")
metadata_df = pd.read_csv(args.metadata, dtype=str)

# Ensure metadata contains at least two columns
if metadata_df.shape[1] < 2:
    raise ValueError("Metadata must have at least two columns: 'Sample' and at least one metadata column.")

metadata_columns = metadata_df.columns[1:].tolist() # Exclude "Sample"
combat_column = metadata_columns[0]
harmony_columns = metadata_columns

log_message(f"Detected metadata columns: {metadata_columns}", "INFO")
log_message(f"Using '{combat_column}' for ComBat correction", "INFO")
log_message(f"Using {harmony_columns} for Harmony correction", "INFO")

# Merge metadata into adata.obs
# Sample column already exists in obs, so merge will align correctly
adata.obs = adata.obs.merge(metadata_df, on="Sample", how="left")
adata.obs.index = unique_cell_ids # Restore the UniqueCellId index lost during merge

# Ensure metadata columns are categorical for Harmony
for col in metadata_columns:
    adata.obs[col] = adata.obs[col].fillna("Unknown").astype("category")

# # **Branch 1: Apply ComBat for batch correction on raw counts**
# adata_combat = adata.copy()  # Separate object for ComBat correction

# # Remove genes with zero variance within any batch
# print("Filtering genes with zero variance within any batch before ComBat...")
# dense_matrix = adata_combat.X.toarray() if isinstance(adata_combat.X, csr_matrix) else adata_combat.X
# zero_variance_per_batch = [
#     (dense_matrix[adata_combat.obs[combat_column] == batch, :].var(axis=0) == 0)
#     for batch in adata_combat.obs[combat_column].unique()
# ]
# zero_variance_in_any_batch = pd.DataFrame(zero_variance_per_batch).any(axis=0)
# adata_combat = adata_combat[:, ~zero_variance_in_any_batch]
# print(f"🔍 Removed {zero_variance_in_any_batch.sum()} genes with zero variance in at least one batch.")

# # Apply ComBat
# print(f"Applying ComBat batch correction using '{combat_column}'...")
# sc.pp.combat(adata_combat, key=combat_column)

# # Prevent NaN/negative values before normalization
# adata_combat.X = np.nan_to_num(adata_combat.X, nan=0.0, posinf=0.0, neginf=0.0)
# adata_combat.X[adata_combat.X < 0] = 0

# # **Fix: Save batch-corrected counts BEFORE normalization**
# print("Saving batch-corrected count matrix...")
# corrected_counts_df = pd.DataFrame(adata_combat.X, index=adata_combat.obs.index, columns=adata_combat.var_names)
# corrected_counts_df["Sample"] = adata_combat.obs["Sample"].values
# corrected_counts_df["Cell Barcode"] = adata_combat.obs["Cell Barcode"].values
# corrected_counts_df = corrected_counts_df.melt(id_vars=["Sample", "Cell Barcode"], var_name="Ensembl Id", value_name="Raw gene expression")
# corrected_counts_df.to_csv(os.path.join(args.output, "batch_corrected_counts.csv"), index=False)

# # **Fix: Create a copy before normalization**
# adata_norm = adata_combat.copy()

# # Normalize batch-corrected counts
# print("Normalizing batch-corrected counts...")
# sc.pp.normalize_total(adata_norm, target_sum=1e4)
# sc.pp.log1p(adata_norm)

# # Save batch-corrected normalized counts
# print("Saving batch-corrected normalized count matrix...")
# normalized_counts_df = pd.DataFrame(adata_norm.X, index=adata_norm.obs.index, columns=adata_norm.var_names)
# normalized_counts_df["Sample"] = adata_norm.obs["Sample"].values
# normalized_counts_df["Cell Barcode"] = adata_norm.obs["Cell Barcode"].values
# normalized_counts_df = normalized_counts_df.melt(id_vars=["Sample", "Cell Barcode"], var_name="Ensembl Id", value_name="Normalized gene expression")
# normalized_counts_df.to_csv(os.path.join(args.output, "batch_corrected_normalized_counts.csv"), index=False)

# **Branch 2: Apply Harmony for embedding correction**
adata_harmony = adata.copy()  # Separate object for Harmony correction

# Run PCA on raw counts (no ComBat)
print("Running PCA on raw counts...")
sc.pp.normalize_total(adata_harmony, target_sum=1e4)
sc.pp.log1p(adata_harmony)
sc.pp.scale(adata_harmony)
sc.tl.pca(adata_harmony, n_comps=100, random_state=0)

# Apply Harmony with multiple metadata columns
print(f"Applying Harmony batch correction using {harmony_columns}...")
# Harmony expects cells as columns, so we transpose the PCA matrix (n_cells x n_pcs -> n_pcs x n_cells)
harmony_results = hm.run_harmony(np.array(adata_harmony.obsm["X_pca"]).T, adata_harmony.obs, harmony_columns, max_iter_harmony=50)
adata_harmony.obsm["X_pca_harmony"] = harmony_results.Z_corr.T  # Transpose back to n_cells x n_pcs

# Save Harmony-corrected PCA embeddings
print("Saving Harmony PCA embeddings...")
harmony_df = pd.DataFrame(adata_harmony.obsm["X_pca_harmony"], index=adata_harmony.obs.index)
harmony_df.columns = [f"PC{i+1}" for i in range(harmony_df.shape[1])]
harmony_df = harmony_df.reset_index().rename(columns={'index': 'UniqueCellId'})
harmony_df = harmony_df.melt(id_vars=["UniqueCellId"], var_name="PC", value_name="value")

# Use SEPARATOR and regex=False for metadata integrity
split_harmony = harmony_df["UniqueCellId"].str.split(SEPARATOR, n=1, expand=True, regex=False)
harmony_df["Sample"] = split_harmony[0]
harmony_df["Cell Barcode"] = split_harmony[1]
harmony_df = harmony_df[["Sample", "Cell Barcode", "PC", "value"]]
harmony_df.to_csv(os.path.join(args.output, "harmony_results.csv"), index=False)

# Run UMAP and tSNE
print("Running UMAP...")
sc.pp.neighbors(adata_harmony, use_rep="X_pca_harmony", random_state=0)
sc.tl.umap(adata_harmony, random_state=0)

print("Running tSNE...")
sc.tl.tsne(adata_harmony, use_rep="X_pca_harmony", random_state=0)

# **Save Outputs**
print("Saving UMAP coordinates...")
umap_df = pd.DataFrame(adata_harmony.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata_harmony.obs.index)
umap_df["Sample"] = adata_harmony.obs["Sample"].values
umap_df["Cell Barcode"] = adata_harmony.obs["Cell Barcode"].values
umap_df = umap_df[["Sample", "Cell Barcode", "UMAP1", "UMAP2"]]  # Ensure correct order
umap_df.to_csv(os.path.join(args.output, "umap_dimensions.csv"), index=False)

print("Saving tSNE coordinates...")
tsne_df = pd.DataFrame(adata_harmony.obsm["X_tsne"], columns=["tSNE1", "tSNE2"], index=adata_harmony.obs.index)
tsne_df["Sample"] = adata_harmony.obs["Sample"].values
tsne_df["Cell Barcode"] = adata_harmony.obs["Cell Barcode"].values
tsne_df = tsne_df[["Sample", "Cell Barcode", "tSNE1", "tSNE2"]]  # Ensure correct order
tsne_df.to_csv(os.path.join(args.output, "tsne_dimensions.csv"), index=False)
print("✅ All outputs saved successfully!")
