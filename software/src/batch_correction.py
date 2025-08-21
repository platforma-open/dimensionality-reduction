import pandas as pd
import scanpy as sc
import harmonypy as hm
import argparse
import os
import numpy as np
from scipy.sparse import csr_matrix

# Argument parsing
parser = argparse.ArgumentParser(description="Batch correction for scRNA-seq using ComBat (counts) and Harmony (embeddings)")
parser.add_argument("--counts", help="Path to raw counts CSV file", required=True)
parser.add_argument("--metadata", help="Path to metadata CSV file", required=True)
parser.add_argument("--output", help="Path to output directory", required=True)
args = parser.parse_args()

# Create output directory if it doesn't exist
os.makedirs(args.output, exist_ok=True)

# Load raw count data (long format)
print("Loading raw counts...")
counts_df = pd.read_csv(args.counts, dtype={"Sample": str, "Cell Barcode": str, "Ensembl Id": str})

# Load metadata
print("Loading metadata...")
metadata_df = pd.read_csv(args.metadata, dtype=str)

# Ensure metadata contains at least two columns
if metadata_df.shape[1] < 2:
    raise ValueError("Metadata must have at least two columns: 'Sample' and at least one metadata column.")

# Detect metadata columns
metadata_columns = metadata_df.columns[1:].tolist()  # Exclude "Sample"
combat_column = metadata_columns[0]  # First column for ComBat
harmony_columns = metadata_columns  # Use all columns for Harmony

print(f"✅ Detected metadata columns: {metadata_columns}")
print(f"✅ Using '{combat_column}' for ComBat correction")
print(f"✅ Using {harmony_columns} for Harmony correction")

# Merge metadata into count matrix
counts_df = counts_df.merge(metadata_df, on="Sample", how="left")

# Convert long format to cell x gene matrix
print("Converting long format to cell-by-gene matrix...")
cell_gene_matrix = counts_df.pivot_table(
    index=["Sample", "Cell Barcode"], columns="Ensembl Id", values="Raw gene expression", aggfunc="sum", fill_value=0
)

# Convert to AnnData
print("Creating AnnData object...")
adata = sc.AnnData(X=csr_matrix(cell_gene_matrix.values))  # Convert to sparse matrix
adata.obs = cell_gene_matrix.index.to_frame(index=False)  # Preserve Sample & Cell Barcode
adata.var_names = cell_gene_matrix.columns  # Gene names

# Merge metadata into AnnData object
adata.obs = adata.obs.merge(metadata_df, on="Sample", how="left")

# Ensure metadata columns are categorical
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

# Ensure UniqueCellId is set before running Harmony
adata_harmony.obs["UniqueCellId"] = adata_harmony.obs["Sample"] + "_" + adata_harmony.obs["Cell Barcode"]
adata_harmony.obs.set_index("UniqueCellId", inplace=True)

# Run PCA on raw counts (no ComBat)
print("Running PCA on raw counts...")
sc.pp.normalize_total(adata_harmony, target_sum=1e4)
sc.pp.log1p(adata_harmony)
sc.pp.scale(adata_harmony)
sc.tl.pca(adata_harmony, n_comps=100)

# Apply Harmony with multiple metadata columns
print(f"Applying Harmony batch correction using {harmony_columns}...")
harmony_results = hm.run_harmony(np.array(adata_harmony.obsm["X_pca"]).T, adata_harmony.obs, harmony_columns, max_iter_harmony=50)
adata_harmony.obsm["X_pca_harmony"] = harmony_results.Z_corr.T  # Transpose back

# Save Harmony-corrected PCA embeddings
print("Saving Harmony PCA embeddings...")
harmony_df = pd.DataFrame(adata_harmony.obsm["X_pca_harmony"], index=adata_harmony.obs.index)
harmony_df.columns = [f"PC{i+1}" for i in range(harmony_df.shape[1])]
harmony_df = harmony_df.reset_index().melt(id_vars=["UniqueCellId"], var_name="PC", value_name="value")

harmony_df["Sample"] = harmony_df["UniqueCellId"].apply(lambda x: x.split('_')[0])
harmony_df["Cell Barcode"] = harmony_df["UniqueCellId"].apply(lambda x: x.split('_')[1])
harmony_df = harmony_df[["Sample", "Cell Barcode", "PC", "value"]]
harmony_df.to_csv(os.path.join(args.output, "harmony_results.csv"), index=False)

# Run UMAP and tSNE
print("Running UMAP...")
sc.pp.neighbors(adata_harmony, use_rep="X_pca_harmony")
sc.tl.umap(adata_harmony)

print("Running tSNE...")
sc.tl.tsne(adata_harmony, use_rep="X_pca_harmony")

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