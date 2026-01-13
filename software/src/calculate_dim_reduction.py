import warnings
import sys

# Remove all existing warning filters
warnings.resetwarnings()

# Force Python to ignore FutureWarnings globally
warnings.simplefilter("ignore", category=FutureWarning)

# Ensure all warnings are ignored, even if a library resets them
if not sys.warnoptions:
    import os
    os.environ["PYTHONWARNINGS"] = "ignore::FutureWarning"


import pandas as pd
import polars as pl
import scanpy as sc
import argparse
import os
import numpy as np
import time
from scipy.sparse import csr_matrix
import pyarrow

warnings.simplefilter(action='ignore', category=FutureWarning)


def log_message(message, status="INFO"):
    """Logs messages in a structured format."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{status}] {message}")


def load_and_process_data(file_path, hvg_count=0):
    log_message("Starting data loading and preprocessing with Polars", "STEP")
    
    # MEMORY OPTIMIZATION: Read string columns with repeated values as categoricals
    # Sample,Cell Barcode/Cell ID, Ensembl Id
    
    ## First, peek at the schema to see which columns exist (lazy operation, doesn't load data)
    log_message("Scanning CSV schema to identify columns", "STEP")
    temp_scan = pl.scan_csv(file_path)
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
    
    log_message(f"Reading CSV with categorical types for: {list(schema_overrides.keys())}", "STEP")
    
    ## Finally, load the data with categorical types for repeated string columns
    raw_data_long_pl = pl.read_csv(file_path, schema_overrides=schema_overrides)
    log_message(f"Loaded data from {file_path}, shape: {raw_data_long_pl.shape}")

    # Normalize to legacy internal name 'Cell Barcode'
    if "Cell ID" in raw_data_long_pl.columns and "Cell Barcode" not in raw_data_long_pl.columns:
        raw_data_long_pl = raw_data_long_pl.rename({"Cell ID": "Cell Barcode"})

    # Create a unique identifier for each cell
    SEPARATOR = '|||'
    raw_data_long_pl = raw_data_long_pl.with_columns(
        (pl.col('Sample').cast(str) + pl.lit(SEPARATOR) + pl.col('Cell Barcode').cast(str))
        .cast(pl.Categorical)
        .alias('UniqueCellId')
    )
    
    log_message("Creating sparse matrix from long format data", "STEP")
    
    # Extract integer codes directly from categorical columns (much faster than join)
    row_codes_raw = raw_data_long_pl['UniqueCellId'].to_physical().to_numpy()
    col_codes_raw = raw_data_long_pl['Ensembl Id'].to_physical().to_numpy()
    expression_values = raw_data_long_pl['Raw gene expression'].to_numpy()
    
    # Remap codes to 0-indexed contiguous using np.unique (efficient integer-based mapping)
    u_row_phys, row_idx = np.unique(row_codes_raw, return_inverse=True)
    u_col_phys, col_idx = np.unique(col_codes_raw, return_inverse=True)

    # Map labels to sorted ranks and get sorted unique IDs (efficiently processing only unique labels)
    unique_cell_ids, row_map = np.unique(raw_data_long_pl['UniqueCellId'].cat.get_categories().gather(u_row_phys).to_numpy(), return_inverse=True)
    unique_gene_ids, col_map = np.unique(raw_data_long_pl['Ensembl Id'].cat.get_categories().gather(u_col_phys).to_numpy(), return_inverse=True)

    # Final row and column codes are the mapped indices
    row_codes = row_map[row_idx].astype(np.int32)
    col_codes = col_map[col_idx].astype(np.int32)
    
    # Delete intermediate objects to free memory
    del raw_data_long_pl, row_codes_raw, col_codes_raw, u_row_phys, u_col_phys, row_idx, col_idx, row_map, col_map
    
    n_unique_cells = len(unique_cell_ids)
    n_unique_genes = len(unique_gene_ids)
    
    # Create the sparse matrix
    sparse_matrix = csr_matrix(
        (expression_values, (row_codes, col_codes)),
        shape=(n_unique_cells, n_unique_genes),
        dtype=np.float32
    )

    # Create AnnData object, which requires pandas DataFrames for obs and var
    # Use to_pandas() instead of to_list() - more efficient for categoricals
    adata = sc.AnnData(
        sparse_matrix,
        obs=pd.DataFrame(index=unique_cell_ids),
        var=pd.DataFrame(index=unique_gene_ids)
    )
    log_message("Sparse matrix and AnnData object created", "DONE")

    # Add SampleId and CellId metadata using vectorized operations (much faster than list comprehensions)
    obs_names_series = pd.Series(adata.obs_names)
    adata.obs['Sample'] = obs_names_series.str.split(SEPARATOR, n=1, expand=True)[0].values
    adata.obs['Cell Barcode'] = obs_names_series.str.split(SEPARATOR, n=1, expand=True)[1].values

    # Preprocessing steps: normalization, log transformation, finding HVGs, and scaling
    log_message("Starting data normalization and transformation", "STEP")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    if hvg_count > 0:
        log_message(f"Finding {hvg_count} highly variable genes", "STEP")
        sc.pp.highly_variable_genes(adata, n_top_genes=hvg_count)
        log_message("Highly variable genes found", "DONE")

        # Subset the data to the most variable genes
        adata.raw = adata  # Store the full dataset
        adata = adata[:, adata.var.highly_variable].copy()
    else:
        log_message("Skipping highly variable genes selection (using all genes)", "STEP")

    log_message("Scaling data to unit variance and zero mean", "STEP")
    sc.pp.scale(adata, max_value=10)
    log_message("Scaling completed", "DONE")

    log_message("Normalization and transformation completed", "DONE")
    
    return adata


def save_pca(adata, output_dir):
    """Saves PCA results in long format: UniqueCellId, SampleId, CellId, PC, value"""
    # Get PCA matrix and cell names
    log_message("Saving PCA results", "STEP")
    pca_matrix = adata.obsm['X_pca']
    cell_ids = adata.obs_names
    n_cells, n_pcs = pca_matrix.shape

    # MEMORY OPTIMIZATION: Build long format directly instead of using melt()
    # This avoids creating a wide DataFrame which can be memory-intensive for many PCs
    log_message(f"Building long format for {n_cells} cells and {n_pcs} PCs", "STEP")
    
    # Pre-allocate arrays for better memory efficiency
    unique_cell_ids = np.repeat(cell_ids, n_pcs)
    pc_indices = np.tile(np.arange(n_pcs), n_cells)
    pc_values = pca_matrix.flatten()
    
    # Create DataFrame directly in long format
    pca_df = pd.DataFrame({
        'UniqueCellId': unique_cell_ids,
        'PC': pc_indices,
        'value': pc_values
    })

    # Extract SampleId and CellId from UniqueCellId
    SEPARATOR = '|||'
    # Ensure UniqueCellId is string type and split explicitly
    # This step is important since we load some string columns as categoricals
    unique_cell_id_series = pca_df['UniqueCellId'].astype(str)
    split_result = unique_cell_id_series.str.split(pat=SEPARATOR, n=1, expand=True, regex=False)
    pca_df['SampleId'] = split_result[0]
    pca_df['CellId'] = split_result[1]

    # Convert PC column to proper naming (e.g., "PC1", "PC2", etc.) using vectorized operations
    pca_df['PC'] = 'PC' + (pca_df['PC'].astype(int) + 1).astype(str)

    # Adjust column order to match other outputs
    pca_df = pca_df[['UniqueCellId', 'SampleId', 'CellId', 'PC', 'value']]

    # Save to CSV
    pca_df.to_csv(os.path.join(output_dir, "pca_results.csv"), index=False)
    log_message("PCA results saved", "DONE")


def save_formatted_output(adata, embedding, embedding_name, output_dir):
    # Ensure embedding has exactly 3 dimensions
    log_message(f"Saving {embedding_name} results", "STEP")
    if embedding.shape[1] < 3:
        padding = np.zeros((embedding.shape[0], 3 - embedding.shape[1]))
        embedding = np.hstack((embedding, padding))
    elif embedding.shape[1] > 3:
        embedding = embedding[:, :3]

    # Create DataFrame with the correct format
    embedding_df = pd.DataFrame(
        embedding,
        index=adata.obs_names,
        columns=[f'{embedding_name}{i+1}' for i in range(3)]
    )
    embedding_df = embedding_df.reset_index().rename(columns={'index': 'UniqueCellId'})
    SEPARATOR = '|||'
    # Ensure UniqueCellId is string type and split explicitly
    # This step is important since we load some string columns as categoricals
    unique_cell_id_series = embedding_df['UniqueCellId'].astype(str)
    split_result = unique_cell_id_series.str.split(pat=SEPARATOR, n=1, expand=True, regex=False)
    embedding_df['SampleId'] = split_result[0]
    embedding_df['CellId'] = split_result[1]

    # Reorder columns
    embedding_df = embedding_df[['UniqueCellId', 'SampleId', 'CellId', f'{embedding_name}1', f'{embedding_name}2', f'{embedding_name}3']]

    # Save to CSV
    output_path = os.path.join(output_dir, f'{embedding_name.lower()}_results.csv')
    embedding_df.to_csv(output_path, index=False)
    log_message(f"{embedding_name} results saved", "DONE")


def run_dimensionality_reduction(adata, output_dir, n_pcs, n_neighbors):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Run PCA
    log_message("Running PCA", "STEP")
    sc.tl.pca(adata, n_comps=n_pcs, random_state=0)

    # Save PCA results in the required long format
    save_pca(adata, output_dir)

    # Adaptive t-SNE parameters
    log_message("Running t-SNE", "STEP")
    n_cells = adata.n_obs
    perplexity = min(30, max(5, (n_cells - 1) // 3))  # Adjust perplexity
    learning_rate = max(50, min(200, n_cells / 12))   # Adaptive learning rate

    # Run t-SNE with adaptive parameters
    sc.tl.tsne(adata, n_pcs=n_pcs, perplexity=perplexity, early_exaggeration=12.0, learning_rate=learning_rate, random_state=0)
    tsne_results = adata.obsm['X_tsne']
    save_formatted_output(adata, tsne_results, 'TSNE', output_dir)

    # Adaptive UMAP parameters
    log_message("Running UMAP", "STEP")
    if n_cells < 1000:
        umap_n_neighbors = max(5, n_neighbors // 4)  # More local structure
        min_dist = 0.01
        spread = 1.0
    else:
        umap_n_neighbors = min(50, n_neighbors)  # Emphasize global structure
        min_dist = 0.3
        spread = 1.5

    # Run UMAP with adaptive parameters
    sc.pp.neighbors(adata, n_neighbors=umap_n_neighbors, n_pcs=n_pcs, random_state=0)
    sc.tl.umap(adata, n_components=3, min_dist=min_dist, spread=spread, random_state=0)
    umap_results = adata.obsm['X_umap']
    save_formatted_output(adata, umap_results, 'UMAP', output_dir)
    log_message("Dimensionality reduction completed", "DONE")


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process counts in long format and perform dimensionality reduction.')
    parser.add_argument('--file_path', type=str, help='Path to the counts CSV file.')
    parser.add_argument('--output_dir', type=str, help='Directory to store the output CSV files.')
    parser.add_argument('--n_pcs', type=int, default=50, help='Number of principal components (default: 50).')
    parser.add_argument('--n_neighbors', type=int, default=15, help='Number of neighbors for UMAP (default: 15).')
    parser.add_argument('--hvg_count', type=int, default=0, help='Number of highly variable genes to use (0 to disable).')
    return parser.parse_args()


def main():
    log_message("Starting block execution", "INFO")
    np.random.seed(0)
    args = parse_arguments()
    adata = load_and_process_data(args.file_path, hvg_count=args.hvg_count)
    run_dimensionality_reduction(adata, args.output_dir, args.n_pcs, args.n_neighbors)
    log_message("Block execution finished", "INFO")


if __name__ == '__main__':
    main()
