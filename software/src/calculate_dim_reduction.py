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
import scanpy as sc
import argparse
import os
import numpy as np
import time

warnings.simplefilter(action='ignore', category=FutureWarning)


def log_message(message, status="INFO"):
    """Logs messages in a structured format."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{status}] {message}")


def load_and_process_data(file_path):
    log_message("Starting data loading and preprocessing", "STEP")
    # Load the data from the CSV file
    raw_data_long = pd.read_csv(file_path)
    log_message(f"Loaded data from {file_path}, shape: {raw_data_long.shape}")

    # Validate and normalize column headers to support both legacy and new names
    base_required = {"Sample", "Ensembl Id", "Raw gene expression"}
    cell_headers = {"Cell Barcode", "Cell ID"}

    missing_base = base_required - set(raw_data_long.columns)
    has_cell_header = any(h in raw_data_long.columns for h in cell_headers)
    if missing_base or not has_cell_header:
        expected_desc = f"{sorted(base_required)} and one of {sorted(cell_headers)}"
        raise KeyError(f"Counts CSV must contain columns: {expected_desc}. Found: {list(raw_data_long.columns)}")

    # Normalize to legacy internal name 'Cell Barcode'
    if "Cell ID" in raw_data_long.columns and "Cell Barcode" not in raw_data_long.columns:
        raw_data_long = raw_data_long.rename(columns={"Cell ID": "Cell Barcode"})

    # Create a unique identifier for each cell
    raw_data_long['UniqueCellId'] = raw_data_long['Sample'] + '_' + raw_data_long['Cell Barcode']

    # Pivot the data to have genes as columns and UniqueCellId as rows
    raw_data = raw_data_long.pivot_table(
        index='UniqueCellId', 
        columns='Ensembl Id', 
        values='Raw gene expression', 
        fill_value=0
    )

    # Create AnnData object
    adata = sc.AnnData(raw_data)

    # Add SampleId and CellId metadata
    adata.obs['Sample'] = [uid.split('_')[0] for uid in adata.obs_names]
    adata.obs['Cell Barcode'] = [uid.split('_')[1] for uid in adata.obs_names]

    # Preprocessing steps: normalization, log transformation, and scaling
    log_message("Starting data normalization and transformation", "STEP")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    log_message("Normalization and normalization completed", "DONE")

    return adata


def save_pca(adata, output_dir):
    """Saves PCA results in long format: UniqueCellId, SampleId, CellId, PC, value"""
    # Get PCA matrix and cell names
    log_message("Saving PCA results", "STEP")
    pca_matrix = adata.obsm['X_pca']
    cell_ids = adata.obs_names

    # Convert to DataFrame and rename index
    pca_df = pd.DataFrame(pca_matrix, index=cell_ids)
    pca_df.index.name = 'UniqueCellId'  # Ensure the index has a name

    # Reset index before melting
    pca_df = pca_df.reset_index()

    # Melt into long format
    pca_df = pca_df.melt(id_vars=['UniqueCellId'], var_name='PC', value_name='value')

    # Extract SampleId and CellId from UniqueCellId
    pca_df['SampleId'] = pca_df['UniqueCellId'].apply(lambda x: x.split('_')[0])
    pca_df['CellId'] = pca_df['UniqueCellId'].apply(lambda x: x.split('_')[1])

    # Convert PC column to proper naming (e.g., "PC1", "PC2", etc.)
    pca_df['PC'] = pca_df['PC'].astype(int) + 1  # Convert to integer and shift index
    pca_df['PC'] = pca_df['PC'].apply(lambda x: f'PC{x}')

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
    embedding_df['SampleId'] = embedding_df['UniqueCellId'].apply(lambda x: x.split('_')[0])
    embedding_df['CellId'] = embedding_df['UniqueCellId'].apply(lambda x: x.split('_')[1])

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
    sc.tl.pca(adata, n_comps=n_pcs)

    # Save PCA results in the required long format
    save_pca(adata, output_dir)

    # Adaptive t-SNE parameters
    log_message("Running t-SNE", "STEP")
    n_cells = adata.n_obs
    perplexity = min(30, max(5, (n_cells - 1) // 3))  # Adjust perplexity
    learning_rate = max(50, min(200, n_cells / 12))   # Adaptive learning rate

    # Run t-SNE with adaptive parameters
    sc.tl.tsne(adata, n_pcs=n_pcs, perplexity=perplexity, early_exaggeration=12.0, learning_rate=learning_rate)
    tsne_results = adata.obsm['X_tsne']
    save_formatted_output(adata, tsne_results, 'TSNE', output_dir)

    # Adaptive UMAP parameters
    log_message("Running UMAP", "STEP")
    if n_cells < 1000:
        umap_n_neighbors = max(5, n_neighbors // 2)  # More local structure
        min_dist = 0.01
        spread = 1.0
    else:
        umap_n_neighbors = min(50, n_neighbors * 2)  # Emphasize global structure
        min_dist = 0.3
        spread = 1.5

    # Run UMAP with adaptive parameters
    sc.pp.neighbors(adata, n_neighbors=umap_n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata, n_components=3, min_dist=min_dist, spread=spread)
    umap_results = adata.obsm['X_umap']
    save_formatted_output(adata, umap_results, 'UMAP', output_dir)
    log_message("Dimensionality reduction completed", "DONE")


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process counts in long format and perform dimensionality reduction.')
    parser.add_argument('--file_path', type=str, help='Path to the counts CSV file.')
    parser.add_argument('--output_dir', type=str, help='Directory to store the output CSV files.')
    parser.add_argument('--n_pcs', type=int, default=50, help='Number of principal components (default: 50).')
    parser.add_argument('--n_neighbors', type=int, default=15, help='Number of neighbors for UMAP (default: 15).')
    return parser.parse_args()


def main():
    log_message("Starting block execution", "INFO")
    args = parse_arguments()
    adata = load_and_process_data(args.file_path)
    run_dimensionality_reduction(adata, args.output_dir, args.n_pcs, args.n_neighbors)
    log_message("Block execution finished", "INFO")


if __name__ == '__main__':
    main()
