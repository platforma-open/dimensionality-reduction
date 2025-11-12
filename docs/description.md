# Overview

Takes high-dimensional gene expression data from scRNA-seq preprocessing blocks (e.g., Cell Ranger) and transforms it into a lower-dimensional space while preserving biological variation using three complementary methods. Principal Component Analysis (PCA) reduces the data to a configurable number of principal components capturing major sources of variation. t-distributed Stochastic Neighbor Embedding (t-SNE) creates a two-dimensional embedding optimized for local structure preservation, while Uniform Manifold Approximation and Projection (UMAP) generates a three-dimensional embedding balancing local and global structure. Both t-SNE and UMAP use adaptive parameters that adjust based on dataset size, and all methods operate on the PCA space for computational efficiency.

When metadata covariates are provided, the block optionally performs batch correction using Harmony, which integrates cells across batches while preserving biological variation. Harmony correction is applied to the PCA space, and both UMAP and t-SNE embeddings are recomputed using the corrected principal components.

The resulting dimension values are used by downstream blocks: PCA components for clustering and pseudotime inference, and UMAP/t-SNE embeddings for visualization. The block provides both standard and batch-corrected embeddings for comparison.

The block uses scanpy v1.10.1 for dimensionality reduction algorithms and preprocessing. When using this block in your research, cite the scanpy publication (Wolf et al. 2018) listed below.

The following publications describe the methodologies used:

> Wolf, F. A., Angerer, P., & Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. _Genome Biology_ **19**, 15 (2018). [https://doi.org/10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)

> McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. _Journal of Open Source Software_ **3**(29), 861. [https://doi.org/10.21105/joss.00861](https://doi.org/10.21105/joss.00861)

> van der Maaten, L., & Hinton, G. (2008). Visualizing data using t-SNE. _Journal of Machine Learning Research_ **9**, 2579-2605. [http://www.jmlr.org/papers/v9/vandermaaten08a.html](http://www.jmlr.org/papers/v9/vandermaaten08a.html)

> Korsunsky, I., Millard, N., Fan, J. et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. _Nature Methods_ **16**, 1289–1296 (2019). [https://doi.org/10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0)
