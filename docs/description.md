# Overview

Takes high-dimensional gene expression data and transforms it into a lower-dimensional space while preserving the most important biological variation. The block takes the output from any scRNA-seq preprocessing block (e.g., Cell Ranger) as input. It generates plots with tSNE and UMAP projections to aid dataset exploration and outputs dimension values to be used by downstream blocks (e.g. Cell Browser).
