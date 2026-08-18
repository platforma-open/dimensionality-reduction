---
'@platforma-open/milaboratories.dimensionality-reduction.model': patch
'@platforma-open/milaboratories.dimensionality-reduction.ui': patch
'@platforma-open/milaboratories.dimensionality-reduction.workflow': patch
'@platforma-open/milaboratories.dimensionality-reduction.software': patch
'@platforma-open/milaboratories.dimensionality-reduction': patch
---

Export the counts matrix as Parquet instead of CSV, and give that conversion its own budget.

The `xsv.exportFrame` step that materialises the long-format counts p-column for the Python tools
was OOM-killed (exit -1) on large datasets — the p-frame reached 240M+ rows and the step ran with a
hardcoded 16 GiB / 1 CPU. Counts are one row per non-zero (cell, gene) pair with three
high-cardinality string columns, which is the worst case for CSV text encoding; upstream blocks
already store the column as Parquet, so the export was inflating columnar data into text for no
benefit.

- workflow: `xsv.exportFrame([rawCounts], "parquet", ...)` with a dedicated 32 GiB / 4 CPU budget,
  separate from the (much smaller) UMAP/t-SNE/PCA result imports which keep 16 GiB / 1 CPU. The
  intermediate is renamed `csvCounts` -> `countsParquet` and `rawCounts.csv` -> `rawCounts.parquet`
  end-to-end through `dim-reduction-calculation` and `batch-correction`. The covariates export stays
  CSV.
- software: `calculate_dim_reduction.py` and `batch_correction.py` read the counts via
  `pl.scan_parquet`, casting the repeated string columns to `Categorical` inside the scan plan
  rather than through `read_csv(schema_overrides=...)`.

Migrate block onto the structurer (block-tools 2.12.13) — full SDK upgrade: model/ui-vue 1.81.1, workflow-tengo 6.8.2, tengo-builder 4.0.22, package-builder 3.14.2, test 1.81.2. Adopts the canonical tool-managed layout (oxlint/oxfmt, tsconfig, turbo, CI workflows, managed package.json + catalog) and the slim facade for the root block package. Author-code fixes for the SDK majors: explicit type argument on the `isPColumn` filter feeding `createPFrame`, removal of the retired `@platforma-sdk/ui-vue/styles` import, and the model export renamed `model` -> `platforma` for the facade.
