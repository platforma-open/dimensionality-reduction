# @platforma-open/milaboratories.dimensionality-reduction.ui

## 1.11.1

### Patch Changes

- 11a9a50: Export the counts matrix as Parquet instead of CSV, and give that conversion its own budget.

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

- Updated dependencies [11a9a50]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.10.1

## 1.11.0

### Minor Changes

- a04ff01: Migrate to latest layout and improve memory efficiency

### Patch Changes

- Updated dependencies [a04ff01]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.10.0

## 1.10.5

### Patch Changes

- 39d2f9b: technical release
- b7b0041: technical release
- 925cfba: technical release
- cdbee0e: technical release
- Updated dependencies [39d2f9b]
- Updated dependencies [b7b0041]
- Updated dependencies [925cfba]
- Updated dependencies [cdbee0e]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.9.5

## 1.10.4

### Patch Changes

- Updated dependencies [7e91dd9]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.9.4

## 1.10.3

### Patch Changes

- Updated dependencies [8c1e2f0]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.9.3

## 1.10.2

### Patch Changes

- Updated dependencies [0893c92]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.9.2

## 1.10.1

### Patch Changes

- Updated dependencies [1759fcb]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.9.1

## 1.10.0

### Minor Changes

- 078da45: Add batch correction and update to latest SDK

### Patch Changes

- Updated dependencies [078da45]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.9.0

## 1.9.2

### Patch Changes

- ed32e66: Fixed github workflow build
- Updated dependencies [ed32e66]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.8.2

## 1.9.1

### Patch Changes

- Updated dependencies [cd0b06b]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.8.1

## 1.9.0

### Minor Changes

- 205d015: Fixed script version. Added tabs for plot visualization. Minor fixes and reduced redundancy.

### Patch Changes

- Updated dependencies [205d015]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.8.0

## 1.8.0

### Minor Changes

- 13f4377: Migrate to PlTabs to switch between UMAP and tSNE

### Patch Changes

- Updated dependencies [13f4377]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.7.0

## 1.7.0

### Minor Changes

- 3128a46: Update SDK packages

## 1.6.0

### Minor Changes

- 7b57a6d: Updated UI

### Patch Changes

- Updated dependencies [7b57a6d]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.6.0

## 1.5.0

### Minor Changes

- 7aa1003: Improved visualization

### Patch Changes

- Updated dependencies [7aa1003]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.5.0

## 1.4.2

### Patch Changes

- 655312c: Updated dependencies
- Updated dependencies [655312c]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.4.2

## 1.4.1

### Patch Changes

- 2dfec16: Updated dependencies
- Updated dependencies [2dfec16]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.4.1

## 1.4.0

### Minor Changes

- a26208b: Updated UI

### Patch Changes

- Updated dependencies [a26208b]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.4.0

## 1.3.0

### Minor Changes

- 8aeec7f: Compatibility with Leiden Clustering

### Patch Changes

- Updated dependencies [8aeec7f]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.3.0

## 1.2.0

### Minor Changes

- 5c88156: Include exports

### Patch Changes

- Updated dependencies [5c88156]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.2.0

## 1.1.0

### Minor Changes

- d755fa3: Minimum viable block.
- 800493f: Updated UMAP/tSNE visualization

### Patch Changes

- Updated dependencies [d755fa3]
- Updated dependencies [800493f]
  - @platforma-open/milaboratories.dimensionality-reduction.model@1.1.0
