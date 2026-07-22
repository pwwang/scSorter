# scSorter

> **Fork of [cran/scSorter](https://github.com/hyguo2/scSorter)** — adds Seurat v5 object support via `RunScSorter()`.

scSorter is an R package that implements a semi-supervised algorithm for assigning cells to known cell types in single-cell RNA sequencing (scRNA-seq) data, based on user-provided marker genes.

The method is described in:

> Guo, H., and Li, J. (2021). **scSorter: assigning cells to known cell types according to known marker genes.** *Genome Biology*, 22, 11.

## What's New in This Fork

The original scSorter requires raw expression matrices. This fork adds **native Seurat v5 support** via `RunScSorter()`, so you can run the algorithm directly on a Seurat object without manually extracting and preprocessing the expression matrix.

| Feature | Upstream | This Fork |
|---|---|---|
| Raw matrix input (`scSorter()`) | ✅ | ✅ |
| Seurat v5 object input (`RunScSorter()`) | ❌ | ✅ |
| Automatic variable feature selection | ❌ | ✅ |
| Sets Seurat identity class automatically | ❌ | ✅ |
| Expression filtering by detection rate (`min_pct`) | ❌ | ✅ |
| Allowing all markers to be missed in object for certain cell types | ❌ | ✅ |

## Installation

```r
# install.packages("remotes")
remotes::install_github("pwwang/scSorter")
```

## Quick Start — Seurat (Recommended)

```r
library(Seurat)
library(scSorter)

# Your Seurat object with FindVariableFeatures already run
# anno: data frame with columns "Type" and "Marker"
result_obj <- RunScSorter(
  object = seurat_obj,
  anno = anno,
  top_vf = 2000,     # use top 2000 variable features
  min_pct = 0.1      # gene must be detected in ≥10% of cells
)

# Cell type predictions are stored in the metadata and set as the active identity
table(result_obj$scSorter_celltype)
DimPlot(result_obj, label = TRUE)
```

### `RunScSorter()` Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `object` | A Seurat object (must have `FindVariableFeatures` run) | — |
| `anno` | Marker gene annotation with columns `Type`, `Marker`, plus optional `Weight` | — |
| `layer` | Seurat layer for expression data (`NULL` = default layer) | `NULL` |
| `assay` | Seurat assay to use (`NULL` = default assay) | `NULL` |
| `top_vf` | Number of top variable features to include (`NULL` = all) | `2000` |
| `min_pct` | Minimum fraction of cells expressing a gene to retain it | `0.1` |
| `set_ident` | Set predicted cell types as active Seurat identity | `TRUE` |
| `name` | Metadata column name for predictions | `"scSorter_celltype"` |
| `...` | Additional arguments passed to `scSorter()` | — |


