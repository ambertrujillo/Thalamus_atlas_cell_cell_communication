# Cell-cell communication Pipeline - Thalamus atlas
## Overview

This pipeline infers and compares cell-cell communication networks between two conditions
(e.g., FTD vs. Control) using `CellChat`. It consists of three sequential scripts:

1. **`CellChat_FTD.R`** — runs CellChat on the FTD dataset, inferring ligand-receptor (L-R)
   interactions and aggregating into pathway-level communication probabilities
2. **`CellChat_Control.R`** — runs the same pipeline on the Control dataset
3. **`CellChat_Merge.R`** — merges the two CellChat objects, computes differential
   interactions, and produces comparative visualizations

Scripts 1 and 2 are structurally identical and only differ in their input `.h5ad` file —
they're kept as separate scripts (rather than parameterized) so each condition has a clean
audit trail and intermediate `.rds` output.

---

## Pipeline Flow

```
FTD_adata.h5ad ───────► CellChat_FTD.R ───────► cellChat_FTD.rds
                                                       │
Control_adata.h5ad ───► CellChat_Control.R ──► cellChat_Control.rds
                                                       │
                                                       ▼
                                          CellChat_Merge.R
                                                       │
                                                       ▼
                              cellchat_merged_FTD_vs_Control.RData
                              cellchat_object.list_FTD_vs_Control.RData
                              heatmap_cellchat_FTDvControl_strength.svg
                              GABA_InNeu.pdf
                              (differential interaction plots)
```

---

## Dependencies

```{r dependencies, eval=FALSE}
# Bioconductor
BiocManager::install("ComplexHeatmap")

# CRAN
install.packages(c("anndata", "Seurat", "dplyr", "patchwork", "ggplot2", "svglite"))

# CellChat from GitHub
devtools::install_github("jinworks/CellChat")
```

`ComplexHeatmap` is a required dependency for CellChat's `netVisual_heatmap()` but does not
need an explicit `library()` call — it loads internally with CellChat.

---

## Script 1: `CellChat_FTD.R`

### Purpose

Builds a CellChat object from the FTD dataset, identifies overexpressed ligand-receptor
interactions, and infers cell-cell communication probabilities at both the L-R pair level
and the pathway level.

### Inputs

| File | Description |
|---|---|
| `FTD_adata.h5ad` | AnnData object containing FTD cells. Must include raw counts in `layers["counts"]` and cell-level metadata in `obs` with at least `SampleID` and an `annotation` column |

### User Configuration

The script assumes the cell type annotation column in `obs` is named `annotation`. Adjust
the `group.by` argument and the `colnames(metadata)` rename step if your column names differ.

For mouse data, change:

```{r mouse-db, eval=FALSE}
CellChatDB <- CellChatDB.human
```

to:

```{r mouse-db-2, eval=FALSE}
CellChatDB <- CellChatDB.mouse
```

### Pipeline Steps

1. Load AnnData and extract raw counts from `layers["counts"]`
2. Log-normalize counts with `Seurat::NormalizeData(scale.factor = 10000)`
3. Build CellChat object with cell annotations as identities
4. Set CellChat database (`CellChatDB.human` by default)
5. Identify overexpressed genes and L-R interactions
6. Compute communication probabilities (L-R pair level)
7. Compute pathway-level communication probabilities
8. Aggregate the cell-cell communication network

### Outputs

| File | Description |
|---|---|
| `cellChat_FTD.rds` | Saved CellChat object with inferred communication networks for the FTD dataset |

### Dependencies

```
library(anndata)
library(Seurat)
library(CellChat)
library(dplyr)
```

---

## Script 2: `CellChat_Control.R`

### Purpose

Identical pipeline to Script 1, applied to the Control dataset.

### Inputs

| File | Description |
|---|---|
| `Control_adata.h5ad` | AnnData object containing Control cells, with the same structure as the FTD input |

### Outputs

| File | Description |
|---|---|
| `cellChat_Control.rds` | Saved CellChat object for the Control dataset |

### Dependencies

Same as Script 1.

### Note on Keeping Scripts Separate

`CellChat_FTD.R` and `CellChat_Control.R` are intentionally kept as two separate scripts
rather than parameterized into one. This:

- Allows each condition to be run independently (e.g., on different cluster nodes)
- Preserves an audit trail per condition
- Makes the merge step in Script 3 robust to one condition being re-run without re-running
  the other

---

## Script 3: `CellChat_Merge.R`

### Purpose

Loads the two CellChat objects produced by Scripts 1 and 2, merges them for comparative
analysis, and generates differential interaction visualizations including aggregate
comparisons, differential interaction circle plots, an interaction strength heatmap, and
pathway-specific circle plots.

### Inputs

| File | Description |
|---|---|
| `cellChat_FTD.rds` | Output of Script 1 |
| `cellChat_Control.rds` | Output of Script 2 |

### Outputs

| File | Description |
|---|---|
| `cellchat_object.list_FTD_vs_Control.RData` | List object containing both individual CellChat objects, for downstream reuse |
| `cellchat_merged_FTD_vs_Control.RData` | Merged CellChat object for direct comparative analysis |
| `heatmap_cellchat_FTDvControl_strength.svg` | Heatmap showing differential interaction strength across cell type pairs |
| `GABA_InNeu.pdf` | Pathway-specific (GABA-A) circle plot showing inhibitory neuron communication in each condition |
| Inline plots (Rstudio) | Aggregate interaction count/strength bar plots and differential interaction circle plots per source cell type |

### Key Comparative Visualizations

**1. `compareInteractions()`** — bar plot of total interaction count and strength,
side by side for the two conditions.

**2. `netVisual_diffInteraction()`** — circle plot showing the difference in interaction
strength between conditions, with red/blue edges denoting up/down regulation. Run separately
per source cell type of interest (the script runs it for the four neuronal subtypes).

**3. `netVisual_heatmap()`** — comprehensive heatmap of differential interaction strength
across all cell type pairs, with rows and columns ordered by `cell_order`.

**4. `netVisual_aggregate()` per pathway** — pathway-specific (e.g. GABA-A) circle plots
side by side for the two conditions, with `edge.weight.max` fixed across both plots for
visual comparability.

### Visualization Comparability

`getMaxWeight(object.list, ...)` is used to compute the **global maximum edge weight**
across the two conditions, then passed to `edge.weight.max` in every subsequent plot. This
ensures that edge thickness encodes the same absolute scale across both conditions, making
differences in edge thickness directly interpretable.

### Dependencies

```
library(CellChat)
library(patchwork)
library(ggplot2)
library(svglite)
```

---

## Directory Structure

```
cellchat_pipeline/
├── README.Rmd
    ├── scripts
        ├── CellChat_FTD.R
        ├── CellChat_Control.R
        ├── CellChat_Merge.R

    ├── outputs/
        ├── cellChat_FTD.rds
        ├── cellChat_Control.rds
        ├── cellchat_object.list_FTD_vs_Control.RData
        ├── cellchat_merged_FTD_vs_Control.RData
        ├── heatmap_cellchat_FTDvControl_strength.svg
        └── GABA_InNeu.pdf
```

---

## Running the Full Pipeline

```{bash run-pipeline, eval=FALSE}
Rscript CellChat_FTD.R       # produces cellChat_FTD.rds
Rscript CellChat_Control.R   # produces cellChat_Control.rds
Rscript CellChat_Merge.R     # merges and generates comparative plots
```

Scripts 1 and 2 are independent and can be run in parallel. Script 3 requires both `.rds`
outputs to exist before running.

---

## Statistical and Methodological Notes

- **Communication probability**: CellChat uses a mass action model based on average
  expression of ligand (sender) and receptor subunits (receiver). Significance is assessed
  by permutation testing.
- **Interaction strength vs. interaction number**: "Strength" refers to the sum of
  communication probabilities (magnitude of signaling), while "number" refers to the count
  of statistically significant L-R pairs (diversity of signaling). Both measures can tell
  different biological stories and are often examined together.
- **Visual comparability**: Comparability across merged conditions is achieved by
  normalizing edge weights to a **global maximum** via `getMaxWeight()` and passing this
  to `edge.weight.max` in every plot. Without this, each condition's edges would be scaled
  independently and not directly comparable.
- **Mass action model assumption**: CellChat assumes that communication probability is
  proportional to the product of ligand and receptor expression. This treats expression
  as a proxy for signaling activity and does not account for post-translational regulation
  or spatial constraints.

---

## Common Issues

| Error | Likely Cause | Fix |
|---|---|---|
| `count_matrix dimensions don't match metadata` | Cells were filtered between count extraction and metadata extraction | Confirm `identical(colnames(norm_count_matrix), rownames(metadata))` returns `TRUE` before creating the CellChat object |
| `annotation column not found` | The metadata column name differs from `"annotation"` | Update `group.by` and `setIdent(ident.use = ...)` calls to match your actual column name |
| `CellChatDB.human not found` | CellChat installed incompletely | Reinstall via `devtools::install_github("jinworks/CellChat")` ensuring all dependencies install successfully |
| Pathway shown is empty in one condition | That pathway has no significant interactions in that condition | Verify with `cellchat@netP$pathways` per condition; if absent, that condition lacks signaling for that pathway |
| Heatmap throws cell order error | Cell types in `cell_order` don't match identity levels exactly | Run `levels(cellchat@idents)` and confirm spelling/spacing matches `cell_order` |

---

## Citation

Jin S, Plikus MV, Nie Q (2025). *CellChat for systematic analysis of cell-cell
communication from single-cell transcriptomics.* Nature Protocols. doi:10.1038/s41596-024-01045-4