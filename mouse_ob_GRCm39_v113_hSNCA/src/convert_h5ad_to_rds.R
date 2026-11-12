# Convert merged .h5ad to a Seurat .rds and attach MapMyCells annotations.
#
#   1. Read `outputs/h5ad/aso_wt_data.h5ad` (raw int UMI in `.X`,
#      `hSNCA_transgene` preserved) via anndata.
#   2. Transpose `.X` (cells x genes -> genes x cells) and build a Seurat
#      object with the AnnData `.obs` as cell metadata. Seurat replaces
#      underscores in feature names with dashes -> `hSNCA-transgene`.
#   3. Sanity-check that `hSNCA_transgene` made it through.
#   4. Auto-discover the single MapMyCells CSV under
#      `outputs/mapmycells_csv/v1/` and join its per-cell labels via
#      `AddMetaData`.
#   5. Write `outputs/rds/aso_wt_data_annotated.rds` and print its size.

library(anndata)
library(Seurat)
library(Matrix)

paste0("Started running the script...")


in_h5ad  <- "outputs/h5ad/aso_wt_data.h5ad"
csv_dir  <- "outputs/mapmycells_csv/v1"
# Expect exactly one MapMyCells CSV under csv_dir/<subdir>/*.csv
csv_path <- list.files(csv_dir, pattern = "\\.csv$",
                       recursive = TRUE, full.names = TRUE)
stopifnot(length(csv_path) == 1)
out_rds  <- "outputs/rds/aso_wt_data_annotated.rds"

ad <- read_h5ad(in_h5ad)

# AnnData is cells x genes; Seurat wants genes x cells
counts <- as(Matrix::t(ad$X), "CsparseMatrix")
rownames(counts) <- ad$var_names
colnames(counts) <- rownames(ad$obs)

stopifnot("hSNCA_transgene" %in% rownames(counts))

dataSeurat <- CreateSeuratObject(
  counts    = counts,
  meta.data = as.data.frame(ad$obs),
  project   = "OB_ASO_WT_GRCm39_v113_hSNCA"
)

# --- Merge MapMyCells labels from the single CSV ---
mapping <- read.csv(csv_path, row.names = 1, comment.char = "#")

stopifnot(setequal(rownames(mapping), colnames(dataSeurat)))
dataSeurat <- AddMetaData(dataSeurat, metadata = mapping)

dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(dataSeurat, out_rds)
paste0("wrote:", out_rds, "  (", round(file.size(out_rds) / 1e9, 2), " GB)\n", sep = "")
