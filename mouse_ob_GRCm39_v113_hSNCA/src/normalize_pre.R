# Pre-filter Seurat pipeline: normalize + dim reduction.
#
#   1. Read `outputs/rds/aso_wt_data_annotated.rds`.
#   2. NormalizeData (LogNormalize, scale.factor = 1e4 default).
#   3. FindVariableFeatures (vst, 2000 features default).
#   4. ScaleData on the variable features.
#   5. RunPCA (50 PCs default).
#   6. RunUMAP on PCs 1:30.
#   7. Write `outputs/rds/aso_wt_data_annotated_normalized.rds`.

in_rds <- "outputs/rds/aso_wt_data_annotated.rds"
out_rds <- "outputs/rds/aso_wt_data_annotated_normalized.rds"

paste0("Started running the script...")
library(Seurat)      
dataSeurat <- readRDS(in_rds)

# Standard Seurat pipeline
paste0("Running NormalizeData...")
dataSeurat <- NormalizeData(dataSeurat, verbose = TRUE)

paste0("Running FindVariableFeatures...")
dataSeurat <- FindVariableFeatures(dataSeurat, verbose = TRUE)

paste0("Running ScaleData...")
dataSeurat <- ScaleData(dataSeurat, verbose = TRUE)

paste0("Running RunPCA...")
dataSeurat <- RunPCA(dataSeurat, verbose = TRUE)

paste0("Running RunUMAP...")
dataSeurat <- RunUMAP(dataSeurat, dims = 1:30, verbose = TRUE)

saveRDS(dataSeurat, out_rds)
paste0("wrote:", out_rds, "  (", round(file.size(out_rds) / 1e9, 2), " GB)\n", sep = "")