# QC filtering v2. Basic-QC only (lenient doublet threshold 0.3 vs v1's 0.2):
#   keep: nCount_RNA < 400,000, scDblFinder.score < 0.3, percent.mt < 0.5
# Brain-region (supertype) filtering is deferred to v3 (v1 + region) and v4
# (v2 + region). Writes the filtered Seurat .rds + h5ad. Gene-level (>= 3 nuclei)
# filter is already applied upstream in merge_count_matrices.py.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(anndata)
})

# ---- Filter thresholds ----
thresh_scDblFinder_max <- 0.3
thresh_percent_mt_max  <- 0.5
thresh_nCount_RNA_max  <- 400000

# ---- Inputs / outputs ----
in_rds   <- "outputs/rds/aso_wt_data_annotated_normalized.rds"
out_rds  <- "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v2.rds"
out_h5ad <- "outputs/h5ad/aso_wt_data_annotated_normalized_qc_filtered_v2.h5ad"
tsv_dir  <- "outputs/qc_post_v2/tsv"
dir.create(dirname(out_rds),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_h5ad), recursive = TRUE, showWarnings = FALSE)
dir.create(tsv_dir,           recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
cat("Loading:", in_rds, "\n")
dataSeurat <- readRDS(in_rds)
n_initial  <- ncol(dataSeurat)
cat("Initial nuclei:", n_initial, " Genes:", nrow(dataSeurat), "\n")

# ---- Per-nucleus filter status (priority: doublet > mt > nCount) ----
dataSeurat$filter_status <- "pass"
dataSeurat$filter_status[dataSeurat$scDblFinder.score >= thresh_scDblFinder_max] <- "removed_doublet"
dataSeurat$filter_status[dataSeurat$percent.mt >= thresh_percent_mt_max &
                         dataSeurat$filter_status == "pass"] <- "removed_mt"
dataSeurat$filter_status[dataSeurat$nCount_RNA >= thresh_nCount_RNA_max &
                         dataSeurat$filter_status == "pass"] <- "removed_ncount"

# ---- Apply filters ----
cat("\nApplying filters...\n")
dataSeurat_filtered <- subset(
  dataSeurat,
  subset = scDblFinder.score < thresh_scDblFinder_max &
           percent.mt        < thresh_percent_mt_max &
           nCount_RNA        < thresh_nCount_RNA_max
)
n_final <- ncol(dataSeurat_filtered)
cat("Nuclei after filtering:", n_final, "\n")
cat("Removed:", n_initial - n_final,
    " (", round((n_initial - n_final) / n_initial * 100, 2), "%)\n", sep = "")

# ---- Per-criterion summary ----
filter_summary <- data.frame(
  criterion = c(
    sprintf("scDblFinder.score < %s", thresh_scDblFinder_max),
    sprintf("percent.mt < %s",        thresh_percent_mt_max),
    sprintf("nCount_RNA < %s",        thresh_nCount_RNA_max),
    "all criteria combined"
  ),
  n_nuclei_passing = c(
    sum(dataSeurat$scDblFinder.score < thresh_scDblFinder_max),
    sum(dataSeurat$percent.mt        < thresh_percent_mt_max),
    sum(dataSeurat$nCount_RNA        < thresh_nCount_RNA_max),
    n_final
  ),
  n_nuclei_removed = c(
    sum(dataSeurat$scDblFinder.score >= thresh_scDblFinder_max),
    sum(dataSeurat$percent.mt        >= thresh_percent_mt_max),
    sum(dataSeurat$nCount_RNA        >= thresh_nCount_RNA_max),
    n_initial - n_final
  )
)
print(filter_summary)
write.table(filter_summary, file.path(tsv_dir, "filter_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- Per-sample filter_status tally ----
per_sample_status <- dataSeurat@meta.data %>%
  count(sample, condition, filter_status, name = "n_nuclei") %>%
  arrange(condition, sample, filter_status)
write.table(per_sample_status,
            file.path(tsv_dir, "filter_status_per_sample.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ---- Save filtered RDS ----
cat("\nSaving filtered RDS:", out_rds, "\n")
saveRDS(dataSeurat_filtered, out_rds)

# ---- Export to h5ad (raw counts) ----
cat("Exporting h5ad:", out_h5ad, "\n")
counts_matrix <- dataSeurat_filtered@assays$RNA$counts
obs_meta <- dataSeurat_filtered@meta.data
var_meta <- data.frame(gene = rownames(counts_matrix),
                       row.names = rownames(counts_matrix))
ad <- AnnData(X = Matrix::t(counts_matrix), obs = obs_meta, var = var_meta)
write_h5ad(ad, out_h5ad, compression = "gzip", compression_opts = 7)
cat("h5ad file size:", round(file.size(out_h5ad) / (1024^2), 2), "MB\n")

cat("\nDone.\n")
