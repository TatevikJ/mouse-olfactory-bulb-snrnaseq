# QC filtering v3. Basic-QC (v1, doublet < 0.2) + brain-region filter.
# Reads the v1 filtered RDS (already passes nCount < 400k, scDblFinder < 0.2,
# percent.mt < 0.5) and drops nuclei whose supertype is in the
# TO_REMOVE sheet of meta/qc_post_v1_region_list.xlsx (non-OB supertypes).
# Writes the v3 filtered Seurat .rds + h5ad + per-criterion summary TSVs.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(anndata)
  library(readxl)
})

# ---- Inputs / outputs ----
in_rds       <- "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v1.rds"
in_region    <- "meta/qc_post_v1_region_list.xlsx"
out_rds      <- "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v3.rds"
out_h5ad     <- "outputs/h5ad/aso_wt_data_annotated_normalized_qc_filtered_v3.h5ad"
tsv_dir      <- "outputs/qc_post_v3/tsv"
dir.create(dirname(out_rds),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_h5ad), recursive = TRUE, showWarnings = FALSE)
dir.create(tsv_dir,           recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
cat("Loading:", in_rds, "\n")
dataSeurat <- readRDS(in_rds)
n_initial  <- ncol(dataSeurat)
cat("Initial nuclei (post-basic-QC v1):", n_initial,
    " Genes:", nrow(dataSeurat), "\n")

# ---- Load brain-region removal list (TO_REMOVE sheet) ----
cat("Loading TO_REMOVE supertypes:", in_region, "\n")
remove_df <- suppressMessages(readxl::read_excel(in_region, sheet = "TO_REMOVE"))
remove_supertypes <- unique(remove_df$supertype_name[!is.na(remove_df$supertype_name)])
cat("Supertypes to remove:", length(remove_supertypes), "\n")

# ---- Per-nucleus filter status ----
dataSeurat$filter_status <- "pass"
dataSeurat$filter_status[dataSeurat$supertype_name %in% remove_supertypes] <-
  "removed_region"

# ---- Apply region filter ----
cat("\nApplying brain-region filter...\n")
dataSeurat_filtered <- subset(
  dataSeurat,
  subset = !(supertype_name %in% remove_supertypes)
)
n_final <- ncol(dataSeurat_filtered)
cat("Nuclei after region filter:", n_final, "\n")
cat("Removed:", n_initial - n_final,
    " (", round((n_initial - n_final) / n_initial * 100, 2), "%)\n", sep = "")

# ---- Per-criterion summary ----
filter_summary <- data.frame(
  criterion = c(
    "supertype not in TO_REMOVE (region filter)",
    "all criteria combined"
  ),
  n_nuclei_passing = c(
    sum(!(dataSeurat$supertype_name %in% remove_supertypes)),
    n_final
  ),
  n_nuclei_removed = c(
    sum(dataSeurat$supertype_name %in% remove_supertypes),
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

# ---- Per-removed-supertype tally ----
removed_per_supertype <- dataSeurat@meta.data %>%
  filter(supertype_name %in% remove_supertypes) %>%
  count(class_name, subclass_name, supertype_name, name = "n_nuclei_removed") %>%
  arrange(desc(n_nuclei_removed))
write.table(removed_per_supertype,
            file.path(tsv_dir, "removed_per_supertype.tsv"),
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
