# Post-filter normalization, clustering, and marker discovery (v4).
# Re-runs the standard Seurat pipeline on the v4 filtered RDS (subsetting in
# qc_filtering_v4.R invalidated the HVGs / scale.data / PCA / UMAP from
# normalize_pre.R), then adds Leiden clustering at multiple resolutions and
# FindAllMarkers for the low-resolution sweeps and the Allen taxonomy levels
# (class, subclass, supertype). A protein-coding-only variant
# of the final object is written for downstream DEG / GSEA. The hSNCA
# transgene (dash form) is explicitly retained alongside the protein-coding
# set because org.Mm.eg.db does not classify it.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

# ---- Inputs / outputs ----
in_rds       <- "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v4.rds"
out_rds      <- "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v4_clustered.rds"
out_rds_pc   <- "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v4_clustered_protein_coding.rds"
markers_dir  <- "outputs/qc_post_v4_normalized/tsv"
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
dir.create(markers_dir,      recursive = TRUE, showWarnings = FALSE)

# ---- Resolution sweeps ----
cluster_resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0)
markers_resolutions <- c(0.1, 0.2, 0.3, 0.4)
markers_levels      <- c("class_name", "subclass_name", "supertype_name")

# ---- Load v4 ----
cat("Loading:", in_rds, "\n")
dataSeurat <- readRDS(in_rds)
cat("Nuclei:", ncol(dataSeurat), " Genes:", nrow(dataSeurat), "\n")

# ---- Standard Seurat pipeline (re-run on the v4 subset) ----
cat("\n[1/5] NormalizeData...\n")
dataSeurat <- NormalizeData(dataSeurat, verbose = FALSE)

cat("[2/5] FindVariableFeatures...\n")
dataSeurat <- FindVariableFeatures(dataSeurat, verbose = FALSE)

cat("[3/5] ScaleData...\n")
dataSeurat <- ScaleData(dataSeurat, verbose = FALSE)

cat("[4/5] RunPCA...\n")
dataSeurat <- RunPCA(dataSeurat, verbose = FALSE)

cat("[5/5] RunUMAP (dims=1:30)...\n")
dataSeurat <- RunUMAP(dataSeurat, dims = 1:30, verbose = FALSE)

# ---- Clustering: Leiden over a resolution sweep ----
cat("\nFindNeighbors (dims=1:30)...\n")
dataSeurat <- FindNeighbors(dataSeurat, dims = 1:30, verbose = FALSE)

cat("FindClusters (Leiden, algorithm = 4) over",
    length(cluster_resolutions), "resolutions...\n")
for (res in cluster_resolutions) {
  cat("  resolution =", res, "\n")
  dataSeurat <- FindClusters(dataSeurat,
                             resolution = res,
                             algorithm   = 4,
                             method      = "igraph",
                             verbose     = FALSE)
}

# Summary of cluster counts per resolution
res_cols <- grep("^RNA_snn_res\\.", colnames(dataSeurat@meta.data), value = TRUE)
cat("\nClusters per resolution:\n")
for (col in res_cols) {
  cat("  ", col, ": ", length(unique(dataSeurat@meta.data[[col]])), "\n", sep = "")
}

# ---- Save clustered RDS ----
cat("\nSaving clustered RDS:", out_rds, "\n")
saveRDS(dataSeurat, out_rds)

# ---- FindAllMarkers at low-resolution sweeps ----
cat("\nRunning FindAllMarkers for low-resolution sweeps...\n")
for (res in markers_resolutions) {
  cluster_col <- sprintf("RNA_snn_res.%s", res)
  if (!cluster_col %in% colnames(dataSeurat@meta.data)) {
    cat("  Skipping res =", res, "; column", cluster_col, "not found.\n")
    next
  }
  Idents(dataSeurat) <- dataSeurat@meta.data[[cluster_col]]

  cat("  resolution =", res, "(", cluster_col, ")\n")
  markers <- FindAllMarkers(dataSeurat,
                            only.pos        = TRUE,
                            min.pct         = 0.25,
                            logfc.threshold = 0.25,
                            verbose         = FALSE)
  out_csv <- file.path(markers_dir, sprintf("markers_res_%s.csv", res))
  write.csv(markers, out_csv, row.names = FALSE)
  cat("    found ", nrow(markers), " markers; wrote ", out_csv, "\n", sep = "")

  top10 <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = 10) %>%
    dplyr::arrange(cluster, dplyr::desc(avg_log2FC))
  out_csv_top10 <- file.path(markers_dir, sprintf("markers_top10_res_%s.csv", res))
  write.csv(top10, out_csv_top10, row.names = FALSE)
  cat("    wrote ", out_csv_top10, "\n", sep = "")
}

# ---- FindAllMarkers at Allen taxonomy levels ----
cat("\nRunning FindAllMarkers for taxonomy levels...\n")
for (level in markers_levels) {
  Idents(dataSeurat) <- dataSeurat@meta.data[[level]]

  cat("  level =", level, "\n")
  markers <- FindAllMarkers(dataSeurat,
                            only.pos        = TRUE,
                            min.pct         = 0.25,
                            logfc.threshold = 0.25,
                            verbose         = FALSE)
  out_csv <- file.path(markers_dir, sprintf("markers_%s.csv", level))
  write.csv(markers, out_csv, row.names = FALSE)
  cat("    found ", nrow(markers), " markers; wrote ", out_csv, "\n", sep = "")

  top10 <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = 10) %>%
    dplyr::arrange(cluster, dplyr::desc(avg_log2FC))
  out_csv_top10 <- file.path(markers_dir, sprintf("markers_top10_%s.csv", level))
  write.csv(top10, out_csv_top10, row.names = FALSE)
  cat("    wrote ", out_csv_top10, "\n", sep = "")
}

# ---- Protein-coding-only subset ----
# Filter to protein-coding genes via org.Mm.eg.db, then UNION the hSNCA
# transgene (dash form post-CreateSeuratObject) because the annotation DB
# does not classify it as protein-coding.
cat("\nBuilding protein-coding subset (org.Mm.eg.db + hSNCA transgene)...\n")
all_genes <- rownames(dataSeurat)
cat("  Total genes in object:", length(all_genes), "\n")

gene_info <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys      = all_genes,
  columns   = c("SYMBOL", "GENETYPE"),
  keytype   = "SYMBOL"
)
gene_info <- gene_info[!duplicated(gene_info$SYMBOL), ]
cat("  Gene-type distribution:\n")
print(table(gene_info$GENETYPE, useNA = "ifany"))

protein_coding_genes <- gene_info %>%
  dplyr::filter(GENETYPE == "protein-coding") %>%
  dplyr::pull(SYMBOL)
cat("  Protein-coding gene symbols:", length(protein_coding_genes), "\n")

genes_to_keep <- intersect(all_genes, protein_coding_genes)
if ("hSNCA-transgene" %in% all_genes) {
  genes_to_keep <- union(genes_to_keep, "hSNCA-transgene")
  cat("  hSNCA-transgene retained.\n")
} else {
  warning("hSNCA-transgene not found in object rownames; check upstream merge.")
}
cat("  Genes kept in protein-coding subset:", length(genes_to_keep), "\n")

dataSeurat_pc <- dataSeurat[genes_to_keep, ]
cat("  Object dims after PC filter: ",
    nrow(dataSeurat_pc), " genes x ",
    ncol(dataSeurat_pc), " nuclei\n", sep = "")

cat("\nSaving protein-coding RDS:", out_rds_pc, "\n")
saveRDS(dataSeurat_pc, out_rds_pc)

cat("\nDone.\n")
