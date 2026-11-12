# Post-normalization diagnostics (v4). Reads the v4 clustered protein-coding
# RDS produced by normalize_post_v4.R and writes UMAPs (sample, condition,
# class, subclass, supertype, clusters at each resolution) plus per-cluster
# composition TSVs to outputs/qc_post_v4_normalized/{plots,tsv}/. The marker CSVs are
# already produced by normalize_post_v4.R and not regenerated here.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# ---- Inputs / outputs ----
in_rds    <- "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v4_clustered_protein_coding.rds"
plots_dir <- "outputs/qc_post_v4_normalized/plots"
tsv_dir   <- "outputs/qc_post_v4_normalized/tsv"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsv_dir,   recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
cat("Loading:", in_rds, "\n")
dataSeurat <- readRDS(in_rds)
if (!"condition" %in% colnames(dataSeurat@meta.data)) {
  dataSeurat$condition <- ifelse(grepl("^ASO", dataSeurat$sample), "ASO", "WT")
}
cat("Nuclei:", ncol(dataSeurat), " Genes:", nrow(dataSeurat), "\n")

# Condition palette: WT = teal, ASO = red
condition_colors <- c(WT = "#59b2aa", ASO = "#d43f2e")

# ---- Plot themes ----
project_base_theme <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
}

project_umap_theme <- function(legend_corner = c("top-left")) {
  pos <- switch(legend_corner,
                "top-left"    = list(p = c(0.02, 0.98), j = c(0, 1)),
                "top-right"   = list(p = c(0.98, 0.98), j = c(1, 1)),
                "bottom-right"= list(p = c(0.98, 0.02), j = c(1, 0)),
                "bottom-left" = list(p = c(0.02, 0.02), j = c(0, 0)))
  theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text  = element_blank(),
      aspect.ratio = 1,
      plot.title = element_text(size = 13, hjust = 0),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.position = pos$p,
      legend.justification = pos$j,
      legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
      legend.box.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
}

save_pdf <- function(path, plot, width, height) {
  cat("  saving:", basename(path), sprintf("(%g x %g)\n", width, height))
  ggsave(path, plot = plot, width = width, height = height)
}

# ============================================================
# UMAP plots: sample, condition, class, subclass, supertype
# ============================================================
cat("\n[1/3] UMAPs (sample / condition / class / subclass / supertype)...\n")

p_umap_sample <- DimPlot(dataSeurat, reduction = "umap",
                         group.by = "sample") +
  labs(x = "UMAP1", y = "UMAP2", title = "Sample") +
  project_umap_theme("top-left") +
  theme(legend.position = "right",
        legend.background = element_blank(),
        legend.box.background = element_blank())
save_pdf(file.path(plots_dir, "umap_sample.pdf"),
         plot = p_umap_sample, width = 8, height = 6)

p_umap_condition <- DimPlot(dataSeurat, reduction = "umap",
                            group.by = "condition") +
  scale_color_manual(values = condition_colors) +
  labs(x = "UMAP1", y = "UMAP2", title = "Condition") +
  project_umap_theme("top-left") +
  theme(legend.position = "right",
        legend.background = element_blank(),
        legend.box.background = element_blank())
save_pdf(file.path(plots_dir, "umap_condition.pdf"),
         plot = p_umap_condition, width = 7, height = 6)

collapse_to_top_n <- function(x, n = 15, other_label = "other") {
  x <- as.character(x)
  freq <- sort(table(x[!is.na(x)]), decreasing = TRUE)
  top <- names(freq)[seq_len(min(n, length(freq)))]
  ifelse(is.na(x) | x %in% top, x, other_label)
}

if ("class_name" %in% colnames(dataSeurat@meta.data) &&
    !"class_new" %in% colnames(dataSeurat@meta.data)) {
  dataSeurat$class_new <- collapse_to_top_n(dataSeurat$class_name, n = 15)
}
if ("subclass_name" %in% colnames(dataSeurat@meta.data) &&
    !"subclass_new" %in% colnames(dataSeurat@meta.data)) {
  dataSeurat$subclass_new <- collapse_to_top_n(dataSeurat$subclass_name, n = 15)
}
if ("supertype_name" %in% colnames(dataSeurat@meta.data) &&
    !"supertype_new" %in% colnames(dataSeurat@meta.data)) {
  dataSeurat$supertype_new <- collapse_to_top_n(dataSeurat$supertype_name, n = 20)
}

to_ordered_factor <- function(x) {
  x <- as.character(x)
  lvls <- sort(setdiff(unique(x), c(NA, "other")))
  if ("other" %in% x) lvls <- c(lvls, "other")
  factor(x, levels = lvls)
}
for (col in intersect(c("class_name", "class_new", "subclass_name", "subclass_new",
                        "supertype_new"),
                      colnames(dataSeurat@meta.data))) {
  dataSeurat@meta.data[[col]] <- to_ordered_factor(dataSeurat@meta.data[[col]])
}

add_centroid_labels <- function(p, group_vec) {
  emb <- dataSeurat@reductions$umap@cell.embeddings
  centroids <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                          label = as.character(group_vec)) %>%
    dplyr::filter(!is.na(label)) %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(x = median(UMAP_1), y = median(UMAP_2), .groups = "drop")
  p + geom_text(data = centroids, aes(x = x, y = y, label = label),
                inherit.aes = FALSE,
                size = 2.5, color = "black",
                bg.color = "white", bg.r = 0.15)
}

for (col in intersect(c("class_name", "class_new"), colnames(dataSeurat@meta.data))) {
  p <- DimPlot(dataSeurat, reduction = "umap", group.by = col, label = FALSE) +
    labs(x = "UMAP1", y = "UMAP2",
         title = paste0("Class (", col, if (grepl("_new$", col)) ", top_n = 15" else "", ")")) +
    project_umap_theme("top-left") +
    NoLegend()
  p <- add_centroid_labels(p, dataSeurat@meta.data[[col]])
  save_pdf(file.path(plots_dir, paste0("umap_", col, ".pdf")),
           plot = p, width = 6, height = 6)
}

for (col in intersect(c("subclass_name", "subclass_new"), colnames(dataSeurat@meta.data))) {
  p <- DimPlot(dataSeurat, reduction = "umap", group.by = col,
               label = TRUE, repel = TRUE, label.size = 2.5) +
    labs(x = "UMAP1", y = "UMAP2",
         title = paste0("Subclass (", col, if (grepl("_new$", col)) ", top_n = 15" else "", ")")) +
    project_umap_theme("top-left") +
    NoLegend()
  save_pdf(file.path(plots_dir, paste0("umap_", col, ".pdf")),
           plot = p, width = 6, height = 6)
}

# Supertypes: only the top-20-collapsed supertype_new is plotted (too many
# levels otherwise), once with on-plot labels and once with an outside legend.
if ("supertype_new" %in% colnames(dataSeurat@meta.data)) {
  p <- DimPlot(dataSeurat, reduction = "umap", group.by = "supertype_new",
               label = TRUE, repel = TRUE, label.size = 2.5) +
    labs(x = "UMAP1", y = "UMAP2", title = "Supertype (supertype_new, top_n = 20)") +
    project_umap_theme("top-left") +
    NoLegend()
  save_pdf(file.path(plots_dir, "umap_supertype_new.pdf"),
           plot = p, width = 6, height = 6)

  p_leg <- DimPlot(dataSeurat, reduction = "umap", group.by = "supertype_new",
                   label = TRUE, repel = TRUE, label.size = 2.5) +
    labs(x = "UMAP1", y = "UMAP2", title = "Supertype (supertype_new, top_n = 20)") +
    project_umap_theme("top-left") +
    theme(legend.position = "right",
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "lines")) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 1.5)))
  save_pdf(file.path(plots_dir, "umap_supertype_new_legend.pdf"),
           plot = p_leg, width = 8, height = 6)
}

# hSNCA transgene expression (dash form post-CreateSeuratObject).
if ("hSNCA-transgene" %in% rownames(dataSeurat)) {
  # UMAP: grey -> color sequential ramp (standard for gene expression);
  # order = TRUE so expressing nuclei draw on top of the grey non-expressing bulk.
  p_umap_hsnca <- FeaturePlot(dataSeurat, features = "hSNCA-transgene",
                              reduction = "umap", order = TRUE) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    labs(x = "UMAP1", y = "UMAP2", title = "hSNCA-transgene") +
    project_umap_theme("top-left")
  save_pdf(file.path(plots_dir, "umap_hSNCA_transgene.pdf"),
           plot = p_umap_hsnca, width = 6, height = 6)

  # Violins across cell types, split by condition (WT carries no transgene).
  # Log-normalized expression; scale = "width" keeps zero-inflated violins legible.
  expr_df <- data.frame(
    expr         = FetchData(dataSeurat, vars = "hSNCA-transgene", layer = "data")[, 1],
    condition    = factor(dataSeurat$condition, levels = c("WT", "ASO")),
    class_name   = dataSeurat$class_name,
    subclass_new = dataSeurat$subclass_new
  )
  make_hsnca_violin <- function(group_col) {
    ggplot(expr_df, aes(x = .data[[group_col]], y = expr, fill = condition)) +
      geom_violin(scale = "width", linewidth = 0.3) +
      scale_fill_manual(values = condition_colors) +
      labs(x = NULL, y = "hSNCA-transgene (log-norm)",
           title = paste0("hSNCA-transgene expression by ", group_col)) +
      project_base_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  save_pdf(file.path(plots_dir, "violin_hSNCA_transgene_by_class_name.pdf"),
           plot = make_hsnca_violin("class_name"), width = 7, height = 6)
  save_pdf(file.path(plots_dir, "violin_hSNCA_transgene_by_subclass_new.pdf"),
           plot = make_hsnca_violin("subclass_new"), width = 10, height = 6)
} else {
  warning("hSNCA-transgene not found in object rownames; skipping transgene plots.")
}

# ============================================================
# UMAP plots: clusters at each resolution
# ============================================================
cat("\n[2/3] UMAPs per cluster resolution...\n")

res_cols <- grep("^RNA_snn_res\\.", colnames(dataSeurat@meta.data), value = TRUE)
if (length(res_cols) == 0) {
  warning("No RNA_snn_res.* columns found; skipping cluster UMAPs.")
} else {
  for (col in res_cols) {
    res_tag <- sub("^RNA_snn_res\\.", "", col)
    dataSeurat@meta.data[[col]] <- factor(dataSeurat@meta.data[[col]])
    p <- DimPlot(dataSeurat, reduction = "umap", group.by = col,
                 label = TRUE, label.size = 3) +
      labs(x = "UMAP1", y = "UMAP2",
           title = sprintf("Leiden clusters (resolution = %s)", res_tag)) +
      project_umap_theme("top-left") +
      NoLegend()
    save_pdf(file.path(plots_dir, sprintf("umap_clusters_res_%s.pdf", res_tag)),
             plot = p, width = 7, height = 6)
  }
}

# ============================================================
# Cluster composition TSVs
# ============================================================
cat("\n[3/3] Cluster composition TSVs...\n")

# Per (resolution, cluster, condition) tally + WT/ASO ratio
for (col in res_cols) {
  res_tag <- sub("^RNA_snn_res\\.", "", col)
  comp <- dataSeurat@meta.data %>%
    dplyr::count(.data[[col]], condition, name = "n_nuclei") %>%
    tidyr::pivot_wider(names_from = condition, values_from = n_nuclei,
                       values_fill = 0) %>%
    dplyr::rename(cluster = !!col)
  if (all(c("ASO", "WT") %in% colnames(comp))) {
    comp <- comp %>%
      dplyr::mutate(total = ASO + WT,
                    aso_fraction = ASO / total)
  }
  out_tsv <- file.path(tsv_dir,
                       sprintf("cluster_composition_res_%s.tsv", res_tag))
  write.table(comp, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  wrote:", out_tsv, "\n")
}

# Per (cluster, class, subclass) tally at each resolution (helps map clusters
# to Allen taxonomy)
class_cols <- intersect(c("class_name", "subclass_name"),
                        colnames(dataSeurat@meta.data))
for (col in res_cols) {
  res_tag <- sub("^RNA_snn_res\\.", "", col)
  group_cols <- c(col, class_cols)
  comp <- dataSeurat@meta.data %>%
    dplyr::count(dplyr::across(dplyr::all_of(group_cols)), name = "n_nuclei") %>%
    dplyr::rename(cluster = !!col) %>%
    dplyr::arrange(cluster, dplyr::desc(n_nuclei))
  out_tsv <- file.path(tsv_dir,
                       sprintf("cluster_taxonomy_res_%s.tsv", res_tag))
  write.table(comp, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  wrote:", out_tsv, "\n")
}

# Per (taxonomy level, condition) tally + ASO fraction, to check whether any
# Allen group (class / subclass / supertype) is enriched in ASO vs WT.
tax_levels <- intersect(c("class_name", "subclass_name", "supertype_name"),
                        colnames(dataSeurat@meta.data))
for (lvl in tax_levels) {
  comp <- dataSeurat@meta.data %>%
    dplyr::count(.data[[lvl]], condition, name = "n_nuclei") %>%
    tidyr::pivot_wider(names_from = condition, values_from = n_nuclei,
                       values_fill = 0)
  if (all(c("ASO", "WT") %in% colnames(comp))) {
    comp <- comp %>%
      dplyr::mutate(total = ASO + WT, aso_fraction = ASO / total) %>%
      dplyr::arrange(dplyr::desc(aso_fraction))
  }
  out_tsv <- file.path(tsv_dir, sprintf("composition_by_genotype_%s.tsv", lvl))
  write.table(comp, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  wrote:", out_tsv, "\n")
}

# Per (sample, condition) nuclei count on the post-PC object
qc_summary_per_sample <- dataSeurat@meta.data %>%
  dplyr::group_by(sample, condition) %>%
  dplyr::summarise(n_nuclei = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(condition, sample)
write.table(qc_summary_per_sample,
            file.path(tsv_dir, "nuclei_counts_per_sample.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  wrote:", file.path(tsv_dir, "nuclei_counts_per_sample.tsv"), "\n")

cat("\nDone. Plots in:", plots_dir, "\n")
cat("       TSVs in:",  tsv_dir,   "\n")

# ---- Free memory ----
rm(list = ls())
gc(full = TRUE, reset = TRUE)
