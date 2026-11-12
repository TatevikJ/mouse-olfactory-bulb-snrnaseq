# Post-filter QC diagnostics (v4). Reads the v4 QC-filtered RDS (basic-QC
# v2 with lenient doublet 0.3, plus brain-region exclusion from the
# TO_REMOVE list) and writes violins, scatters, UMAPs, bootstrap boxplots,
# and summary TSVs to outputs/qc_post_v4/{plots,tsv}/. Mirrors qc_post_v2.R.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
  library(scDblFinder)
  library(SingleCellExperiment)
})

# ---- Filter thresholds ----
thresh_scDblFinder_max <- 0.3
thresh_percent_mt_max  <- 0.5
thresh_nCount_RNA_max  <- 400000
# nFeature_RNA min filter dropped: Parse DGE_filtered min >= 260 across samples.

# ---- Inputs / outputs ----
in_rds    <- "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v4.rds"
plots_dir <- "outputs/qc_post_v4/plots"
tsv_dir   <- "outputs/qc_post_v4/tsv"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsv_dir,   recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
cat("Loading:", in_rds, "\n")
dataSeurat <- readRDS(in_rds)
Idents(dataSeurat) <- dataSeurat@meta.data$sample
if (!"condition" %in% colnames(dataSeurat@meta.data)) {
  dataSeurat$condition <- ifelse(grepl("^ASO", dataSeurat$sample), "ASO", "WT")
}
cat("Nuclei:", ncol(dataSeurat), " Genes:", nrow(dataSeurat), "\n")

# Condition palette: WT = teal, ASO = red
condition_colors <- c(WT = "#59b2aa", ASO = "#d43f2e")
# Per-sample palette: levels(factor(...)) is alphabetical, so ASO1..ASO6 come first
cnames <- setNames(
  rep(c(condition_colors[["ASO"]], condition_colors[["WT"]]), each = 6),
  levels(factor(dataSeurat$sample))
)

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

# QC metrics should already be present from qc_pre.R; recompute defensively if not.
needs_save <- FALSE

if (!"percent.mt" %in% colnames(dataSeurat@meta.data)) {
  cat("Computing percent.mt...\n")
  mt_genes <- grep("^mt-", rownames(dataSeurat), ignore.case = TRUE, value = TRUE)
  dataSeurat$percent.mt <- if (length(mt_genes) > 0) {
    PercentageFeatureSet(dataSeurat, features = mt_genes)
  } else 0
  needs_save <- TRUE
}

if (!"scDblFinder.score" %in% colnames(dataSeurat@meta.data)) {
  cat("Running scDblFinder (expensive: copies to SCE, simulates doublets)...\n")
  sce <- as.SingleCellExperiment(dataSeurat)
  sce <- scDblFinder(sce, samples = dataSeurat$sample)
  dataSeurat$scDblFinder.score <- sce$scDblFinder.score
  dataSeurat$scDblFinder.class <- sce$scDblFinder.class
  rm(sce); invisible(gc())
  needs_save <- TRUE
}

if (needs_save) {
  cat("Persisting QC metrics back to:", in_rds, "\n")
  saveRDS(dataSeurat, in_rds)
}

# Print available QC ranges so threshold choices are anchored to data
cat("\nQC metric ranges:\n")
cat("  nCount_RNA:        "); print(summary(dataSeurat$nCount_RNA))
cat("  nFeature_RNA:      "); print(summary(dataSeurat$nFeature_RNA))
cat("  percent.mt:        "); print(summary(dataSeurat$percent.mt))
cat("  scDblFinder.score: "); print(summary(dataSeurat$scDblFinder.score))

# ============================================================
# Per-sample violin plots
# ============================================================
cat("\n[1/4] Per-sample violin plots...\n")

format_thresh <- function(x) format(x, big.mark = ",", scientific = FALSE,
                                    trim = TRUE)

make_vln <- function(feature, thresh = NULL) {
  p <- VlnPlot(dataSeurat, features = feature, layer = "counts",
               group.by = "sample", raster = FALSE,
               alpha = 0.9, pt.size = 0) +
    geom_jitter(width = 0.3, size = 0.4, alpha = 0.18,
                color = "gray15", shape = 16) +
    scale_fill_manual(values = cnames) +
    labs(x = "Sample", y = feature) +
    project_base_theme() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  if (!is.null(thresh)) {
    p <- p +
      geom_hline(yintercept = thresh, linetype = "dashed",
                 color = "black", alpha = 0.8, linewidth = 0.5) +
      annotate("text", x = Inf, y = thresh,
               label = format_thresh(thresh),
               hjust = 1.05, vjust = -0.3, size = 3.2, color = "black")
  }
  p
}

p_nCount_RNA    <- make_vln("nCount_RNA",        thresh = thresh_nCount_RNA_max)
p_nFeature_RNA  <- make_vln("nFeature_RNA")
p_percent_mt    <- make_vln("percent.mt",        thresh = thresh_percent_mt_max)
p_doubletscore  <- make_vln("scDblFinder.score",
                            thresh = thresh_scDblFinder_max)

ggsave(file.path(plots_dir, "nCount_RNA_violin.pdf"),    plot = p_nCount_RNA,   width = 7, height = 5)
ggsave(file.path(plots_dir, "nFeature_RNA_violin.pdf"),  plot = p_nFeature_RNA, width = 7, height = 5)
ggsave(file.path(plots_dir, "percent_mt_violin.pdf"),    plot = p_percent_mt,   width = 7, height = 5)
ggsave(file.path(plots_dir, "doublet_score_violin.pdf"), plot = p_doubletscore, width = 7, height = 5)

# ============================================================
# Per-sample nuclei-count bar plot (filled by condition)
# ============================================================
cat("[1b/4] Per-sample nuclei-count bar plot...\n")

counts_df <- dataSeurat@meta.data %>%
  dplyr::count(sample, condition, name = "n_nuclei") %>%
  dplyr::arrange(condition, sample)

p_nucleicounts <- ggplot(counts_df,
                         aes(x = sample, y = n_nuclei, fill = condition)) +
  geom_col(width = 0.9, alpha = 1, color = "gray30", linewidth = 0.3) +
  geom_text(aes(label = format(n_nuclei, big.mark = ",")),
            vjust = 1.5, size = 2.8, color = "black") +
  scale_fill_manual(values = condition_colors) +
  labs(x = "Sample", y = "Nuclei (n)") +
  project_base_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plots_dir, "nuclei_counts_per_sample.pdf"),
       plot = p_nucleicounts, width = 7, height = 5)

# ============================================================
# Per-sample doublet-score histograms (one figure, faceted)
# ============================================================
cat("[1c/4] Per-sample doublet-score histograms...\n")

doublet_thresh_df <- data.frame(
  xint  = thresh_scDblFinder_max,
  label = format_thresh(thresh_scDblFinder_max)
)

p_doublet_hist <- ggplot(dataSeurat@meta.data,
                         aes(x = scDblFinder.score, fill = condition)) +
  geom_histogram(bins = 60, color = "gray30", linewidth = 0.15) +
  geom_vline(data = doublet_thresh_df,
             aes(xintercept = xint),
             linetype = "dashed", color = "black",
             alpha = 0.8, linewidth = 0.5,
             inherit.aes = FALSE) +
  geom_text(data = doublet_thresh_df,
            aes(x = xint, y = Inf, label = label),
            hjust = -0.15, vjust = 1.4, size = 3, color = "black",
            inherit.aes = FALSE) +
  facet_wrap(~ sample, ncol = 4) +
  scale_fill_manual(values = condition_colors) +
  labs(x = "scDblFinder.score", y = "Nuclei") +
  project_base_theme() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "gray92", color = NA),
        strip.text = element_text(size = 11))
ggsave(file.path(plots_dir, "doublet_score_histogram_per_sample.pdf"),
       plot = p_doublet_hist, width = 10, height = 7)

# ============================================================
# Feature scatter plots
# ============================================================
cat("[2/4] Feature scatter plots...\n")

p_percent_mt_nFeature_RNA <- FeatureScatter(
  dataSeurat, feature1 = "percent.mt", feature2 = "nFeature_RNA",
  group.by = "sample", raster = FALSE
) +
  geom_vline(xintercept = thresh_percent_mt_max,   linetype = "dashed",
             color = "black", alpha = 0.8, linewidth = 0.5) +
  annotate("text", x = thresh_percent_mt_max,   y = Inf,
           label = format_thresh(thresh_percent_mt_max),
           hjust = -0.15, vjust = 1.4, size = 3.2, color = "black") +
  project_base_theme() +
  theme(plot.title = element_blank())

p_nCount_RNA_nFeature_RNA <- FeatureScatter(
  dataSeurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
  group.by = "sample", raster = TRUE
) +
  geom_vline(xintercept = thresh_nCount_RNA_max,   linetype = "dashed",
             color = "black", alpha = 0.8, linewidth = 0.5) +
  annotate("text", x = thresh_nCount_RNA_max,   y = Inf,
           label = format_thresh(thresh_nCount_RNA_max),
           hjust = -0.15, vjust = 1.4, size = 3.2, color = "black") +
  project_base_theme() +
  theme(plot.title = element_blank())

ggsave(file.path(plots_dir, "percent_mt_vs_nFeature_scatter.pdf"),
       plot = p_percent_mt_nFeature_RNA, width = 7, height = 6)
ggsave(file.path(plots_dir, "nCount_vs_nFeature_scatter.pdf"),
       plot = p_nCount_RNA_nFeature_RNA, width = 7, height = 6)

# ============================================================
# UMAP plots: QC metrics + annotations + sample/condition
# ============================================================
cat("[3/4] UMAP feature plots...\n")

# Logs each filename before ggsave so patchwork failures pinpoint the file.
save_pdf <- function(path, plot, width, height) {
  cat("  saving:", basename(path), sprintf("(%g x %g)\n", width, height))
  ggsave(path, plot = plot, width = width, height = height)
}

make_feature_umap <- function(feature, title = feature) {
  FeaturePlot(dataSeurat, features = feature, reduction = "umap") +
    scale_color_viridis_c(option = "plasma") +
    labs(x = "UMAP1", y = "UMAP2", title = title) +
    project_umap_theme("top-left")
}

p_nCount_umap        <- make_feature_umap("nCount_RNA")
p_nFeature_umap      <- make_feature_umap("nFeature_RNA")
p_percent_mt_umap    <- make_feature_umap("percent.mt")
p_doubletscore_umap  <- make_feature_umap("scDblFinder.score")

save_pdf(file.path(plots_dir, "umap_nCount_RNA.pdf"),    plot = p_nCount_umap,       width = 6, height = 6)
save_pdf(file.path(plots_dir, "umap_nFeature_RNA.pdf"),  plot = p_nFeature_umap,     width = 6, height = 6)
save_pdf(file.path(plots_dir, "umap_percent_mt.pdf"),    plot = p_percent_mt_umap,   width = 6, height = 6)
save_pdf(file.path(plots_dir, "umap_doublet_score.pdf"), plot = p_doubletscore_umap, width = 6, height = 6)

p_umap_sample <- DimPlot(dataSeurat, reduction = "umap",
                         group.by = "sample") +
  labs(x = "UMAP1", y = "UMAP2", title = "Sample") +
  project_umap_theme("top-left") +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )
p_umap_condition <- DimPlot(dataSeurat, reduction = "umap",
                            group.by = "condition") +
  scale_color_manual(values = condition_colors) +
  labs(x = "UMAP1", y = "UMAP2", title = "Condition") +
  project_umap_theme("top-left") +
  theme(
    legend.position = "right",
    legend.justification = "center",
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )
save_pdf(file.path(plots_dir, "umap_sample.pdf"),    plot = p_umap_sample,    width = 8, height = 6)
save_pdf(file.path(plots_dir, "umap_condition.pdf"), plot = p_umap_condition, width = 7, height = 6)

# Top-15 + "other" collapse for class_name / subclass_name.
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

class_cols    <- intersect(c("class_name",    "class_new"),
                           colnames(dataSeurat@meta.data))
subclass_cols <- intersect(c("subclass_name", "subclass_new"),
                           colnames(dataSeurat@meta.data))

# Cast to factor so DimPlot picks a discrete palette; place "other" last.
to_ordered_factor <- function(x) {
  x <- as.character(x)
  lvls <- sort(setdiff(unique(x), c(NA, "other")))
  if ("other" %in% x) lvls <- c(lvls, "other")
  factor(x, levels = lvls)
}
for (col in c(class_cols, subclass_cols)) {
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

for (col in class_cols) {
  is_full <- identical(col, "class_name")
  legend_text_theme <- if (is_full) {
    theme(legend.position = "right",
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
  } else {
    theme(legend.position = "right",
          legend.background = element_blank(),
          legend.box.background = element_blank())
  }
  # Wider canvas + 2-col legend for full taxonomy (else patchwork errors).
  legend_width  <- if (is_full) 12 else 8
  legend_height <- if (is_full) 8  else 6
  legend_guide  <- if (is_full) {
    guides(color = guide_legend(ncol = 2, byrow = TRUE,
                                override.aes = list(size = 2)))
  } else NULL

  p_no_legend <- DimPlot(dataSeurat, reduction = "umap",
                         group.by = col, label = FALSE) +
    labs(x = "UMAP1", y = "UMAP2", title = paste0("Class (", col, ")")) +
    project_umap_theme("top-left") +
    NoLegend()
  p_no_legend <- add_centroid_labels(p_no_legend, dataSeurat@meta.data[[col]])
  save_pdf(file.path(plots_dir, paste0("umap_", col, ".pdf")),
           plot = p_no_legend, width = 6, height = 6)

  p_with_legend <- DimPlot(dataSeurat, reduction = "umap",
                           group.by = col, label = FALSE) +
    labs(x = "UMAP1", y = "UMAP2", title = paste0("Class (", col, ")")) +
    project_umap_theme("top-left") +
    legend_text_theme +
    legend_guide
  p_with_legend <- add_centroid_labels(p_with_legend, dataSeurat@meta.data[[col]])
  save_pdf(file.path(plots_dir, paste0("umap_", col, "_legend.pdf")),
           plot = p_with_legend, width = legend_width, height = legend_height)
}

for (col in subclass_cols) {
  is_collapsed <- identical(col, "subclass_new")
  # Multi-col legend + bigger canvas for full taxonomy (~100 levels).
  legend_width     <- if (is_collapsed) 8   else 14
  legend_height    <- if (is_collapsed) 6   else 9
  legend_text_size <- if (is_collapsed) 8   else 6
  legend_key_size  <- if (is_collapsed) 0.4 else 0.32
  legend_ncol      <- if (is_collapsed) 1   else 3

  p_no_legend <- DimPlot(dataSeurat, reduction = "umap",
                         group.by = col, label = TRUE, repel = TRUE,
                         label.size = 2.5) +
    labs(x = "UMAP1", y = "UMAP2", title = paste0("Subclass (", col, ")")) +
    project_umap_theme("top-left") +
    NoLegend()
  save_pdf(file.path(plots_dir, paste0("umap_", col, ".pdf")),
           plot = p_no_legend, width = 6, height = 6)

  p_with_legend <- DimPlot(dataSeurat, reduction = "umap",
                           group.by = col, label = TRUE, repel = TRUE,
                           label.size = 2.5) +
    labs(x = "UMAP1", y = "UMAP2", title = paste0("Subclass (", col, ")")) +
    project_umap_theme("top-left") +
    guides(color = guide_legend(ncol = legend_ncol, byrow = TRUE,
                                override.aes = list(size = 2))) +
    theme(legend.position = "right",
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.text = element_text(size = legend_text_size),
          legend.key.size = unit(legend_key_size, "cm"))
  save_pdf(file.path(plots_dir, paste0("umap_", col, "_legend.pdf")),
           plot = p_with_legend, width = legend_width, height = legend_height)
}

# Bootstrap groups by full taxonomy ("other" would hide low-confidence subclasses).
class_col    <- if ("class_name"    %in% colnames(dataSeurat@meta.data)) "class_name"    else NA_character_
subclass_col <- if ("subclass_name" %in% colnames(dataSeurat@meta.data)) "subclass_name" else NA_character_

# ============================================================
# Bootstrapping probability diagnostics
# ============================================================
cat("[4/4] Bootstrapping probability diagnostics...\n")

bootstrap_cols <- c("class_bootstrapping_probability",
                    "subclass_bootstrapping_probability")
bootstrap_present <- bootstrap_cols[bootstrap_cols %in% colnames(dataSeurat@meta.data)]

bootstrap_group_col <- c(
  class    = class_col,
  subclass = subclass_col
)

for (col in bootstrap_present) {
  level_name <- sub("_bootstrapping_probability$", "", col)
  pretty     <- paste0(toupper(substring(level_name, 1, 1)), substring(level_name, 2))
  group_col  <- bootstrap_group_col[[level_name]]

  if (!is.null(group_col) && group_col %in% colnames(dataSeurat@meta.data)) {
    is_subclass <- identical(level_name, "subclass")
    x_text_size <- if (is_subclass) 6 else 9
    x_text_angle <- if (is_subclass) 90 else 45
    x_text_vjust <- if (is_subclass) 0.5 else 1
    plot_width  <- if (is_subclass) 16 else 10

    md <- dataSeurat@meta.data
    md <- md[!is.na(md[[group_col]]) & !is.na(md[[col]]), , drop = FALSE]

    p_box <- ggplot(md, aes(x = reorder(.data[[group_col]], .data[[col]], median),
                            y = .data[[col]])) +
      geom_boxplot(fill = "#4A90E2", alpha = 0.9, outlier.size = 0.1,
                   color = "gray30", linewidth = 0.4) +
      geom_hline(yintercept = c(0.5, 0.7), linetype = "dashed", alpha = 0.5) +
      labs(x = group_col,
           y = paste(pretty, "Bootstrapping Probability")) +
      project_base_theme() +
      theme(
        axis.text.x = element_text(angle = x_text_angle, hjust = 1,
                                   vjust = x_text_vjust, size = x_text_size),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
      )
    ggsave(file.path(plots_dir,
                     paste0("boxplot_", level_name, "_bootstrap_by_", group_col, ".pdf")),
           plot = p_box, width = plot_width, height = 6)
  } else {
    cat("  Skipping ", level_name, " per-group boxplot: column '",
        group_col, "' not in metadata.\n", sep = "")
  }

  p_umap_boot <- FeaturePlot(dataSeurat, features = col, reduction = "umap") +
    scale_color_viridis_c(option = "plasma") +
    labs(x = "UMAP1", y = "UMAP2",
         title = paste(pretty, "Bootstrapping Probability")) +
    project_umap_theme("top-left")
  ggsave(file.path(plots_dir, paste0("umap_", level_name, "_bootstrap_prob.pdf")),
         plot = p_umap_boot, width = 6, height = 6)
}

# Supertype-level bootstrap UMAP only (no boxplot: too many supertypes for that view).
if ("supertype_bootstrapping_probability" %in% colnames(dataSeurat@meta.data)) {
  p_umap_super_boot <- FeaturePlot(dataSeurat,
                                   features = "supertype_bootstrapping_probability",
                                   reduction = "umap") +
    scale_color_viridis_c(option = "plasma") +
    labs(x = "UMAP1", y = "UMAP2",
         title = "Supertype Bootstrapping Probability") +
    project_umap_theme("top-left")
  ggsave(file.path(plots_dir, "umap_supertype_bootstrap_prob.pdf"),
         plot = p_umap_super_boot, width = 6, height = 6)
}

# ============================================================
# Summary TSVs
# ============================================================
cat("Writing summary TSVs...\n")

qc_summary_per_sample <- dataSeurat@meta.data %>%
  group_by(sample, condition) %>%
  summarise(
    n_nuclei              = n(),
    min_nCount_RNA        = min(nCount_RNA),
    median_nCount_RNA     = median(nCount_RNA),
    max_nCount_RNA        = max(nCount_RNA),
    min_nFeature_RNA      = min(nFeature_RNA),
    median_nFeature_RNA   = median(nFeature_RNA),
    max_nFeature_RNA      = max(nFeature_RNA),
    min_percent_mt        = min(percent.mt),
    median_percent_mt     = median(percent.mt),
    max_percent_mt        = max(percent.mt),
    n_nuclei_doublet_ge_0.3    = sum(scDblFinder.score >= thresh_scDblFinder_max),
    n_nuclei_percent_mt_ge_0.5 = sum(percent.mt >= thresh_percent_mt_max),
    n_nuclei_nCount_ge_400k    = sum(nCount_RNA >= thresh_nCount_RNA_max),
    .groups = "drop"
  ) %>%
  arrange(condition, sample)
write.table(qc_summary_per_sample,
            file.path(tsv_dir, "qc_metrics_summary_per_sample.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  ", file.path(tsv_dir, "qc_metrics_summary_per_sample.tsv"), "\n")

# Per-(condition, sample, class, subclass, supertype) nuclei counts
nuclei_count_cols <- intersect(
  c("condition", "sample", "class_name", "subclass_name", "supertype_name"),
  colnames(dataSeurat@meta.data)
)
nuclei_counts <- dataSeurat@meta.data %>%
  dplyr::count(across(all_of(nuclei_count_cols)), name = "n_nuclei") %>%
  arrange(condition, sample)
write.table(nuclei_counts,
            file.path(tsv_dir, "nuclei_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  ", file.path(tsv_dir, "nuclei_counts.tsv"), "\n")

# ---- Mapping-stats TSVs (kept nuclei only; per Allen-taxonomy level) ----
# Operates on the post-filter object, so n_nuclei is the kept count and the
# bootstrap aggregates describe the mapping confidence of the retained nuclei.
mapping_probability_stats_class <- dataSeurat@meta.data %>%
  group_by(class_name) %>%
  summarise(
    n_nuclei                               = n(),
    mean_class_bootstrapping_probability   = mean(class_bootstrapping_probability,   na.rm = TRUE),
    median_class_bootstrapping_probability = median(class_bootstrapping_probability, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(class_name)
write.table(mapping_probability_stats_class,
            file.path(tsv_dir, "mapping_probability_stats_class.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  ", file.path(tsv_dir, "mapping_probability_stats_class.tsv"), "\n")

mapping_probability_stats_subclass <- dataSeurat@meta.data %>%
  group_by(class_name, subclass_name) %>%
  summarise(
    n_nuclei                                  = n(),
    mean_subclass_bootstrapping_probability   = mean(subclass_bootstrapping_probability,   na.rm = TRUE),
    median_subclass_bootstrapping_probability = median(subclass_bootstrapping_probability, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(class_name, subclass_name)
write.table(mapping_probability_stats_subclass,
            file.path(tsv_dir, "mapping_probability_stats_subclass.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  ", file.path(tsv_dir, "mapping_probability_stats_subclass.tsv"), "\n")

mapping_probability_stats_supertype <- dataSeurat@meta.data %>%
  group_by(class_name, subclass_name, supertype_name) %>%
  summarise(
    n_nuclei                                   = n(),
    mean_supertype_bootstrapping_probability   = mean(supertype_bootstrapping_probability,   na.rm = TRUE),
    median_supertype_bootstrapping_probability = median(supertype_bootstrapping_probability, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(class_name, subclass_name, supertype_name)
write.table(mapping_probability_stats_supertype,
            file.path(tsv_dir, "mapping_probability_stats_supertype.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  ", file.path(tsv_dir, "mapping_probability_stats_supertype.tsv"), "\n")

cat("\nDone. Plots in:", plots_dir, "\n")
cat("       TSVs in:",  tsv_dir,   "\n")

# ---- Free memory ----
rm(list = ls())
gc(full = TRUE, reset = TRUE)
