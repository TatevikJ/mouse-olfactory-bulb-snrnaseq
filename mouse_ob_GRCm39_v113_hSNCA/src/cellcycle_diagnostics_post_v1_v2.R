# Diagnostic: are high-nFeature_RNA nuclei enriched in cycling (S/G2M) phases?
# Run CellCycleScoring on the filtered (normalized) object, flag high_nFeature
# (nFeature_RNA > 8000), and for each subclass cross-tabulate Phase against
# high_nFeature with both raw counts and within-group proportions. Run for both
# qc_post_v1 and qc_post_v2. Mature-neuron subclasses (041, 035, ...) scoring
# G2M among high_nFeature nuclei would suggest missed doublets; G1 supports
# biologically complex nuclei.
#
# Tirosh cell-cycle genes in title-case mouse symbols (Mcm5, Pcna, ...), matching
# GRCm39 var_names. Pre-converted list:
#   https://github.com/David-R-bird/Seurat_ccgenes/blob/main/mouse_cell_cycle_genes_2019.rds
# Place it at meta/mouse_cell_cycle_genes_2019.rds before running.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

HIGH_NFEATURE_THRESH <- 8000
phase_levels <- c("G1", "S", "G2M")

# Shared UMAP theme: blank square panel, so the UMAP is the same size with or
# without an outside legend.
umap_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text  = element_blank(),
    aspect.ratio = 1,
    axis.title.y = element_text(margin = margin(r = 8))
  )

# Score UMAPs: colorbar inside, top-left. The Phase UMAP keeps its legend outside.
legend_inside <- theme(
  legend.position = c(0.02, 0.98),
  legend.justification = c(0, 1),
  legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
  legend.title = element_blank()
)

# Top-15 + "other" collapse for the MapMyCells label columns (as in
# qc_post_v1.R); used to build the collapsed _new columns for a legible dotplot.
collapse_to_top_n <- function(x, n = 15, other_label = "other") {
  x <- as.character(x)
  freq <- sort(table(x[!is.na(x)]), decreasing = TRUE)
  top <- names(freq)[seq_len(min(n, length(freq)))]
  ifelse(is.na(x) | x %in% top, x, other_label)
}
to_ordered_factor <- function(x) {
  x <- as.character(x)
  lvls <- sort(setdiff(unique(x), c(NA, "other")))
  if ("other" %in% x) lvls <- c(lvls, "other")
  factor(x, levels = lvls)
}

cc <- readRDS("meta/mouse_cell_cycle_genes_2019.rds")
# Drop "N/A" placeholders: human cc genes with no mouse ortholog in the
# converted list (2 in s.genes). AddModuleScore would drop them anyway, this
# just silences the "features not present: N/A" warning.
s_genes   <- cc$s.genes[cc$s.genes != "N/A" & !is.na(cc$s.genes)]
g2m_genes <- cc$m.genes[cc$m.genes != "N/A" & !is.na(cc$m.genes)]
stopifnot(length(s_genes) > 0, length(g2m_genes) > 0)

pipelines <- list(
  v1 = list(
    rds      = "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v1.rds",
    tsv_dir  = "outputs/qc_post_v1/tsv",
    plot_dir = "outputs/qc_post_v1/plots/cellcycle_diagnostics"
  ),
  v2 = list(
    rds      = "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v2.rds",
    tsv_dir  = "outputs/qc_post_v2/tsv",
    plot_dir = "outputs/qc_post_v2/plots/cellcycle_diagnostics"
  )
)

run_pipeline <- function(tag, paths) {
  cat("\n==== Pipeline:", tag, "====\n")
  stopifnot(file.exists(paths$rds))
  dir.create(paths$tsv_dir, recursive = TRUE, showWarnings = FALSE)

  cat("Loading:", paths$rds, "\n")
  obj <- readRDS(paths$rds)
  DefaultAssay(obj) <- "RNA"

  obj <- CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes,
                          set.ident = FALSE)
  obj$Phase <- factor(obj$Phase, levels = phase_levels)

  # UMAPs colored by the cell-cycle scores (continuous) and Phase
  # (discrete). UMAP carried from normalize_pre.R (subset, not recomputed).
  dir.create(paths$plot_dir, recursive = TRUE, showWarnings = FALSE)
  for (metric in c("S.Score", "G2M.Score")) {
    # Scores are signed and centered near 0: diverging palette, grey at 0,
    # blue below, red above; nonzero nuclei drawn on top.
    p <- FeaturePlot(obj, features = metric, reduction = "umap", order = TRUE) +
      scale_color_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                            midpoint = 0) +
      labs(x = "UMAP1", y = "UMAP2", title = metric) + umap_theme + legend_inside
    ggsave(file.path(paths$plot_dir, paste0("umap_", metric, ".pdf")),
           p, width = 6, height = 6)
  }
  p <- DimPlot(obj, group.by = "Phase", reduction = "umap") +
    labs(x = "UMAP1", y = "UMAP2", title = "Phase") + umap_theme
  ggsave(file.path(paths$plot_dir, "umap_Phase.pdf"), p, width = 7.5, height = 6)
  cat("wrote UMAPs to:", paths$plot_dir, "\n")

  # Dot plots of cell-cycle genes by subclass (full and top-15 collapsed) and by
  # predicted Phase. Color = mean expression (unscaled), size = % detected.
  for (lvl in c("class", "subclass", "supertype")) {
    name_col <- paste0(lvl, "_name"); new_col <- paste0(lvl, "_new")
    if (name_col %in% colnames(obj@meta.data) &&
        !new_col %in% colnames(obj@meta.data)) {
      obj@meta.data[[new_col]] <-
        to_ordered_factor(collapse_to_top_n(obj@meta.data[[name_col]], n = 15))
    }
  }

  cc_feat <- list("S genes"   = intersect(s_genes,   rownames(obj)),
                  "G2M genes" = intersect(g2m_genes, rownames(obj)))
  n_genes <- length(unlist(cc_feat))
  for (grp_col in c("subclass_name", "subclass_new", "Phase")) {
    p <- DotPlot(obj, features = cc_feat, group.by = grp_col,
                 scale = FALSE, dot.scale = 5) +
      scale_color_viridis_c(option = "plasma") +
      labs(x = NULL, y = NULL, color = "Mean expr", size = "% detected",
           title = paste("Cell-cycle gene expression by", grp_col)) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 8),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", color = NA)
      )
    n_grp <- length(unique(as.character(obj@meta.data[[grp_col]])))
    ggsave(file.path(paths$plot_dir, paste0("dotplot_cellcycle_gene_by_", grp_col, ".pdf")),
           p, width = max(12, n_genes * 0.18), height = max(4, n_grp * 0.25))
  }
  cat("wrote cc-gene dotplots to:", paths$plot_dir, "\n")

  md <- obj@meta.data
  rm(obj); invisible(gc())
  md$Phase <- factor(md$Phase, levels = phase_levels)
  md$grp <- factor(ifelse(md$nFeature_RNA > HIGH_NFEATURE_THRESH,
                          "high_nFeature_RNA", "low_nFeature_RNA"),
                   levels = c("high_nFeature_RNA", "low_nFeature_RNA"))

  # Per-subclass Phase x high_nFeature: counts + within-group proportions.
  tab <- md %>%
    group_by(subclass_name, grp, Phase, .drop = FALSE) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(subclass_name, grp) %>%
    mutate(prop = if (sum(n) > 0) n / sum(n) else 0) %>%
    ungroup() %>%
    select(subclass_name, Phase, grp, n, prop) %>%
    pivot_wider(names_from = grp, values_from = c(n, prop), values_fill = 0) %>%
    arrange(subclass_name, Phase)
  out_tab <- file.path(paths$tsv_dir, "cellcycle_phase_by_high_nfeature_per_subclass.tsv")
  write.table(tab, out_tab, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("wrote:", out_tab, "\n")

  rm(md); invisible(gc())
}

for (tag in names(pipelines)) run_pipeline(tag, pipelines[[tag]])

cat("\nDone.\n")
