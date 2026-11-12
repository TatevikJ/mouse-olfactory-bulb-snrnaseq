# Over-representation analysis (ORA, GO:BP) of the ASO vs WT DEG results
# from deg_post_v4_normalized.R. clusterProfiler::enrichGO per cell type on
# the up- and down-regulated DEGs separately.
#
#   Tracks:   pseudobulk (DEG: FDR < 0.05 & |logFC| > 0.5; primary) and
#             singlecell_mast (DEG: Bonferroni < 0.01 & |coef| > 0.5;
#             nucleus-level, supplementary).
#   Variants: all_genes, and no_Pcdh11x (Pcdh11x removed from DEGs and universe).
#   Levels:   class_name, subclass_name, supertype_name.
#
# DEGs are the `significant` column of the DEG tables. Universe = genes tested
# in that cell type (non-NA adjusted p-value), not the whole genome. Hypergeometric
# test; GO terms of 10-500 genes; BH across terms within each cell type and
# direction. All terms are saved (pvalueCutoff = qvalueCutoff = 1);
# significance is p.adjust < 0.05.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(dplyr)
  library(ggplot2)
})

# ---- Inputs / outputs ----
in_root  <- "outputs/deg_post_v4_normalized"
out_root <- "outputs/ora_post_v4_normalized"
# Track -> DEG file prefix, effect-size column (direction), adjusted p column (universe).
tracks <- list(
  pseudobulk      = list(prefix = "pseudobulk", fc_col = "logFC", p_col = "FDR"),
  singlecell_mast = list(prefix = "mast",       fc_col = "coef",  p_col = "padjust")
)
# Variant -> genes removed before analysis.
variants <- list(all_genes = character(0), no_Pcdh11x = "Pcdh11x")
# Each (track, variant) writes to its own subtree with rds/, tsv/, plots/.
track_dirs <- function(track, variant) {
  d <- list(rds   = file.path(out_root, track, variant, "rds"),
            tsv   = file.path(out_root, track, variant, "tsv"),
            plots = file.path(out_root, track, variant, "plots"))
  for (x in d) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  d
}

# ---- ORA params ----
levels_to_run   <- c("class_name", "subclass_name", "supertype_name")
min_gs_size     <- 10    # GO term size bounds (clusterProfiler defaults)
max_gs_size     <- 500
min_degs        <- 5     # skip (cell type, direction) with fewer DEGs
padj_cutoff     <- 0.05  # significant term: BH p.adjust < 0.05
n_dot_terms     <- 20    # dotplot: top terms
n_heatmap_terms <- 5     # heatmap: top terms per cell type and direction

sanitize <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)

# ============================================================
# ORA per cell type and direction
# ============================================================
run_ora_level <- function(deg_list, fc_col, p_col, drop_genes) {
  out <- list()
  for (ct in names(deg_list)) {
    df <- as.data.frame(deg_list[[ct]])
    df <- df[!is.na(df[[p_col]]) & !df$gene %in% drop_genes, ]
    universe <- df$gene
    deg_sets <- list(
      up   = df$gene[df$significant %in% TRUE & df[[fc_col]] > 0],
      down = df$gene[df$significant %in% TRUE & df[[fc_col]] < 0]
    )
    for (dir in names(deg_sets)) {
      genes <- deg_sets[[dir]]
      if (length(genes) < min_degs) {
        cat("  Skipping '", ct, "' ", dir, " (", length(genes), " DEGs)\n", sep = "")
        next
      }
      ego <- tryCatch({
        clusterProfiler::enrichGO(
          gene          = genes,
          universe      = universe,
          OrgDb         = org.Mm.eg.db,
          keyType       = "SYMBOL",
          ont           = "BP",
          minGSSize     = min_gs_size,
          maxGSSize     = max_gs_size,
          pAdjustMethod = "BH",
          pvalueCutoff  = 1,
          qvalueCutoff  = 1
        )
      }, error = function(e) {
        cat("  enrichGO failed for '", ct, "' ", dir, ": ", conditionMessage(e), "\n", sep = "")
        NULL
      })
      if (is.null(ego)) next
      ego@geneSets <- ego@geneSets[ego@result$ID]  # keep tested terms only (RDS size)
      out[[ct]][[dir]] <- ego
      cat(sprintf("  %s %s: %d DEGs / %d universe, %d terms, %d significant\n",
                  ct, dir, length(genes), length(universe),
                  nrow(ego@result), sum(ego@result$p.adjust < padj_cutoff)))
    }
  }
  out
}

# ============================================================
# Plots
# ============================================================
# Dotplot of significant terms for one (cell type, direction).
plot_dot <- function(ego, ct, dir, title_suffix, out_dir) {
  sig <- ego
  sig@result <- ego@result[ego@result$p.adjust < padj_cutoff, ]
  if (nrow(sig@result) == 0) return(invisible(NULL))
  p <- enrichplot::dotplot(sig, showCategory = n_dot_terms) +
    ggtitle(sprintf("%s, %s in ASO (GO:BP ORA, %s)", ct, dir, title_suffix)) +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(out_dir, sprintf("%s__%s__dotplot.pdf", sanitize(ct), dir)), p, width = 10, height = 8)
}

# Heatmap across cell types: -log10 p.adjust of top significant terms, faceted by direction.
plot_padj_heatmap <- function(tsv_df, level_col, title_suffix, out_pdf) {
  sig <- tsv_df[tsv_df$p.adjust < padj_cutoff, ]
  if (nrow(sig) == 0) return(invisible(NULL))
  terms <- sig %>%
    group_by(cell_type, direction) %>%
    slice_min(p.adjust, n = n_heatmap_terms, with_ties = FALSE) %>%
    pull(Description) %>%
    unique()
  hm <- sig %>%
    filter(Description %in% terms) %>%
    mutate(neg_log10_padj = -log10(p.adjust),
           direction = factor(direction, levels = c("up", "down")))
  p <- ggplot(hm, aes(cell_type, Description, fill = neg_log10_padj)) +
    geom_tile(color = "white") +
    facet_grid(. ~ direction, scales = "free_x", space = "free_x") +
    scale_fill_viridis_c(option = "plasma", name = "-log10\np.adjust") +
    labs(x = NULL, y = NULL,
         title = sprintf("GO:BP ORA by %s (%s), significant terms only", level_col, title_suffix)) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(),
          plot.title.position = "plot",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(out_pdf, p, limitsize = FALSE,
         width  = max(8, 4 + 0.25 * nrow(distinct(hm, cell_type, direction))),
         height = max(6, 2 + 0.15 * length(terms)))
  cat("wrote:", out_pdf, "\n")
}

# Save per-level RDS (cell type -> list(up, down) of enrichResult), combined TSV, plots.
save_level <- function(res_list, level_col, dirs, title_suffix) {
  if (length(res_list) == 0) return(invisible(NULL))
  out <- file.path(dirs$rds, sprintf("ora_results_%s.rds", level_col))
  saveRDS(res_list, out); cat("wrote:", out, "\n")

  tsv_df <- dplyr::bind_rows(lapply(names(res_list), function(ct) {
    dplyr::bind_rows(lapply(names(res_list[[ct]]), function(dir) {
      df <- res_list[[ct]][[dir]]@result
      if (nrow(df) == 0) return(NULL)
      df$cell_type <- ct
      df$direction <- dir
      df
    }))
  }))
  out_tsv <- file.path(dirs$tsv, sprintf("ora_%s.tsv", level_col))
  write.table(tsv_df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("wrote:", out_tsv, "\n")

  plot_dir <- file.path(dirs$plots, level_col)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  for (ct in names(res_list)) {
    for (dir in names(res_list[[ct]])) {
      plot_dot(res_list[[ct]][[dir]], ct, dir, title_suffix, plot_dir)
    }
  }
  plot_padj_heatmap(tsv_df, level_col, title_suffix, file.path(plot_dir, "padj_heatmap.pdf"))
}

# ============================================================
# Run all tracks x levels x variants
# ============================================================
for (tr in names(tracks)) {
  for (lvl in levels_to_run) {
    in_rds <- file.path(in_root, tr, "rds",
                        sprintf("de_results_%s_%s.rds", tracks[[tr]]$prefix, lvl))
    cat("\nLoading:", in_rds, "\n")
    deg_list <- readRDS(in_rds)
    for (vr in names(variants)) {
      cat(sprintf("\n--- ORA: %s / %s / %s ---\n", tr, vr, lvl))
      res <- run_ora_level(deg_list, tracks[[tr]]$fc_col, tracks[[tr]]$p_col, variants[[vr]])
      save_level(res, lvl, track_dirs(tr, vr), sprintf("%s, %s", tr, vr))
    }
  }
}

cat("\nDone. Outputs in:", out_root, "\n")
print(sessionInfo())  # package + GO annotation versions for the methods

# ---- Free memory ----
rm(list = ls())
gc(full = TRUE, reset = TRUE)
