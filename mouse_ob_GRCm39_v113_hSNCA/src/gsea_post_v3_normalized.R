# Gene set enrichment analysis (GSEA, GO:BP) of the ASO vs WT DEG results
# from deg_post_v3_normalized.R. clusterProfiler::gseGO per cell type on the
# full ranked list of tested genes.
#
#   Tracks:   pseudobulk (edgeR logFC / PValue; primary) and singlecell_mast
#             (MAST hurdle coef / pvalue; nucleus-level, supplementary).
#   Rankings: logfc (effect size; primary) and
#             signed_logp (sign(effect) * -log10(p); Reimand 2019, sensitivity).
#   Variants: all_genes, and no_Pcdh11x (Pcdh11x removed before ranking).
#   Levels:   class_name, subclass_name, supertype_name.
#
# Genes not tested in the DEG step (NA rows, e.g. pseudobulk filterByExpr) are
# excluded. p = 0 is set to the smallest non-zero p before -log10. fgsea
# multilevel with eps = 0 (no floor on small p-values) and a fixed seed; GO
# terms of 10-500 genes; BH across terms within each cell type. All terms are
# saved (pvalueCutoff = 1); significance is p.adjust < 0.05. hSNCA-transgene
# has no GO annotation and drops out.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
})

# ---- Inputs / outputs ----
in_root  <- "outputs/deg_post_v3_normalized"
out_root <- "outputs/gsea_post_v3_normalized"
# Track -> DEG file prefix, effect-size column, raw p-value column.
tracks <- list(
  pseudobulk      = list(prefix = "pseudobulk", fc_col = "logFC", p_col = "PValue"),
  singlecell_mast = list(prefix = "mast",       fc_col = "coef",  p_col = "pvalue")
)
rankings <- c("logfc", "signed_logp")
# Variant -> genes removed before analysis.
variants <- list(all_genes = character(0), no_Pcdh11x = "Pcdh11x")
# Each (track, ranking, variant) writes to its own subtree with rds/, tsv/, plots/.
track_dirs <- function(track, ranking, variant) {
  d <- list(rds   = file.path(out_root, track, ranking, variant, "rds"),
            tsv   = file.path(out_root, track, ranking, variant, "tsv"),
            plots = file.path(out_root, track, ranking, variant, "plots"))
  for (x in d) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  d
}

# ---- GSEA params ----
levels_to_run    <- c("class_name", "subclass_name", "supertype_name")
min_gs_size      <- 10    # GO term size bounds (clusterProfiler defaults)
max_gs_size      <- 500
min_ranked_genes <- 10    # skip cell types with fewer tested genes
padj_cutoff      <- 0.05  # significant term: BH p.adjust < 0.05
n_dot_terms      <- 20    # dotplot: top terms per NES sign
n_running_terms  <- 5     # running-score plot: top terms
n_heatmap_terms  <- 5     # NES heatmap: top terms per cell type
gsea_seed        <- 2010
gsea_cores       <- 16    # match sbatch --cpus-per-task
BiocParallel::register(BiocParallel::MulticoreParam(gsea_cores))

go_symbols <- AnnotationDbi::keys(org.Mm.eg.db, keytype = "SYMBOL")
sanitize   <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)

# ============================================================
# GSEA per cell type
# ============================================================
# Tested genes, named by symbol, sorted decreasing by the ranking score.
make_ranked <- function(df, tr, ranking, drop_genes) {
  df <- df[!is.na(df[[tr$fc_col]]) & !is.na(df[[tr$p_col]]) & !df$gene %in% drop_genes, ]
  score <- df[[tr$fc_col]]
  if (ranking == "signed_logp") {
    p <- pmax(df[[tr$p_col]], min(df[[tr$p_col]][df[[tr$p_col]] > 0]))
    score <- sign(score) * -log10(p)
  }
  sort(setNames(score, df$gene), decreasing = TRUE)
}

run_gsea_level <- function(deg_list, tr, ranking, drop_genes) {
  out <- list()
  for (ct in names(deg_list)) {
    gl <- make_ranked(as.data.frame(deg_list[[ct]]), tr, ranking, drop_genes)
    if (length(gl) < min_ranked_genes) {
      cat("  Skipping '", ct, "' (", length(gl), " ranked genes)\n", sep = "")
      next
    }
    set.seed(gsea_seed)
    gse <- tryCatch({
      clusterProfiler::gseGO(
        geneList      = gl,
        ont           = "BP",
        OrgDb         = org.Mm.eg.db,
        keyType       = "SYMBOL",
        minGSSize     = min_gs_size,
        maxGSSize     = max_gs_size,
        eps           = 0,
        pvalueCutoff  = 1,
        pAdjustMethod = "BH",
        seed          = TRUE,
        verbose       = FALSE
      )
    }, error = function(e) {
      cat("  gseGO failed for '", ct, "': ", conditionMessage(e), "\n", sep = "")
      NULL
    })
    if (is.null(gse)) next
    gse@geneSets <- gse@geneSets[gse@result$ID]  # keep tested terms only (RDS size)
    out[[ct]] <- gse
    cat(sprintf("  %s: %d ranked genes (%d GO symbols), %d terms, %d significant\n",
                ct, length(gl), sum(names(gl) %in% go_symbols),
                nrow(gse@result), sum(gse@result$p.adjust < padj_cutoff)))
  }
  out
}

# ============================================================
# Plots
# ============================================================
# Dotplot of significant terms, split by NES sign.
plot_dot <- function(gse, ct, title_suffix, out_dir) {
  sig <- gse
  sig@result <- gse@result[gse@result$p.adjust < padj_cutoff, ]
  if (nrow(sig@result) == 0) return(invisible(NULL))
  p <- enrichplot::dotplot(sig, showCategory = n_dot_terms, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle(sprintf("%s (GO:BP GSEA, %s)", ct, title_suffix)) +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(out_dir, sprintf("%s__dotplot.pdf", sanitize(ct))), p, width = 12, height = 8)
}

# Running enrichment score of the top significant terms.
plot_running <- function(gse, ct, title_suffix, out_dir) {
  res <- gse@result[order(gse@result$p.adjust), ]
  ids <- head(res$ID[res$p.adjust < padj_cutoff], n_running_terms)
  if (length(ids) == 0) return(invisible(NULL))
  p <- enrichplot::gseaplot2(gse, geneSetID = ids,
                             title = sprintf("%s (top %d terms, %s)", ct, length(ids), title_suffix))
  pdf(file.path(out_dir, sprintf("%s__running_score.pdf", sanitize(ct))), width = 10, height = 7)
  print(p)
  dev.off()
}

# NES heatmap across cell types: top significant terms per cell type, * = significant.
plot_nes_heatmap <- function(tsv_df, level_col, title_suffix, out_pdf) {
  sig <- tsv_df[tsv_df$p.adjust < padj_cutoff, ]
  if (nrow(sig) == 0) return(invisible(NULL))
  terms <- sig %>%
    group_by(cell_type) %>%
    slice_min(p.adjust, n = n_heatmap_terms, with_ties = FALSE) %>%
    pull(Description) %>%
    unique()
  hm <- tsv_df %>%
    filter(Description %in% terms) %>%
    mutate(Description = reorder(Description, NES, mean),
           sig = p.adjust < padj_cutoff)
  p <- ggplot(hm, aes(cell_type, Description, fill = NES)) +
    geom_tile(color = "white") +
    geom_text(data = filter(hm, sig), label = "*", size = 3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
    labs(x = NULL, y = NULL, caption = sprintf("* p.adjust < %s", padj_cutoff),
         title = sprintf("GO:BP GSEA NES by %s (%s)", level_col, title_suffix)) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(),
          plot.title.position = "plot",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(out_pdf, p, limitsize = FALSE,
         width  = max(8, 4 + 0.25 * n_distinct(hm$cell_type)),
         height = max(6, 2 + 0.15 * length(terms)))
  cat("wrote:", out_pdf, "\n")
}

# Save per-level RDS (named list of gseaResult), combined TSV, plots.
save_level <- function(res_list, level_col, dirs, title_suffix) {
  if (length(res_list) == 0) return(invisible(NULL))
  out <- file.path(dirs$rds, sprintf("gsea_results_%s.rds", level_col))
  saveRDS(res_list, out); cat("wrote:", out, "\n")

  tsv_df <- dplyr::bind_rows(lapply(names(res_list), function(ct) {
    df <- res_list[[ct]]@result
    if (nrow(df) == 0) return(NULL)
    df$cell_type <- ct
    df
  }))
  out_tsv <- file.path(dirs$tsv, sprintf("gsea_%s.tsv", level_col))
  write.table(tsv_df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("wrote:", out_tsv, "\n")

  plot_dir <- file.path(dirs$plots, level_col)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  for (ct in names(res_list)) {
    plot_dot(res_list[[ct]], ct, title_suffix, plot_dir)
    plot_running(res_list[[ct]], ct, title_suffix, plot_dir)
  }
  plot_nes_heatmap(tsv_df, level_col, title_suffix, file.path(plot_dir, "nes_heatmap.pdf"))
}

# ============================================================
# Run all tracks x levels x rankings x variants
# ============================================================
for (tr in names(tracks)) {
  for (lvl in levels_to_run) {
    in_rds <- file.path(in_root, tr, "rds",
                        sprintf("de_results_%s_%s.rds", tracks[[tr]]$prefix, lvl))
    cat("\nLoading:", in_rds, "\n")
    deg_list <- readRDS(in_rds)
    for (rk in rankings) {
      for (vr in names(variants)) {
        cat(sprintf("\n--- GSEA: %s / %s / %s / %s ---\n", tr, rk, vr, lvl))
        res <- run_gsea_level(deg_list, tracks[[tr]], rk, variants[[vr]])
        save_level(res, lvl, track_dirs(tr, rk, vr), sprintf("%s, %s, %s", tr, rk, vr))
      }
    }
  }
}

cat("\nDone. Outputs in:", out_root, "\n")
print(sessionInfo())  # package + GO annotation versions for the methods

# ---- Free memory ----
rm(list = ls())
gc(full = TRUE, reset = TRUE)
