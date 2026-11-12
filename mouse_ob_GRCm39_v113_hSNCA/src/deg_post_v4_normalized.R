# Differential expression analysis (ASO vs WT) on the v4 clustered protein-
# coding RDS. Runs three tracks, each toggled by a run_* switch below:
#
#   1. Pseudobulk DEG via scuttle::aggregateAcrossCells + scran::pseudoBulkDGE
#      (edgeR QL backend, filterByExpr gene filter). One pseudobulk sample per
#      (cell-type level, sample); pseudobulk samples with < 20 nuclei are
#      dropped. Design: ~ condition, coef "conditionASO" (WT reference).
#      Primary DEG result (mice are the replicates).
#      DEG: FDR < 0.05 & |logFC| > 0.5.
#
#   2. Single-nucleus Seurat::FindMarkers (Wilcoxon) within each cell type,
#      ASO vs WT, full gene table (logfc.threshold = 0, min.pct = 0.01).
#      Descriptive only: nuclei are treated as independent, so p-values are
#      not valid for genotype inference and no DEG calls are made.
#
#   3. Single-nucleus MAST hurdle model, Allen Institute mouse aging style
#      (github.com/AllenInstitute/mouse_aging_scRNAseq, 09a + functions/
#      degenes_age_mast_simplified_qcscore.R): downsample to <= 2000 nuclei
#      per condition, genes detected in > 10% of nuclei,
#      zlm(~ condition + log.gene.counts.0) on log2 data (Seurat data / log(2)),
#      LRT on conditionASO, Bonferroni. Nucleus-level screen, not a
#      mouse-level test. DEG: Bonferroni < 0.01 & |coef| > 0.5.
#
# All tracks run at class_name, subclass_name, supertype_name. Inputs are
# restricted to protein-coding genes (filtered upstream by
# normalize_post_v4.R, which also retains hSNCA-transgene).

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(SingleCellExperiment)
  library(scuttle)
  library(scran)
  library(edgeR)
  library(MAST)
  library(ggplot2)
  library(EnhancedVolcano)
  library(tibble)
})

# ---- Tracks to run ----
run_pseudobulk  <- TRUE
run_findmarkers <- TRUE
run_mast        <- TRUE

# ---- Inputs / outputs ----
in_rds   <- "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v4_clustered_protein_coding.rds"
out_root <- "outputs/deg_post_v4_normalized"
# Each track writes to its own subtree with rds/, tsv/, plots/.
track_dirs <- function(track) {
  d <- list(rds   = file.path(out_root, track, "rds"),
            tsv   = file.path(out_root, track, "tsv"),
            plots = file.path(out_root, track, "plots"))
  for (x in d) dir.create(x, recursive = TRUE, showWarnings = FALSE)
  d
}
pb_dirs   <- track_dirs("pseudobulk")
fm_dirs   <- track_dirs("singlecell_findmarkers")
mast_dirs <- track_dirs("singlecell_mast")

# ---- DEG params ----
levels_to_run           <- c("class_name", "subclass_name", "supertype_name")
min_cells_per_pb_sample <- 20    # drop (celltype, sample) pseudobulks below this
pb_fdr_cutoff           <- 0.05  # pseudobulk DEG: FDR < 0.05
fc_cutoff               <- 0.5   # |log2 FC| (pseudobulk logFC, MAST coef); ~1.41-fold
min_cells_per_condition <- 50    # single-nucleus tracks: skip cell type below this
fm_min_pct              <- 0.01  # FindMarkers: Seurat v5 default
mast_max_cells          <- 2000  # MAST: downsample cap per condition (Allen)
mast_freq_expressed     <- 0.1   # MAST: gene detected in > 10% of nuclei (Allen)
mast_padj_cutoff        <- 0.01  # MAST DEG: Bonferroni < 0.01 (Allen)
mast_seed               <- 2010
mast_cores              <- 8     # match sbatch --cpus-per-task
options(mc.cores = mast_cores)

# ---- Load v4 protein-coding RDS ----
cat("Loading:", in_rds, "\n")
dataSeurat <- readRDS(in_rds)
dataSeurat <- JoinLayers(dataSeurat)  # no-op if already joined; $counts/$data need one layer
if (!"condition" %in% colnames(dataSeurat@meta.data)) {
  dataSeurat$condition <- ifelse(grepl("^ASO", dataSeurat$sample), "ASO", "WT")
}
dataSeurat$condition <- factor(dataSeurat$condition, levels = c("WT", "ASO"))
cat("Nuclei:", ncol(dataSeurat), " Genes:", nrow(dataSeurat), "\n")

sanitize <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)

# Volcano shared by all tracks; cutoff lines match each track's DEG definition.
make_volcano <- function(res_df, cell_type, x, y, p_cutoff, title_suffix, out_dir) {
  res_df <- as.data.frame(res_df)
  res_df <- res_df[!is.na(res_df[[x]]) & !is.na(res_df[[y]]), ]
  if (nrow(res_df) < 2) return(invisible(NULL))
  p <- EnhancedVolcano(
    res_df,
    lab             = res_df$gene,
    x               = x,
    y               = y,
    title           = sprintf("%s (ASO vs WT, %s)", cell_type, title_suffix),
    pCutoff         = p_cutoff,
    FCcutoff        = fc_cutoff,
    pointSize       = 2.0,
    labSize         = 3.5,
    drawConnectors  = TRUE,
    widthConnectors = 0.5,
    max.overlaps    = 40,
    colAlpha        = 0.7,
    legendPosition  = "right"
  )
  out_pdf <- file.path(out_dir, sprintf("%s__%s.pdf", sanitize(cell_type), sanitize(title_suffix)))
  ggsave(out_pdf, plot = p, width = 10, height = 8)
  cat("  wrote:", out_pdf, "\n")
}

# Save per-level RDS (named list of per-cell-type tables), combined TSV, volcanoes.
save_track <- function(res_list, level_col, dirs, prefix, x, y, p_cutoff, title_suffix) {
  if (is.null(res_list) || length(res_list) == 0) return(invisible(NULL))
  out <- file.path(dirs$rds, sprintf("de_results_%s_%s.rds", prefix, level_col))
  saveRDS(res_list, out); cat("wrote:", out, "\n")

  all_df <- dplyr::bind_rows(lapply(names(res_list), function(ct) {
    df <- as.data.frame(res_list[[ct]])
    if (nrow(df) == 0) return(NULL)
    df$cell_type <- ct
    df
  }))
  out_tsv <- file.path(dirs$tsv, sprintf("deg_%s_%s.tsv", prefix, level_col))
  write.table(all_df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("wrote:", out_tsv, "\n")

  plot_dir <- file.path(dirs$plots, level_col)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  for (ct in names(res_list)) {
    make_volcano(res_list[[ct]], ct, x, y, p_cutoff, title_suffix, plot_dir)
  }
}

# ============================================================
# 1. Pseudobulk DEG
# ============================================================
run_pseudobulk_deg <- function(sce, level_col) {
  cat("\n--- Pseudobulk DEG by", level_col, "---\n")
  pb <- scuttle::aggregateAcrossCells(
    sce,
    ids = colData(sce)[, c(level_col, "sample")]
  )
  cat("  initial pseudobulk: ", ncol(pb), " (celltype, sample) groups\n")
  pb <- pb[, pb$ncells >= min_cells_per_pb_sample]
  cat("  after ncells >= ", min_cells_per_pb_sample, " filter: ",
      ncol(pb), " groups\n", sep = "")
  if (ncol(pb) == 0) {
    cat("  No groups left; skipping.\n")
    return(NULL)
  }
  cat("  condition split:\n")
  print(table(pb$condition, pb[[level_col]]))

  de <- scran::pseudoBulkDGE(
    pb,
    label  = pb[[level_col]],
    design = ~ condition,
    coef   = "conditionASO"
  )
  cat("  cell types analyzed: ", length(de), "\n", sep = "")
  failed <- metadata(de)$failed
  if (length(failed) > 0) {
    cat("  failed for: ", paste(failed, collapse = ", "), "\n", sep = "")
  }

  # Plain data.frames with gene + significant columns (NA rows = filterByExpr).
  res <- lapply(as.list(de), function(df) {
    df <- as.data.frame(df)
    df$gene <- rownames(df)
    df$significant <- !is.na(df$FDR) & df$FDR < pb_fdr_cutoff & abs(df$logFC) > fc_cutoff
    df
  })
  for (ct in names(res)) {
    cat(sprintf("  %s: %d DEGs\n", ct, sum(res[[ct]]$significant)))
  }
  res
}

if (run_pseudobulk) {
  cat("\n[1/3] Pseudobulk...\n")
  metadata <- dataSeurat@meta.data
  for (lvl in levels_to_run) metadata[[lvl]] <- factor(metadata[[lvl]])
  metadata$sample <- factor(metadata$sample)
  sce <- SingleCellExperiment(
    assays  = list(counts = dataSeurat@assays$RNA$counts),
    colData = metadata
  )
  cat("  SCE dims:", nrow(sce), "genes x", ncol(sce), "nuclei\n")

  for (lvl in levels_to_run) {
    res <- run_pseudobulk_deg(sce, lvl)
    save_track(res, lvl, pb_dirs, "pseudobulk", "logFC", "FDR", pb_fdr_cutoff, "pseudobulk")
  }
  rm(sce, metadata); gc()
}

# ============================================================
# 2. Single-nucleus FindMarkers (Wilcoxon), full table, descriptive
# ============================================================
# Seurat uses presto for a fast Wilcoxon when installed; without it this is slow.
run_findmarkers_deg <- function(seu, level_col) {
  cat("\n--- FindMarkers by", level_col, "---\n")
  out <- list()
  cell_types <- sort(unique(as.character(seu@meta.data[[level_col]])))
  for (ct in cell_types) {
    keep <- which(seu@meta.data[[level_col]] == ct)
    sub  <- seu[, keep]
    tab  <- table(sub$condition)
    if (any(tab < min_cells_per_condition)) {
      cat("  Skipping '", ct, "' (nuclei per condition: ",
          paste(names(tab), tab, sep = "=", collapse = ", "), ")\n", sep = "")
      next
    }
    res <- tryCatch({
      Seurat::FindMarkers(
        sub,
        ident.1         = "ASO",
        ident.2         = "WT",
        group.by        = "condition",
        logfc.threshold = 0,
        min.pct         = fm_min_pct,
        verbose         = FALSE
      )
    }, error = function(e) {
      cat("  FindMarkers failed for '", ct, "': ", conditionMessage(e), "\n", sep = "")
      NULL
    })
    if (is.null(res)) next
    res$gene <- rownames(res)
    out[[ct]] <- res
    cat(sprintf("  %s: %d genes (nuclei ASO=%d WT=%d)\n",
                ct, nrow(res), tab[["ASO"]], tab[["WT"]]))
  }
  out
}

if (run_findmarkers) {
  cat("\n[2/3] FindMarkers (Wilcoxon)...\n")
  cat("  presto installed:", requireNamespace("presto", quietly = TRUE), "\n")
  for (lvl in levels_to_run) {
    res <- run_findmarkers_deg(dataSeurat, lvl)
    save_track(res, lvl, fm_dirs, "findmarkers", "avg_log2FC", "p_val_adj", 0.05, "findmarkers")
  }
}

# ============================================================
# 3. Single-nucleus MAST (Allen style)
# ============================================================
run_mast_deg <- function(seu, level_col) {
  cat("\n--- MAST by", level_col, "---\n")
  out <- list()
  cell_types <- sort(unique(as.character(seu@meta.data[[level_col]])))
  for (ct in cell_types) {
    keep <- which(seu@meta.data[[level_col]] == ct)
    tab  <- table(seu$condition[keep])
    if (any(tab < min_cells_per_condition)) {
      cat("  Skipping '", ct, "' (nuclei per condition: ",
          paste(names(tab), tab, sep = "=", collapse = ", "), ")\n", sep = "")
      next
    }

    # Downsample each condition to mast_max_cells; seed per cell type so the
    # draw does not depend on loop order.
    set.seed(mast_seed)
    cells_ct <- colnames(seu)[keep]
    cond_ct  <- seu$condition[keep]
    sampled  <- unlist(lapply(c("WT", "ASO"), function(g) {
      x <- cells_ct[cond_ct == g]
      sample(x, min(length(x), mast_max_cells))
    }))
    sub <- seu[, sampled]

    # Seurat data is ln(1 + CP10k); divide by log(2) so coef is in log2 units.
    expr <- as.matrix(sub[["RNA"]]$data) / log(2)
    cdat <- data.frame(
      wellKey   = colnames(sub),
      condition = factor(sub$condition, levels = c("WT", "ASO")),
      log.gene.counts.0 = as.numeric(scale(log(sub$nFeature_RNA))),
      row.names = colnames(sub)
    )
    fdat <- data.frame(primerid = rownames(expr), row.names = rownames(expr))
    sca  <- MAST::FromMatrix(exprsArray = expr, cData = cdat, fData = fdat)
    sca  <- sca[MAST::freq(sca) > mast_freq_expressed, ]
    rm(expr); gc()

    res <- tryCatch({
      zfit <- MAST::zlm(~ condition + log.gene.counts.0, sca)
      zsum <- summary(zfit, doLRT = "conditionASO")$datatable
      zsum <- as.data.frame(zsum)
      p_part  <- zsum[zsum$contrast == "conditionASO" & zsum$component == "H",
                      c("primerid", "Pr(>Chisq)")]
      fc_part <- zsum[zsum$contrast == "conditionASO" & zsum$component == "logFC",
                      c("primerid", "coef", "ci.lo", "ci.hi")]
      dplyr::inner_join(p_part, fc_part, by = "primerid") %>%
        dplyr::rename(gene = primerid, pvalue = `Pr(>Chisq)`) %>%
        dplyr::mutate(
          padjust     = p.adjust(pvalue, method = "bonferroni"),
          significant = !is.na(padjust) & padjust < mast_padj_cutoff & abs(coef) > fc_cutoff,
          n_wt        = sum(sub$condition == "WT"),
          n_aso       = sum(sub$condition == "ASO")
        ) %>%
        dplyr::arrange(padjust, dplyr::desc(abs(coef)))
    }, error = function(e) {
      cat("  MAST failed for '", ct, "': ", conditionMessage(e), "\n", sep = "")
      NULL
    })
    if (is.null(res)) next
    out[[ct]] <- res
    cat(sprintf("  %s: %d genes tested, %d DEGs (nuclei ASO=%d WT=%d)\n",
                ct, nrow(res), sum(res$significant), res$n_aso[1], res$n_wt[1]))
    rm(sub, sca); gc()
  }
  out
}

if (run_mast) {
  cat("\n[3/3] MAST...\n")
  for (lvl in levels_to_run) {
    res <- run_mast_deg(dataSeurat, lvl)
    save_track(res, lvl, mast_dirs, "mast", "coef", "padjust", mast_padj_cutoff, "mast")
  }
}

cat("\nDone. Outputs in:", out_root, "\n")

# ---- Free memory ----
rm(list = ls())
gc(full = TRUE, reset = TRUE)
