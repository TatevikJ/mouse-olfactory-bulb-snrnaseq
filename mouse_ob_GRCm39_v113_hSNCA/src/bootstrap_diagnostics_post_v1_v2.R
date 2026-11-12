# Diagnostic: nFeature_RNA and nCount_RNA vs supertype_bootstrapping_probability,
# one PDF set per subclass in the TO_KEEP sheet of each pipeline's region-list
# xlsx. Run for both qc_post_v1 and qc_post_v2; intent is to see whether a
# supertype_bootstrapping_probability cutoff cleanly separates the sporadic
# ambiguous nuclei from the on-target ones, across all kept subclasses.

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---- Config ----
pipelines <- list(
  v1 = list(
    rds        = "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v1.rds",
    xlsx       = "meta/qc_post_v1_region_list.xlsx",
    plots_dir  = "outputs/qc_post_v1/plots",
    tsv_dir    = "outputs/qc_post_v1/tsv"
  ),
  v2 = list(
    rds        = "outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_v2.rds",
    xlsx       = "meta/qc_post_v2_region_list.xlsx",
    plots_dir  = "outputs/qc_post_v2/plots",
    tsv_dir    = "outputs/qc_post_v2/tsv"
  )
)

# ---- Plot theme ----
project_base_theme <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.title = element_blank(),
      legend.text = element_text(size = 10)
    )
}

bootstrap_thresh <- c(0.5, 0.7)

# Canonical condition palette (matches qc_pre.R)
condition_colors <- c(WT = "#59b2aa", ASO = "#d43f2e")

qc_metrics <- list(
  nFeature_RNA = list(col = "nFeature_RNA", label = "nFeature_RNA", short = "nFeature"),
  nCount_RNA   = list(col = "nCount_RNA",   label = "nCount_RNA",   short = "nCount")
)

# "041 OB-in Frmd7 Gaba" -> "subclass041" ; fallback for names without a
# numeric prefix: alphanumeric-only sanitized form.
subclass_slug <- function(name) {
  m <- regmatches(name, regexpr("^\\d+", name))
  if (length(m) == 0) {
    sprintf("subclass_%s", gsub("[^A-Za-z0-9]+", "_", name))
  } else {
    sprintf("subclass%s", m[1])
  }
}

# Pearson + Spearman, formatted as a single annotation string
fmt_corr <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) return(sprintf("n = %d (no corr)", length(x)))
  pe <- suppressWarnings(cor(x, y, method = "pearson"))
  sp <- suppressWarnings(cor(x, y, method = "spearman"))
  sprintf("n = %s\nPearson r = %.3f\nSpearman rho = %.3f",
          format(length(x), big.mark = ","), pe, sp)
}

# Per-supertype correlation table for the TSV side-output
per_supertype_corr <- function(df, x_col, y_col, group_col = "supertype_name") {
  median_y_col <- paste0("median_", y_col)
  df %>%
    dplyr::filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]])) %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      n_nuclei  = dplyr::n(),
      pearson_r = ifelse(n_nuclei >= 3 & sd(.data[[x_col]]) > 0 & sd(.data[[y_col]]) > 0,
                         suppressWarnings(cor(.data[[x_col]], .data[[y_col]], method = "pearson")),
                         NA_real_),
      spearman_rho = ifelse(n_nuclei >= 3 & sd(.data[[x_col]]) > 0 & sd(.data[[y_col]]) > 0,
                            suppressWarnings(cor(.data[[x_col]], .data[[y_col]], method = "spearman")),
                            NA_real_),
      !!median_y_col   := median(.data[[y_col]]),
      median_bootstrap  = median(.data[[x_col]]),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_nuclei))
}

# ---- Scatter recipe ----
make_scatter <- function(df, x_col, x_label, y_col, y_label, pipeline_tag,
                         subclass_label) {
  ann <- fmt_corr(df[[x_col]], df[[y_col]])
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = supertype_name)) +
    geom_point(size = 0.6, alpha = 0.35) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                color = "black", linewidth = 0.6, linetype = "solid") +
    geom_vline(xintercept = bootstrap_thresh, linetype = "dashed",
               color = "black", alpha = 0.6, linewidth = 0.4) +
    annotate("text", x = bootstrap_thresh, y = Inf, label = bootstrap_thresh,
             hjust = -0.15, vjust = 1.4, size = 3, color = "black") +
    annotate("text", x = -Inf, y = Inf, label = ann,
             hjust = -0.05, vjust = 1.15, size = 3.6, color = "black",
             lineheight = 0.95) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.02))) +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(
      x = x_label,
      y = y_label,
      title = sprintf("Subclass %s - %s", subclass_label, pipeline_tag),
      subtitle = "Restricted to TO_KEEP supertypes"
    ) +
    project_base_theme() +
    theme(
      legend.position = "right",
      plot.subtitle = element_text(size = 11, color = "gray30")
    )
}

# ---- Faceted scatter: one panel per supertype, per-facet LM + corr ----
# Styling mirrors qc_pre.R's doublet_score_histogram_per_sample: shared
# threshold data frame for vline+label, gray92 strip background, ncol=4 facets.
make_facet_scatter <- function(df, x_col, x_label, y_col, y_label, pipeline_tag,
                               subclass_label) {
  ann_df <- df %>%
    dplyr::filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]])) %>%
    dplyr::group_by(supertype_name) %>%
    dplyr::summarise(
      n_nuclei = dplyr::n(),
      pe = ifelse(n_nuclei >= 3 & sd(.data[[x_col]]) > 0 & sd(.data[[y_col]]) > 0,
                  suppressWarnings(cor(.data[[x_col]], .data[[y_col]], method = "pearson")),
                  NA_real_),
      sp = ifelse(n_nuclei >= 3 & sd(.data[[x_col]]) > 0 & sd(.data[[y_col]]) > 0,
                  suppressWarnings(cor(.data[[x_col]], .data[[y_col]], method = "spearman")),
                  NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = ifelse(is.na(pe),
                     sprintf("n = %s (no corr)", format(n_nuclei, big.mark = ",")),
                     sprintf("n = %s\nr = %.2f, rho = %.2f",
                             format(n_nuclei, big.mark = ","), pe, sp))
    )

  thresh_df <- data.frame(
    xint  = bootstrap_thresh,
    label = format(bootstrap_thresh, nsmall = 1)
  )

  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = condition)) +
    geom_point(size = 0.5, alpha = 0.4) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                color = "black", linewidth = 0.5, linetype = "solid") +
    geom_vline(data = thresh_df,
               aes(xintercept = xint),
               linetype = "dashed", color = "black",
               alpha = 0.8, linewidth = 0.5,
               inherit.aes = FALSE) +
    geom_text(data = thresh_df,
              aes(x = xint, y = Inf, label = label),
              hjust = -0.15, vjust = 1.4, size = 3, color = "black",
              inherit.aes = FALSE) +
    geom_text(data = ann_df,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.15, size = 3, color = "black",
              lineheight = 0.95, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.02))) +
    scale_color_manual(values = condition_colors) +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    facet_wrap(~ supertype_name, ncol = 4) +
    labs(
      x = x_label,
      y = y_label,
      title = sprintf("Subclass %s - %s, per supertype",
                      subclass_label, pipeline_tag),
      subtitle = "Restricted to TO_KEEP supertypes"
    ) +
    project_base_theme() +
    theme(
      legend.position = "bottom",
      plot.subtitle = element_text(size = 11, color = "gray30"),
      strip.background = element_rect(fill = "gray92", color = NA),
      strip.text = element_text(size = 11),
      panel.spacing.x = unit(1.5, "lines"),
      panel.spacing.y = unit(1, "lines")
    )
}

run_pipeline <- function(tag, cfg) {
  cat("\n========================================\n")
  cat("Pipeline:", tag, "\n")
  cat("========================================\n")
  stopifnot(file.exists(cfg$rds), file.exists(cfg$xlsx))
  dir.create(cfg$plots_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$tsv_dir,   recursive = TRUE, showWarnings = FALSE)

  cat("Reading TO_KEEP supertypes:", cfg$xlsx, "\n")
  keep_df <- suppressMessages(readxl::read_excel(cfg$xlsx, sheet = "TO_KEEP"))
  keep_supertypes <- unique(keep_df$supertype_name[!is.na(keep_df$supertype_name)])
  keep_subclasses <- keep_df %>%
    dplyr::filter(!is.na(subclass_name), !is.na(supertype_name)) %>%
    dplyr::distinct(subclass_name) %>%
    dplyr::pull(subclass_name) %>%
    sort()
  cat("  TO_KEEP supertypes:", length(keep_supertypes), "\n")
  cat("  TO_KEEP subclasses:", length(keep_subclasses), "\n")

  cat("Loading Seurat rds:", cfg$rds, "\n")
  obj <- readRDS(cfg$rds)
  md <- obj@meta.data
  rm(obj); invisible(gc())

  required_cols <- c("supertype_name", "subclass_name", "class_name",
                     "condition",
                     "nFeature_RNA", "nCount_RNA",
                     "supertype_bootstrapping_probability",
                     "subclass_bootstrapping_probability",
                     "class_bootstrapping_probability")
  missing_cols <- setdiff(required_cols, colnames(md))
  if (length(missing_cols)) stop("Missing metadata columns: ",
                                  paste(missing_cols, collapse = ", "))

  df_keep_all <- md %>%
    dplyr::filter(supertype_name %in% keep_supertypes)

  cat("Nuclei in object:               ", format(nrow(md), big.mark = ","), "\n")
  cat("Nuclei in TO_KEEP supertypes:   ", format(nrow(df_keep_all), big.mark = ","), "\n")

  if (nrow(df_keep_all) == 0) {
    cat("No nuclei in TO_KEEP supertypes; skipping pipeline\n")
    return(invisible(NULL))
  }

  # Subclass-scoped plots go into a dedicated subdir to keep plots/ uncluttered.
  subclass_plots_dir <- file.path(cfg$plots_dir, "subclass_bootstrap_diagnostics")
  dir.create(subclass_plots_dir, recursive = TRUE, showWarnings = FALSE)

  # Per-metric corr rows accumulate across subclasses, then write one TSV each.
  corr_rows <- setNames(vector("list", length(qc_metrics)), names(qc_metrics))

  for (subclass in keep_subclasses) {
    slug <- subclass_slug(subclass)
    df_sub <- df_keep_all %>%
      dplyr::filter(subclass_name == subclass)

    cat(sprintf("\n-- Subclass '%s' (%s): %s nuclei\n",
                subclass, slug, format(nrow(df_sub), big.mark = ",")))

    if (nrow(df_sub) == 0) {
      cat("   (empty; skipping)\n")
      next
    }
    df_sub$supertype_name <- factor(df_sub$supertype_name,
                                    levels = sort(unique(df_sub$supertype_name)))

    for (m_name in names(qc_metrics)) {
      m <- qc_metrics[[m_name]]
      p_supertype <- make_scatter(df_sub, "supertype_bootstrapping_probability",
                                  "Supertype bootstrapping probability",
                                  m$col, m$label, tag, subclass)
      out_pdf_supertype <- file.path(subclass_plots_dir,
        sprintf("%s_%s_vs_supertype_bootstrap_scatter.pdf", slug, m$short))
      ggsave(out_pdf_supertype, plot = p_supertype, width = 9, height = 6)
      cat("   wrote:", out_pdf_supertype, "\n")

      p_super_facet <- make_facet_scatter(df_sub, "supertype_bootstrapping_probability",
                                          "Supertype bootstrapping probability",
                                          m$col, m$label, tag, subclass)
      out_pdf_super_facet <- file.path(subclass_plots_dir,
        sprintf("%s_%s_vs_supertype_bootstrap_scatter_per_supertype.pdf",
                slug, m$short))
      ggsave(out_pdf_super_facet, plot = p_super_facet, width = 10, height = 7)
      cat("   wrote:", out_pdf_super_facet, "\n")

      corr_sub <- per_supertype_corr(df_sub, "supertype_bootstrapping_probability", m$col)
      corr_sub$subclass_name <- subclass
      corr_rows[[m_name]] <- dplyr::bind_rows(corr_rows[[m_name]], corr_sub)
    }
  }

  # One consolidated per-metric TSV with a subclass_name column.
  for (m_name in names(qc_metrics)) {
    m <- qc_metrics[[m_name]]
    corr_out <- corr_rows[[m_name]]
    if (is.null(corr_out) || nrow(corr_out) == 0) next
    corr_out <- corr_out %>% dplyr::relocate(subclass_name)
    out_tsv <- file.path(cfg$tsv_dir,
      sprintf("%s_vs_supertype_bootstrap_corr_per_supertype.tsv", m$short))
    write.table(corr_out, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\nwrote:", out_tsv, "\n")
  }

  rm(md, df_keep_all); invisible(gc())
}

# ---- Run both pipelines ----
for (tag in names(pipelines)) {
  run_pipeline(tag, pipelines[[tag]])
}

cat("\nDone.\n")
