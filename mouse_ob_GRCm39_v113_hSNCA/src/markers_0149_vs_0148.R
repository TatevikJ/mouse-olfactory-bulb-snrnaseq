# FindMarkers between supertypes 0149 OB Trdn Gaba_2 and 0148 OB Trdn Gaba_1
# (both subclass 040 OB Trdn Gaba) on the v3 and v4 clustered objects,
# to find a staining marker separating them.
# Positive avg_log2FC = higher in 0149, negative = higher in 0148.
# Same test and thresholds as markers_levels_post_v{3,4}.R.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

st_1 <- "0149 OB Trdn Gaba_2"
st_2 <- "0148 OB Trdn Gaba_1"

# TRUE = protein-coding object (markers_0149_vs_0148_protein_coding.csv);
# FALSE = full object incl. lncRNA/Gm genes (markers_0149_vs_0148.csv)
protein_coding <- FALSE
suffix <- if (protein_coding) "_protein_coding" else ""

for (v in c("v3", "v4")) {
  in_rds  <- sprintf("outputs/rds/aso_wt_data_annotated_normalized_qc_filtered_%s_clustered%s.rds", v, suffix)
  out_csv <- sprintf("outputs/qc_post_%s_normalized/tsv/markers_0149_vs_0148%s.csv", v, suffix)

  cat("Loading:", in_rds, "\n")
  dataSeurat <- readRDS(in_rds)
  dataSeurat <- subset(dataSeurat, subset = supertype_name %in% c(st_1, st_2))
  gc()
  Idents(dataSeurat) <- dataSeurat$supertype_name
  print(table(Idents(dataSeurat)))

  markers <- FindMarkers(dataSeurat,
                         ident.1         = st_1,
                         ident.2         = st_2,
                         only.pos        = FALSE,
                         min.pct         = 0.25,
                         logfc.threshold = 0.25,
                         verbose         = FALSE)
  markers <- markers %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::arrange(dplyr::desc(avg_log2FC))
  write.csv(markers, out_csv, row.names = FALSE)
  cat("  found ", nrow(markers), " markers; wrote ", out_csv, "\n", sep = "")

  rm(dataSeurat, markers)
  gc()
}

cat("\nDone.\n")
