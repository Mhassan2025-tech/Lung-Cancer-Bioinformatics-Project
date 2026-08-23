# ============================================================
# Lung Cancer Bioinformatics Project
# Step 6: Gene Annotation + DEG Heatmap
# Dataset: GSE10072
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(tidyverse)
library(AnnotationDbi)
library(hgu133a.db)
library(pheatmap)

# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------

expr <- readRDS(
  "data/raw/GSE10072_expression.rds"
)

group <- readRDS(
  "data/processed/sample_groups.rds"
)

deg <- read.csv(
  "results/tables/all_differential_expression_results.csv",
  stringsAsFactors = FALSE
)

sample_info <- readRDS(
  "data/processed/sample_information.rds"
)

# ------------------------------------------------------------
# 3. Verify data
# ------------------------------------------------------------

cat("============================================\n")
cat("STEP 6: GENE ANNOTATION + HEATMAP\n")
cat("============================================\n")

cat("\nExpression dimensions:\n")
print(dim(expr))

cat("\nDEG dimensions:\n")
print(dim(deg))

# ------------------------------------------------------------
# 4. Probe ID → Gene Symbol
# ------------------------------------------------------------

cat("\nAnnotating probes...\n")

gene_symbols <- mapIds(
  hgu133a.db,
  keys = deg$Probe_ID,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

deg$Gene_Symbol <- unname(gene_symbols)

# ------------------------------------------------------------
# 5. Remove unannotated probes
# ------------------------------------------------------------

deg_annotated <- deg %>%
  
  filter(
    !is.na(Gene_Symbol),
    Gene_Symbol != ""
  )

cat(
  "\nAnnotated probes:",
  nrow(deg_annotated),
  "\n"
)

# ------------------------------------------------------------
# 6. Select significant DEGs
# ------------------------------------------------------------

significant_deg <- deg_annotated %>%
  
  filter(
    adj.P.Val < 0.05,
    abs(logFC) >= 1
  ) %>%
  
  arrange(
    adj.P.Val
  )

cat(
  "Significant annotated DEGs:",
  nrow(significant_deg),
  "\n"
)

if (nrow(significant_deg) == 0) {
  
  stop(
    "ERROR: No significant annotated DEGs found."
  )
  
}

# ------------------------------------------------------------
# 7. Remove duplicate gene symbols
# ------------------------------------------------------------

significant_deg_unique <- significant_deg %>%
  
  group_by(
    Gene_Symbol
  ) %>%
  
  slice_min(
    order_by = adj.P.Val,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup()

cat(
  "Unique significant genes:",
  nrow(significant_deg_unique),
  "\n"
)

# ------------------------------------------------------------
# 8. Select top 30 genes
# ------------------------------------------------------------

top_n <- min(
  30,
  nrow(significant_deg_unique)
)

top_genes <- significant_deg_unique %>%
  
  arrange(
    adj.P.Val
  ) %>%
  
  slice_head(
    n = top_n
  )

cat(
  "Genes selected for heatmap:",
  top_n,
  "\n"
)

# ------------------------------------------------------------
# 9. Save annotated results
# ------------------------------------------------------------

write.csv(
  deg_annotated,
  "results/tables/annotated_differential_expression_results.csv",
  row.names = FALSE
)

write.csv(
  top_genes,
  "results/tables/top_30_annotated_DEGs.csv",
  row.names = FALSE
)

# ============================================================
# HEATMAP MATRIX
# ============================================================

# ------------------------------------------------------------
# 10. Extract probes
# ------------------------------------------------------------

top_probe_ids <- top_genes$Probe_ID

# Make sure probes exist

top_probe_ids <- top_probe_ids[
  top_probe_ids %in% rownames(expr)
]

if (length(top_probe_ids) < 2) {
  
  stop(
    "ERROR: Not enough valid probes available for heatmap."
  )
  
}

# ------------------------------------------------------------
# 11. Extract expression matrix
# ------------------------------------------------------------

heatmap_expr <- expr[
  top_probe_ids,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 12. Match gene symbols to probes
# ------------------------------------------------------------

probe_to_gene <- top_genes$Gene_Symbol[
  match(
    top_probe_ids,
    top_genes$Probe_ID
  )
]

# Remove missing gene symbols

valid <- !is.na(probe_to_gene) &
  probe_to_gene != ""

heatmap_expr <- heatmap_expr[
  valid,
  ,
  drop = FALSE
]

probe_to_gene <- probe_to_gene[
  valid
]

# ------------------------------------------------------------
# 13. Make row names unique
# ------------------------------------------------------------

rownames(heatmap_expr) <- make.unique(
  probe_to_gene
)

# ------------------------------------------------------------
# 14. Remove rows with NA / Inf
# ------------------------------------------------------------

valid_rows <- apply(
  heatmap_expr,
  1,
  function(x) {
    all(
      is.finite(x)
    )
  }
)

heatmap_expr <- heatmap_expr[
  valid_rows,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 15. Remove zero-variance genes
# ------------------------------------------------------------

gene_variance <- apply(
  heatmap_expr,
  1,
  var
)

heatmap_expr <- heatmap_expr[
  gene_variance > 0,
  ,
  drop = FALSE
]

cat(
  "\nFinal genes for heatmap:",
  nrow(heatmap_expr),
  "\n"
)

if (nrow(heatmap_expr) < 2) {
  
  stop(
    "ERROR: Fewer than 2 usable genes remain for heatmap."
  )
  
}

# ============================================================
# Z-SCORE
# ============================================================

cat("\nCalculating gene-wise Z-scores...\n")

heatmap_scaled <- t(
  scale(
    t(heatmap_expr)
  )
)

# ------------------------------------------------------------
# Check scaled matrix
# ------------------------------------------------------------

if (
  any(
    !is.finite(
      heatmap_scaled
    )
  )
) {
  
  stop(
    "ERROR: NA/Inf values detected after scaling."
  )
  
}

# ============================================================
# SAMPLE ANNOTATION
# ============================================================

annotation_col <- data.frame(
  Group = group
)

rownames(annotation_col) <- colnames(
  expr
)

# Make sure annotation order matches matrix

annotation_col <- annotation_col[
  colnames(heatmap_scaled),
  ,
  drop = FALSE
]

# ============================================================
# HEATMAP
# ============================================================

cat("\nCreating heatmap...\n")

heatmap_file <- paste0(
  "results/figures/08_top_DEG_heatmap.png"
)

png(
  filename = heatmap_file,
  width = 2800,
  height = 2200,
  res = 300
)

pheatmap(
  mat = heatmap_scaled,
  
  annotation_col = annotation_col,
  
  show_colnames = FALSE,
  
  show_rownames = TRUE,
  
  fontsize_row = 8,
  
  clustering_method = "complete",
  
  border_color = NA,
  
  main = paste0(
    "Top ",
    nrow(heatmap_scaled),
    " Differentially Expressed Genes"
  )
)

dev.off()

# ------------------------------------------------------------
# Verify file
# ------------------------------------------------------------

if (!file.exists(heatmap_file)) {
  
  stop(
    "ERROR: Heatmap file was not created."
  )
  
}

file_size <- file.info(
  heatmap_file
)$size

cat(
  "\nHeatmap file size:",
  file_size,
  "bytes\n"
)

if (file_size == 0) {
  
  stop(
    "ERROR: Heatmap file is 0 bytes."
  )
  
}

# ============================================================
# SAVE MATRIX
# ============================================================

write.csv(
  heatmap_scaled,
  "results/tables/top_DEG_heatmap_matrix.csv"
)

saveRDS(
  heatmap_scaled,
  "data/processed/top_DEG_heatmap_matrix.rds"
)

# ============================================================
# FINAL STATUS
# ============================================================

cat("\n============================================\n")
cat("STEP 6 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat(
  "\nAnnotated probes:",
  nrow(deg_annotated),
  "\n"
)

cat(
  "Significant annotated DEGs:",
  nrow(significant_deg),
  "\n"
)

cat(
  "Final heatmap genes:",
  nrow(heatmap_scaled),
  "\n"
)

cat(
  "Heatmap saved:",
  heatmap_file,
  "\n"
)

cat(
  "File size:",
  file_size,
  "bytes\n"
)

cat("\nNext step: Machine Learning classification.\n")# ============================================================
# Lung Cancer Bioinformatics Project
# Step 6: Gene Annotation + DEG Heatmap
# Dataset: GSE10072
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(tidyverse)
library(AnnotationDbi)
library(hgu133a.db)
library(pheatmap)

# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------

expr <- readRDS(
  "data/raw/GSE10072_expression.rds"
)

group <- readRDS(
  "data/processed/sample_groups.rds"
)

deg <- read.csv(
  "results/tables/all_differential_expression_results.csv",
  stringsAsFactors = FALSE
)

sample_info <- readRDS(
  "data/processed/sample_information.rds"
)

# ------------------------------------------------------------
# 3. Verify data
# ------------------------------------------------------------

cat("============================================\n")
cat("STEP 6: GENE ANNOTATION + HEATMAP\n")
cat("============================================\n")

cat("\nExpression dimensions:\n")
print(dim(expr))

cat("\nDEG dimensions:\n")
print(dim(deg))

# ------------------------------------------------------------
# 4. Probe ID → Gene Symbol
# ------------------------------------------------------------

cat("\nAnnotating probes...\n")

gene_symbols <- mapIds(
  hgu133a.db,
  keys = deg$Probe_ID,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

deg$Gene_Symbol <- unname(gene_symbols)

# ------------------------------------------------------------
# 5. Remove unannotated probes
# ------------------------------------------------------------

deg_annotated <- deg %>%
  
  filter(
    !is.na(Gene_Symbol),
    Gene_Symbol != ""
  )

cat(
  "\nAnnotated probes:",
  nrow(deg_annotated),
  "\n"
)

# ------------------------------------------------------------
# 6. Select significant DEGs
# ------------------------------------------------------------

significant_deg <- deg_annotated %>%
  
  filter(
    adj.P.Val < 0.05,
    abs(logFC) >= 1
  ) %>%
  
  arrange(
    adj.P.Val
  )

cat(
  "Significant annotated DEGs:",
  nrow(significant_deg),
  "\n"
)

if (nrow(significant_deg) == 0) {
  
  stop(
    "ERROR: No significant annotated DEGs found."
  )
  
}

# ------------------------------------------------------------
# 7. Remove duplicate gene symbols
# ------------------------------------------------------------

significant_deg_unique <- significant_deg %>%
  
  group_by(
    Gene_Symbol
  ) %>%
  
  slice_min(
    order_by = adj.P.Val,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup()

cat(
  "Unique significant genes:",
  nrow(significant_deg_unique),
  "\n"
)

# ------------------------------------------------------------
# 8. Select top 30 genes
# ------------------------------------------------------------

top_n <- min(
  30,
  nrow(significant_deg_unique)
)

top_genes <- significant_deg_unique %>%
  
  arrange(
    adj.P.Val
  ) %>%
  
  slice_head(
    n = top_n
  )

cat(
  "Genes selected for heatmap:",
  top_n,
  "\n"
)

# ------------------------------------------------------------
# 9. Save annotated results
# ------------------------------------------------------------

write.csv(
  deg_annotated,
  "results/tables/annotated_differential_expression_results.csv",
  row.names = FALSE
)

write.csv(
  top_genes,
  "results/tables/top_30_annotated_DEGs.csv",
  row.names = FALSE
)

# ============================================================
# HEATMAP MATRIX
# ============================================================

# ------------------------------------------------------------
# 10. Extract probes
# ------------------------------------------------------------

top_probe_ids <- top_genes$Probe_ID

# Make sure probes exist

top_probe_ids <- top_probe_ids[
  top_probe_ids %in% rownames(expr)
]

if (length(top_probe_ids) < 2) {
  
  stop(
    "ERROR: Not enough valid probes available for heatmap."
  )
  
}

# ------------------------------------------------------------
# 11. Extract expression matrix
# ------------------------------------------------------------

heatmap_expr <- expr[
  top_probe_ids,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 12. Match gene symbols to probes
# ------------------------------------------------------------

probe_to_gene <- top_genes$Gene_Symbol[
  match(
    top_probe_ids,
    top_genes$Probe_ID
  )
]

# Remove missing gene symbols

valid <- !is.na(probe_to_gene) &
  probe_to_gene != ""

heatmap_expr <- heatmap_expr[
  valid,
  ,
  drop = FALSE
]

probe_to_gene <- probe_to_gene[
  valid
]

# ------------------------------------------------------------
# 13. Make row names unique
# ------------------------------------------------------------

rownames(heatmap_expr) <- make.unique(
  probe_to_gene
)

# ------------------------------------------------------------
# 14. Remove rows with NA / Inf
# ------------------------------------------------------------

valid_rows <- apply(
  heatmap_expr,
  1,
  function(x) {
    all(
      is.finite(x)
    )
  }
)

heatmap_expr <- heatmap_expr[
  valid_rows,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# 15. Remove zero-variance genes
# ------------------------------------------------------------

gene_variance <- apply(
  heatmap_expr,
  1,
  var
)

heatmap_expr <- heatmap_expr[
  gene_variance > 0,
  ,
  drop = FALSE
]

cat(
  "\nFinal genes for heatmap:",
  nrow(heatmap_expr),
  "\n"
)

if (nrow(heatmap_expr) < 2) {
  
  stop(
    "ERROR: Fewer than 2 usable genes remain for heatmap."
  )
  
}

# ============================================================
# Z-SCORE
# ============================================================

cat("\nCalculating gene-wise Z-scores...\n")

heatmap_scaled <- t(
  scale(
    t(heatmap_expr)
  )
)

# ------------------------------------------------------------
# Check scaled matrix
# ------------------------------------------------------------

if (
  any(
    !is.finite(
      heatmap_scaled
    )
  )
) {
  
  stop(
    "ERROR: NA/Inf values detected after scaling."
  )
  
}

# ============================================================
# SAMPLE ANNOTATION
# ============================================================

annotation_col <- data.frame(
  Group = group
)

rownames(annotation_col) <- colnames(
  expr
)

# Make sure annotation order matches matrix

annotation_col <- annotation_col[
  colnames(heatmap_scaled),
  ,
  drop = FALSE
]

# ============================================================
# HEATMAP
# ============================================================

cat("\nCreating heatmap...\n")

heatmap_file <- paste0(
  "results/figures/08_top_DEG_heatmap.png"
)

png(
  filename = heatmap_file,
  width = 2800,
  height = 2200,
  res = 300
)

pheatmap(
  mat = heatmap_scaled,
  
  annotation_col = annotation_col,
  
  show_colnames = FALSE,
  
  show_rownames = TRUE,
  
  fontsize_row = 8,
  
  clustering_method = "complete",
  
  border_color = NA,
  
  main = paste0(
    "Top ",
    nrow(heatmap_scaled),
    " Differentially Expressed Genes"
  )
)

dev.off()

# ------------------------------------------------------------
# Verify file
# ------------------------------------------------------------

if (!file.exists(heatmap_file)) {
  
  stop(
    "ERROR: Heatmap file was not created."
  )
  
}

file_size <- file.info(
  heatmap_file
)$size

cat(
  "\nHeatmap file size:",
  file_size,
  "bytes\n"
)

if (file_size == 0) {
  
  stop(
    "ERROR: Heatmap file is 0 bytes."
  )
  
}

# ============================================================
# SAVE MATRIX
# ============================================================

write.csv(
  heatmap_scaled,
  "results/tables/top_DEG_heatmap_matrix.csv"
)

saveRDS(
  heatmap_scaled,
  "data/processed/top_DEG_heatmap_matrix.rds"
)

# ============================================================
# FINAL STATUS
# ============================================================

cat("\n============================================\n")
cat("STEP 6 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat(
  "\nAnnotated probes:",
  nrow(deg_annotated),
  "\n"
)

cat(
  "Significant annotated DEGs:",
  nrow(significant_deg),
  "\n"
)

cat(
  "Final heatmap genes:",
  nrow(heatmap_scaled),
  "\n"
)

cat(
  "Heatmap saved:",
  heatmap_file,
  "\n"
)

cat(
  "File size:",
  file_size,
  "bytes\n"
)
