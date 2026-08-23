# ============================================================
# Lung Cancer Bioinformatics Project
# Step 4: Quality Control & PCA
# Dataset: GSE10072
# ============================================================

# ------------------------------------------------------------
# 1. Clean environment
# ------------------------------------------------------------

rm(list = ls())

# ------------------------------------------------------------
# 2. Load packages
# ------------------------------------------------------------

library(tidyverse)
library(pheatmap)

# ------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------

expr <- readRDS(
  "data/raw/GSE10072_expression.rds"
)

group <- readRDS(
  "data/processed/sample_groups.rds"
)

sample_info <- readRDS(
  "data/processed/sample_information.rds"
)

# ------------------------------------------------------------
# 4. Basic information
# ------------------------------------------------------------

cat("============================================\n")
cat("QUALITY CONTROL & PCA\n")
cat("============================================\n")

cat("\nExpression dimensions:\n")
cat("Probes:", nrow(expr), "\n")
cat("Samples:", ncol(expr), "\n")

cat("\nSample groups:\n")
print(table(group))

# ------------------------------------------------------------
# 5. Verify sample order
# ------------------------------------------------------------

if (!all(
  colnames(expr) == sample_info$Sample_ID
)) {
  
  stop(
    "ERROR: Expression and phenotype sample order do not match."
  )
  
}

cat("\n✓ Sample order verified.\n")

# ------------------------------------------------------------
# 6. Check missing values
# ------------------------------------------------------------

missing_values <- sum(is.na(expr))

cat("\nMissing expression values:", missing_values, "\n")

if (missing_values > 0) {
  
  stop(
    "ERROR: Missing expression values detected."
  )
  
}

cat("✓ No missing expression values.\n")


# ============================================================
# GRAPH 1
# EXPRESSION BOXPLOT
# ============================================================

cat("\nCreating expression boxplot...\n")

png(
  filename = "results/figures/01_expression_boxplot.png",
  width = 2800,
  height = 1800,
  res = 300
)

boxplot(
  expr,
  las = 2,
  cex.axis = 0.35,
  main = "Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Expression Value"
)

dev.off()

cat("✓ Boxplot saved.\n")


# ============================================================
# GRAPH 2
# EXPRESSION DENSITY PLOT
# ============================================================

cat("Creating expression density plot...\n")

# Use a subset of probes for density visualization
# to reduce memory usage.

set.seed(123)

density_probes <- sample(
  seq_len(nrow(expr)),
  size = min(5000, nrow(expr))
)

density_matrix <- expr[
  density_probes,
  ,
  drop = FALSE
]

density_df <- as.data.frame(
  density_matrix
)

density_df$Probe <- rownames(
  density_df
)

density_long <- density_df %>%
  
  pivot_longer(
    cols = -Probe,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  
  left_join(
    sample_info %>%
      select(
        Sample_ID,
        Group
      ),
    by = c(
      "Sample" = "Sample_ID"
    )
  )

density_plot <- ggplot(
  density_long,
  aes(
    x = Expression,
    group = Sample
  )
) +
  
  geom_density(
    alpha = 0.25
  ) +
  
  facet_wrap(
    ~ Group
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Expression Density Distribution",
    subtitle = "Random subset of 5,000 probes",
    x = "Expression Value",
    y = "Density"
  )

ggsave(
  filename = "results/figures/02_expression_density.png",
  plot = density_plot,
  width = 9,
  height = 6,
  dpi = 300
)

cat("✓ Density plot saved.\n")


# ============================================================
# SAMPLE CORRELATION
# ============================================================

cat("\nCalculating sample correlation...\n")

# ------------------------------------------------------------
# Select most variable probes
# ------------------------------------------------------------

probe_variance <- apply(
  expr,
  1,
  var
)

top_n <- min(
  500,
  length(probe_variance)
)

top_variable_probes <- order(
  probe_variance,
  decreasing = TRUE
)[1:top_n]

expr_variable <- expr[
  top_variable_probes,
  ,
  drop = FALSE
]

cat(
  "Using",
  top_n,
  "most variable probes for correlation heatmap.\n"
)

# ------------------------------------------------------------
# Calculate Pearson correlation
# ------------------------------------------------------------

sample_cor <- cor(
  expr_variable,
  method = "pearson"
)

# Verify correlation matrix

if (
  nrow(sample_cor) != ncol(expr) ||
  ncol(sample_cor) != ncol(expr)
) {
  
  stop(
    "ERROR: Correlation matrix dimensions are incorrect."
  )
  
}

# ------------------------------------------------------------
# Create annotation
# ------------------------------------------------------------

annotation_col <- data.frame(
  Group = group
)

rownames(annotation_col) <- colnames(
  expr
)

# Ensure annotation order matches correlation matrix

annotation_col <- annotation_col[
  colnames(sample_cor),
  ,
  drop = FALSE
]

# ============================================================
# GRAPH 3
# SAMPLE CORRELATION HEATMAP
# ============================================================

cat("Creating sample correlation heatmap...\n")

heatmap_file <- "results/figures/03_sample_correlation_heatmap.png"

png(
  filename = heatmap_file,
  width = 2600,
  height = 2400,
  res = 300
)

pheatmap(
  mat = sample_cor,
  
  annotation_col = annotation_col,
  
  annotation_row = annotation_col,
  
  show_colnames = FALSE,
  
  show_rownames = FALSE,
  
  clustering_method = "complete",
  
  clustering_distance_rows = "correlation",
  
  clustering_distance_cols = "correlation",
  
  border_color = NA,
  
  main = "Sample-to-Sample Correlation\nTop 500 Variable Probes"
)

dev.off()

cat("✓ Correlation heatmap saved.\n")


# ------------------------------------------------------------
# Save correlation matrix
# ------------------------------------------------------------

write.csv(
  sample_cor,
  "results/tables/sample_correlation_matrix.csv"
)

# ============================================================
# PCA ANALYSIS
# ============================================================

cat("\nRunning PCA...\n")

# Transpose:
# rows = samples
# columns = probes

expr_for_pca <- t(expr)

# PCA

pca <- prcomp(
  expr_for_pca,
  center = TRUE,
  scale. = FALSE
)

# ------------------------------------------------------------
# Variance explained
# ------------------------------------------------------------

variance_explained <- (
  pca$sdev^2 /
    sum(pca$sdev^2)
) * 100

pc1_variance <- round(
  variance_explained[1],
  2
)

pc2_variance <- round(
  variance_explained[2],
  2
)

pc3_variance <- round(
  variance_explained[3],
  2
)

cat("\nVariance explained:\n")

cat(
  "PC1:",
  pc1_variance,
  "%\n"
)

cat(
  "PC2:",
  pc2_variance,
  "%\n"
)

cat(
  "PC3:",
  pc3_variance,
  "%\n"
)

# ------------------------------------------------------------
# PCA dataframe
# ------------------------------------------------------------

pca_df <- data.frame(
  
  Sample_ID = rownames(
    pca$x
  ),
  
  PC1 = pca$x[, 1],
  
  PC2 = pca$x[, 2],
  
  PC3 = pca$x[, 3],
  
  Group = group
  
)

# ============================================================
# GRAPH 4
# PCA PC1 vs PC2
# ============================================================

cat("\nCreating PCA PC1 vs PC2 plot...\n")

pca_plot <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Group
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Principal Component Analysis",
    subtitle = "GSE10072 Lung Cancer Gene Expression",
    x = paste0(
      "PC1 (",
      pc1_variance,
      "%)"
    ),
    y = paste0(
      "PC2 (",
      pc2_variance,
      "%)"
    ),
    color = "Group"
  )

ggsave(
  filename = "results/figures/04_PCA_PC1_PC2.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("✓ PC1 vs PC2 plot saved.\n")


# ============================================================
# GRAPH 5
# PCA PC1 vs PC3
# ============================================================

cat("Creating PCA PC1 vs PC3 plot...\n")

pca_plot_pc3 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC3,
    color = Group
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  theme_minimal() +
  
  labs(
    title = "PCA: PC1 vs PC3",
    subtitle = "GSE10072 Lung Cancer Gene Expression",
    x = paste0(
      "PC1 (",
      pc1_variance,
      "%)"
    ),
    y = paste0(
      "PC3 (",
      pc3_variance,
      "%)"
    ),
    color = "Group"
  )

ggsave(
  filename = "results/figures/05_PCA_PC1_PC3.png",
  plot = pca_plot_pc3,
  width = 8,
  height = 6,
  dpi = 300
)

cat("✓ PC1 vs PC3 plot saved.\n")


# ============================================================
# SAVE PCA RESULTS
# ============================================================

write.csv(
  pca_df,
  "results/tables/PCA_coordinates.csv",
  row.names = FALSE
)

variance_table <- data.frame(
  PC = paste0(
    "PC",
    seq_along(
      variance_explained
    )
  ),
  
  Variance_Explained = variance_explained
)

write.csv(
  variance_table,
  "results/tables/PCA_variance_explained.csv",
  row.names = FALSE
)

saveRDS(
  pca,
  "data/processed/PCA_object.rds"
)

# ============================================================
# FINAL STATUS
# ============================================================

cat("\n============================================\n")
cat("STEP 4 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat("\nGenerated figures:\n")
cat("1. Expression boxplot\n")
cat("2. Expression density plot\n")
cat("3. Sample correlation heatmap\n")
cat("4. PCA PC1 vs PC2\n")
cat("5. PCA PC1 vs PC3\n")

cat("\nPCA variance explained:\n")

cat(
  "PC1:",
  pc1_variance,
  "%\n"
)

cat(
  "PC2:",
  pc2_variance,
  "%\n"
)

cat(
  "PC3:",
  pc3_variance,
  "%\n"
)
