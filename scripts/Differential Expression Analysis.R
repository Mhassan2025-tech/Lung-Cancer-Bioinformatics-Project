# ============================================================
# STEP 5: Differential Expression Analysis
# Dataset: GSE10072
# Comparison: Tumor vs Normal
# Method: limma
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# Load packages
# ------------------------------------------------------------

library(limma)
library(tidyverse)

# ------------------------------------------------------------
# Load data
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
# Verify data
# ------------------------------------------------------------

cat("============================================\n")
cat("DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat("============================================\n")

cat("\nExpression dimensions:\n")
cat("Probes:", nrow(expr), "\n")
cat("Samples:", ncol(expr), "\n")

cat("\nGroups:\n")
print(table(group))

# Check sample order

if (!all(colnames(expr) == sample_info$Sample_ID)) {
  stop("ERROR: Sample order does not match.")
}

cat("\n✓ Sample order verified.\n")

# ------------------------------------------------------------
# Create group factor
# ------------------------------------------------------------

group <- factor(
  group,
  levels = c("Normal", "Tumor")
)

# ------------------------------------------------------------
# Design matrix
# ------------------------------------------------------------

cat("\nCreating design matrix...\n")

design <- model.matrix(
  ~ group
)

colnames(design) <- c(
  "Intercept",
  "Tumor_vs_Normal"
)

print(head(design))

# ------------------------------------------------------------
# Fit limma model
# ------------------------------------------------------------

cat("\nRunning limma analysis...\n")

fit <- lmFit(
  expr,
  design
)

fit <- eBayes(
  fit
)

cat("✓ limma analysis completed.\n")

# ------------------------------------------------------------
# Extract results
# ------------------------------------------------------------

deg <- topTable(
  fit,
  coef = "Tumor_vs_Normal",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Add probe IDs

deg <- deg %>%
  rownames_to_column("Probe_ID")

# ------------------------------------------------------------
# Classify genes
# ------------------------------------------------------------

deg <- deg %>%
  mutate(
    Regulation = case_when(
      
      adj.P.Val < 0.05 &
        logFC >= 1 ~ "Upregulated",
      
      adj.P.Val < 0.05 &
        logFC <= -1 ~ "Downregulated",
      
      TRUE ~ "Not Significant"
    )
  )

# ------------------------------------------------------------
# DEG summary
# ------------------------------------------------------------

total_probes <- nrow(deg)

significant_degs <- sum(
  deg$adj.P.Val < 0.05
)

upregulated <- sum(
  deg$Regulation == "Upregulated"
)

downregulated <- sum(
  deg$Regulation == "Downregulated"
)

cat("\n============================================\n")
cat("DEG SUMMARY\n")
cat("============================================\n")

cat(
  "\nTotal probes:",
  total_probes,
  "\n"
)

cat(
  "Significant DEGs:",
  significant_degs,
  "\n"
)

cat(
  "Upregulated:",
  upregulated,
  "\n"
)

cat(
  "Downregulated:",
  downregulated,
  "\n"
)

# ------------------------------------------------------------
# Top upregulated
# ------------------------------------------------------------

upregulated_genes <- deg %>%
  filter(
    adj.P.Val < 0.05,
    logFC >= 1
  ) %>%
  arrange(
    adj.P.Val
  )

# ------------------------------------------------------------
# Top downregulated
# ------------------------------------------------------------

downregulated_genes <- deg %>%
  filter(
    adj.P.Val < 0.05,
    logFC <= -1
  ) %>%
  arrange(
    adj.P.Val
  )

# ------------------------------------------------------------
# Save tables
# ------------------------------------------------------------

write.csv(
  deg,
  "results/tables/all_differential_expression_results.csv",
  row.names = FALSE
)

write.csv(
  deg %>% filter(adj.P.Val < 0.05),
  "results/tables/significant_DEGs_FDR_0.05.csv",
  row.names = FALSE
)

write.csv(
  upregulated_genes,
  "results/tables/upregulated_DEGs.csv",
  row.names = FALSE
)

write.csv(
  downregulated_genes,
  "results/tables/downregulated_DEGs.csv",
  row.names = FALSE
)

# Top 20

top_20_DEGs <- deg %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  head(20)

write.csv(
  top_20_DEGs,
  "results/tables/top_20_DEGs.csv",
  row.names = FALSE
)

# Summary

deg_summary <- data.frame(
  Metric = c(
    "Total probes",
    "Significant DEGs",
    "Upregulated",
    "Downregulated"
  ),
  
  Count = c(
    total_probes,
    significant_degs,
    upregulated,
    downregulated
  )
)

write.csv(
  deg_summary,
  "results/tables/DEG_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Volcano Plot
# ------------------------------------------------------------

cat("\nCreating volcano plot...\n")

volcano_df <- deg %>%
  mutate(
    neg_log10_FDR =
      -log10(
        pmax(
          adj.P.Val,
          .Machine$double.xmin
        )
      )
  )

volcano_plot <- ggplot(
  volcano_df,
  aes(
    x = logFC,
    y = neg_log10_FDR,
    color = Regulation
  )
) +
  
  geom_point(
    alpha = 0.65,
    size = 1.5
  ) +
  
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Volcano Plot: Tumor vs Normal",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Regulation"
  )

ggsave(
  "results/figures/06_volcano_plot.png",
  volcano_plot,
  width = 9,
  height = 7,
  dpi = 300
)

cat("✓ Volcano plot saved.\n")

# ------------------------------------------------------------
# MA Plot
# ------------------------------------------------------------

cat("Creating MA plot...\n")

png(
  "results/figures/07_MA_plot.png",
  width = 2400,
  height = 1800,
  res = 300
)

plotMA(
  fit,
  coef = "Tumor_vs_Normal",
  main = "MA Plot: Tumor vs Normal",
  ylim = c(-5, 5)
)

abline(
  h = c(-1, 1),
  lty = 2
)

dev.off()

cat("✓ MA plot saved.\n")

# ------------------------------------------------------------
# Save limma objects
# ------------------------------------------------------------

saveRDS(
  fit,
  "data/processed/limma_fit.rds"
)

saveRDS(
  design,
  "data/processed/limma_design.rds"
)

# ------------------------------------------------------------
# Final message
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 5 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat(
  "\nTotal probes:",
  total_probes,
  "\n"
)

cat(
  "Significant DEGs:",
  significant_degs,
  "\n"
)

cat(
  "Upregulated:",
  upregulated,
  "\n"
)

cat(
  "Downregulated:",
  downregulated,
  "\n"
)
