# ============================================================
# Lung Cancer Bioinformatics Project
# STEP 10: Final Results & Biological Interpretation
# Dataset: GSE10072
# ============================================================

rm(list = ls())

# ============================================================
# 1. LOAD PACKAGES
# ============================================================

library(tidyverse)

# ============================================================
# 2. LOAD RESULTS
# ============================================================

deg <- read.csv(
  "results/tables/annotated_differential_expression_results.csv",
  stringsAsFactors = FALSE
)

top_deg <- read.csv(
  "results/tables/top_30_annotated_DEGs.csv",
  stringsAsFactors = FALSE
)

go_results <- read.csv(
  "results/tables/GO_Biological_Process_enrichment.csv",
  stringsAsFactors = FALSE
)

kegg_results <- read.csv(
  "results/tables/KEGG_pathway_enrichment.csv",
  stringsAsFactors = FALSE
)

model_results <- read.csv(
  "results/tables/cross_validated_model_comparison.csv",
  stringsAsFactors = FALSE
)

fold_results <- read.csv(
  "results/tables/fold_level_model_results.csv",
  stringsAsFactors = FALSE
)

# ============================================================
# 3. BASIC DEG SUMMARY
# ============================================================

cat("============================================\n")
cat("STEP 10: FINAL RESULTS\n")
cat("============================================\n")

total_probes <- nrow(deg)

significant_deg <- deg %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) >= 1
  )

upregulated <- significant_deg %>%
  filter(
    logFC >= 1
  )

downregulated <- significant_deg %>%
  filter(
    logFC <= -1
  )

cat("\nDEG SUMMARY\n")
cat("--------------------------------------------\n")

cat(
  "Total probes:",
  total_probes,
  "\n"
)

cat(
  "Significant DEGs:",
  nrow(significant_deg),
  "\n"
)

cat(
  "Upregulated:",
  nrow(upregulated),
  "\n"
)

cat(
  "Downregulated:",
  nrow(downregulated),
  "\n"
)

# ============================================================
# 4. TOP UPREGULATED GENES
# ============================================================

top_up <- upregulated %>%
  
  arrange(
    adj.P.Val
  ) %>%
  
  select(
    Probe_ID,
    Gene_Symbol,
    logFC,
    P.Value,
    adj.P.Val
  ) %>%
  
  head(20)

write.csv(
  top_up,
  "results/tables/final_top_20_upregulated_genes.csv",
  row.names = FALSE
)

# ============================================================
# 5. TOP DOWNREGULATED GENES
# ============================================================

top_down <- downregulated %>%
  
  arrange(
    adj.P.Val
  ) %>%
  
  select(
    Probe_ID,
    Gene_Symbol,
    logFC,
    P.Value,
    adj.P.Val
  ) %>%
  
  head(20)

write.csv(
  top_down,
  "results/tables/final_top_20_downregulated_genes.csv",
  row.names = FALSE
)

# ============================================================
# 6. TOP GO TERMS
# ============================================================

if (
  nrow(go_results) > 0
) {
  
  top_go <- go_results %>%
    
    arrange(
      p.adjust
    ) %>%
    
    head(15)
  
  write.csv(
    top_go,
    "results/tables/final_top_GO_terms.csv",
    row.names = FALSE
  )
  
} else {
  
  top_go <- data.frame()
  
}

# ============================================================
# 7. TOP KEGG PATHWAYS
# ============================================================

if (
  nrow(kegg_results) > 0
) {
  
  top_kegg <- kegg_results %>%
    
    arrange(
      p.adjust
    ) %>%
    
    head(15)
  
  write.csv(
    top_kegg,
    "results/tables/final_top_KEGG_pathways.csv",
    row.names = FALSE
  )
  
} else {
  
  top_kegg <- data.frame()
  
}

# ============================================================
# 8. BEST ML MODEL
# ============================================================

best_model_index <- which.max(
  model_results$Mean_AUC
)

best_model <- model_results$Model[
  best_model_index
]

best_auc <- model_results$Mean_AUC[
  best_model_index
]

cat("\nMACHINE LEARNING RESULTS\n")
cat("--------------------------------------------\n")

print(
  model_results
)

cat(
  "\nBest model:",
  best_model,
  "\n"
)

cat(
  "Mean 5-fold ROC-AUC:",
  round(
    best_auc,
    3
  ),
  "\n"
)

# ============================================================
# 9. CREATE FINAL SUMMARY TABLE
# ============================================================

final_summary <- data.frame(
  
  Analysis = c(
    "Total probes analyzed",
    "Significant DEGs",
    "Upregulated DEGs",
    "Downregulated DEGs",
    "GO enriched terms",
    "KEGG enriched pathways",
    "Best ML model",
    "Best mean 5-fold ROC-AUC"
  ),
  
  Result = c(
    total_probes,
    nrow(significant_deg),
    nrow(upregulated),
    nrow(downregulated),
    nrow(go_results),
    nrow(kegg_results),
    best_model,
    round(
      best_auc,
      3
    )
  )
  
)

write.csv(
  final_summary,
  "results/tables/final_project_summary.csv",
  row.names = FALSE
)

# ============================================================
# 10. DEG BARPLOT
# ============================================================

deg_counts <- data.frame(
  
  Category = c(
    "Upregulated",
    "Downregulated"
  ),
  
  Count = c(
    nrow(upregulated),
    nrow(downregulated)
  )
  
)

deg_plot <- ggplot(
  deg_counts,
  aes(
    x = Category,
    y = Count
  )
) +
  
  geom_col() +
  
  theme_minimal() +
  
  labs(
    title = "Differentially Expressed Genes",
    x = "Category",
    y = "Number of Genes"
  )

ggsave(
  "results/figures/14_DEG_summary.png",
  deg_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 11. TOP GO PLOT
# ============================================================

if (
  nrow(top_go) > 0
) {
  
  go_plot_data <- top_go %>%
    
    arrange(
      desc(
        -log10(p.adjust)
      )
    ) %>%
    
    head(10)
  
  go_plot_data$Description <- factor(
    go_plot_data$Description,
    levels = rev(
      go_plot_data$Description
    )
  )
  
  go_plot <- ggplot(
    go_plot_data,
    aes(
      x = -log10(p.adjust),
      y = Description
    )
  ) +
    
    geom_col() +
    
    theme_minimal() +
    
    labs(
      title = "Top Enriched GO Biological Processes",
      x = "-Log10 Adjusted P-value",
      y = "Biological Process"
    )
  
  ggsave(
    "results/figures/15_top_GO_terms.png",
    go_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
  
}

# ============================================================
# 12. TOP KEGG PLOT
# ============================================================

if (
  nrow(top_kegg) > 0
) {
  
  kegg_plot_data <- top_kegg %>%
    
    arrange(
      desc(
        -log10(p.adjust)
      )
    ) %>%
    
    head(10)
  
  kegg_plot_data$Description <- factor(
    kegg_plot_data$Description,
    levels = rev(
      kegg_plot_data$Description
    )
  )
  
  kegg_plot <- ggplot(
    kegg_plot_data,
    aes(
      x = -log10(p.adjust),
      y = Description
    )
  ) +
    
    geom_col() +
    
    theme_minimal() +
    
    labs(
      title = "Top Enriched KEGG Pathways",
      x = "-Log10 Adjusted P-value",
      y = "Pathway"
    )
  
  ggsave(
    "results/figures/16_top_KEGG_pathways.png",
    kegg_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
  
}

# ============================================================
# 13. PRINT TOP GENES
# ============================================================

cat("\n")
cat("============================================\n")
cat("TOP UPREGULATED GENES\n")
cat("============================================\n")

print(
  top_up
)

cat("\n")
cat("============================================\n")
cat("TOP DOWNREGULATED GENES\n")
cat("============================================\n")

print(
  top_down
)

# ============================================================
# 14. PRINT TOP GO TERMS
# ============================================================

if (
  nrow(top_go) > 0
) {
  
  cat("\n")
  cat("============================================\n")
  cat("TOP GO BIOLOGICAL PROCESSES\n")
  cat("============================================\n")
  
  print(
    top_go[
      ,
      c(
        "Description",
        "GeneRatio",
        "p.adjust"
      )
    ]
  )
  
}

# ============================================================
# 15. PRINT TOP KEGG PATHWAYS
# ============================================================

if (
  nrow(top_kegg) > 0
) {
  
  cat("\n")
  cat("============================================\n")
  cat("TOP KEGG PATHWAYS\n")
  cat("============================================\n")
  
  print(
    top_kegg[
      ,
      c(
        "Description",
        "GeneRatio",
        "p.adjust"
      )
    ]
  )
  
}

# ============================================================
# 16. FINAL PROJECT STATUS
# ============================================================

cat("\n")
cat("============================================\n")
cat("STEP 10 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat(
  "\nTotal probes:",
  total_probes,
  "\n"
)

cat(
  "Significant DEGs:",
  nrow(significant_deg),
  "\n"
)

cat(
  "Upregulated:",
  nrow(upregulated),
  "\n"
)

cat(
  "Downregulated:",
  nrow(downregulated),
  "\n"
)

cat(
  "GO terms:",
  nrow(go_results),
  "\n"
)

cat(
  "KEGG pathways:",
  nrow(kegg_results),
  "\n"
)

cat(
  "Best ML model:",
  best_model,
  "\n"
)

cat(
  "Best mean 5-fold AUC:",
  round(
    best_auc,
    3
  ),
  "\n"
)
