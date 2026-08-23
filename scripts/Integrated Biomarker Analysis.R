# ============================================================
# Lung Cancer Bioinformatics Project
# STEP 11: Integrated Biomarker Analysis
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

rf_importance <- read.csv(
  "results/tables/random_forest_feature_importance.csv",
  stringsAsFactors = FALSE
)

# ============================================================
# 3. BASIC CHECK
# ============================================================

cat("============================================\n")
cat("STEP 11: INTEGRATED BIOMARKER ANALYSIS\n")
cat("============================================\n")

cat(
  "\nDEG records:",
  nrow(deg),
  "\n"
)

cat(
  "Random Forest features:",
  nrow(rf_importance),
  "\n"
)

# ============================================================
# 4. IDENTIFY PROBE ID COLUMN
# ============================================================

if (!"Probe_ID" %in% colnames(deg)) {
  
  stop(
    "ERROR: Probe_ID column not found in DEG table."
  )
  
}

if (!"Feature" %in% colnames(rf_importance)) {
  
  stop(
    "ERROR: Feature column not found in Random Forest table."
  )
  
}

# ============================================================
# 5. MAKE RF PROBE IDs COMPARABLE
# ============================================================

# Step 7 used make.names(), so convert them back
# approximately to original probe IDs.

rf_importance$Probe_ID <- rf_importance$Feature

# Match safe R names to original DEG probe IDs

safe_deg_names <- make.names(
  deg$Probe_ID,
  unique = TRUE
)

rf_importance$Probe_ID <- deg$Probe_ID[
  match(
    rf_importance$Feature,
    safe_deg_names
  )
]

# Remove unmatched features

rf_importance <- rf_importance %>%
  
  filter(
    !is.na(Probe_ID)
  )

cat(
  "\nMatched RF features:",
  nrow(rf_importance),
  "\n"
)

# ============================================================
# 6. MERGE DEG + RANDOM FOREST RESULTS
# ============================================================

integrated <- deg %>%
  
  select(
    Probe_ID,
    Gene_Symbol,
    logFC,
    P.Value,
    adj.P.Val,
    Regulation
  ) %>%
  
  inner_join(
    rf_importance %>%
      select(
        Probe_ID,
        MeanDecreaseGini
      ),
    by = "Probe_ID"
  )

cat(
  "Integrated features:",
  nrow(integrated),
  "\n"
)

# ============================================================
# 7. IDENTIFY SIGNIFICANT ML + DEG FEATURES
# ============================================================

integrated_significant <- integrated %>%
  
  filter(
    adj.P.Val < 0.05,
    abs(logFC) >= 1,
    !is.na(Gene_Symbol),
    Gene_Symbol != ""
  ) %>%
  
  arrange(
    desc(
      MeanDecreaseGini
    )
  )

cat(
  "\nSignificant DEG + ML features:",
  nrow(integrated_significant),
  "\n"
)

# ============================================================
# 8. REMOVE DUPLICATE GENE SYMBOLS
# ============================================================

integrated_unique <- integrated_significant %>%
  
  group_by(
    Gene_Symbol
  ) %>%
  
  slice_max(
    order_by = MeanDecreaseGini,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup()

cat(
  "Unique candidate genes:",
  nrow(integrated_unique),
  "\n"
)

# ============================================================
# 9. CREATE RANKS
# ============================================================

integrated_unique <- integrated_unique %>%
  
  mutate(
    
    DEG_Rank = rank(
      -abs(logFC),
      ties.method = "min"
    ),
    
    ML_Rank = rank(
      -MeanDecreaseGini,
      ties.method = "min"
    )
    
  ) %>%
  
  arrange(
    DEG_Rank + ML_Rank
  )

# ============================================================
# 10. TOP 20 INTEGRATED CANDIDATES
# ============================================================

top_candidates <- integrated_unique %>%
  
  select(
    Gene_Symbol,
    Probe_ID,
    logFC,
    adj.P.Val,
    Regulation,
    MeanDecreaseGini,
    DEG_Rank,
    ML_Rank
  ) %>%
  
  head(20)

cat("\n")
cat("============================================\n")
cat("TOP INTEGRATED CANDIDATE GENES\n")
cat("============================================\n")

print(
  top_candidates
)

# ============================================================
# 11. SAVE RESULTS
# ============================================================

write.csv(
  integrated,
  "results/tables/integrated_DEG_ML_results.csv",
  row.names = FALSE
)

write.csv(
  integrated_unique,
  "results/tables/integrated_candidate_genes.csv",
  row.names = FALSE
)

write.csv(
  top_candidates,
  "results/tables/top_20_integrated_biomarker_candidates.csv",
  row.names = FALSE
)

# ============================================================
# 12. INTEGRATED SCATTER PLOT
# ============================================================

plot_data <- integrated_unique %>%
  
  mutate(
    NegLog10FDR =
      -log10(
        pmax(
          adj.P.Val,
          .Machine$double.xmin
        )
      )
  )

candidate_plot <- ggplot(
  plot_data,
  aes(
    x = logFC,
    y = MeanDecreaseGini
  )
) +
  
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Integrated DEG + Machine Learning Features",
    subtitle = "Significant DEGs ranked by Random Forest importance",
    x = "Log2 Fold Change",
    y = "Random Forest Mean Decrease Gini"
  )

ggsave(
  "results/figures/17_integrated_DEG_ML_features.png",
  candidate_plot,
  width = 9,
  height = 7,
  dpi = 300
)

cat(
  "\n✓ Integrated feature plot saved.\n"
)

# ============================================================
# 13. TOP 20 BIOMARKER PLOT
# ============================================================

top_plot_data <- top_candidates %>%
  
  arrange(
    MeanDecreaseGini
  ) %>%
  
  mutate(
    Gene_Symbol = factor(
      Gene_Symbol,
      levels = Gene_Symbol
    )
  )

biomarker_plot <- ggplot(
  top_plot_data,
  aes(
    x = Gene_Symbol,
    y = MeanDecreaseGini
  )
) +
  
  geom_col() +
  
  coord_flip() +
  
  theme_minimal() +
  
  labs(
    title = "Top Integrated Biomarker Candidates",
    x = "Gene",
    y = "Random Forest Importance"
  )

ggsave(
  "results/figures/18_top_integrated_biomarkers.png",
  biomarker_plot,
  width = 9,
  height = 8,
  dpi = 300
)

cat(
  "✓ Biomarker plot saved.\n"
)

# ============================================================
# 14. FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================\n")
cat("STEP 11 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat(
  "\nIntegrated features:",
  nrow(integrated),
  "\n"
)

cat(
  "Significant DEG + ML features:",
  nrow(integrated_significant),
  "\n"
)

cat(
  "Unique candidate genes:",
  nrow(integrated_unique),
  "\n"
)

cat(
  "Top candidates saved:",
  nrow(top_candidates),
  "\n"
)

cat("\nGenerated figures:\n")

cat(
  "17_integrated_DEG_ML_features.png\n"
)

cat(
  "18_top_integrated_biomarkers.png\n"
)

cat("\nGenerated tables:\n")

cat(
  "integrated_DEG_ML_results.csv\n"
)

cat(
  "integrated_candidate_genes.csv\n"
)

cat(
  "top_20_integrated_biomarker_candidates.csv\n"
)
