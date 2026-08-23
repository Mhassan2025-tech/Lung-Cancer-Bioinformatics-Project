# ============================================================
# Lung Cancer Bioinformatics Project
# Step 3: Phenotype Processing
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

# ------------------------------------------------------------
# 3. Load saved data from Step 2
# ------------------------------------------------------------

expr <- readRDS(
  "data/raw/GSE10072_expression.rds"
)

pheno <- readRDS(
  "data/raw/GSE10072_phenotype.rds"
)

# ------------------------------------------------------------
# 4. Basic verification
# ------------------------------------------------------------

cat("============================================\n")
cat("PHENOTYPE PROCESSING\n")
cat("============================================\n")

cat("\nExpression dimensions:\n")
print(dim(expr))

cat("\nPhenotype dimensions:\n")
print(dim(pheno))

# ------------------------------------------------------------
# 5. Verify sample names match
# ------------------------------------------------------------

expression_samples <- colnames(expr)

phenotype_samples <- rownames(pheno)

if (!all(expression_samples %in% phenotype_samples)) {
  
  stop(
    "ERROR: Expression and phenotype sample IDs do not match."
  )
  
}

cat("\n✓ Expression and phenotype samples match.\n")

# ------------------------------------------------------------
# 6. Identify sample title column
# ------------------------------------------------------------

if (!"title" %in% colnames(pheno)) {
  
  stop(
    "ERROR: 'title' column was not found in phenotype data."
  )
  
}

sample_titles <- as.character(
  pheno$title
)

# ------------------------------------------------------------
# 7. Create Tumor / Normal labels
# ------------------------------------------------------------

group <- case_when(
  
  grepl(
    "lung tumor",
    sample_titles,
    ignore.case = TRUE
  ) ~ "Tumor",
  
  grepl(
    "normal lung",
    sample_titles,
    ignore.case = TRUE
  ) ~ "Normal",
  
  TRUE ~ NA_character_
  
)

# ------------------------------------------------------------
# 8. Check for unidentified samples
# ------------------------------------------------------------

if (any(is.na(group))) {
  
  cat("\nUnidentified samples:\n")
  
  print(
    data.frame(
      Sample = rownames(pheno)[is.na(group)],
      Title = sample_titles[is.na(group)]
    )
  )
  
  stop(
    "ERROR: Some samples could not be classified."
  )
  
}

# ------------------------------------------------------------
# 9. Convert group to factor
# ------------------------------------------------------------

group <- factor(
  group,
  levels = c(
    "Normal",
    "Tumor"
  )
)

# ------------------------------------------------------------
# 10. Count groups
# ------------------------------------------------------------

cat("\n============================================\n")
cat("SAMPLE GROUPS\n")
cat("============================================\n")

group_counts <- table(group)

print(group_counts)

# ------------------------------------------------------------
# 11. Verify expected sample counts
# ------------------------------------------------------------

if (
  as.integer(group_counts["Normal"]) != 49 ||
  as.integer(group_counts["Tumor"]) != 58
) {
  
  stop(
    paste(
      "ERROR: Unexpected group counts.",
      "Expected: 49 Normal and 58 Tumor."
    )
  )
  
}

cat("\n✓ Expected 49 Normal samples found.\n")
cat("✓ Expected 58 Tumor samples found.\n")

# ------------------------------------------------------------
# 12. Create phenotype table
# ------------------------------------------------------------

sample_info <- data.frame(
  
  Sample_ID = colnames(expr),
  
  Title = pheno[
    colnames(expr),
    "title"
  ],
  
  Group = group
  
)

# ------------------------------------------------------------
# 13. Display sample information
# ------------------------------------------------------------

cat("\nFirst 10 samples:\n\n")

print(
  head(
    sample_info,
    10
  )
)

# ------------------------------------------------------------
# 14. Create group summary
# ------------------------------------------------------------

group_summary <- sample_info %>%
  
  count(
    Group
  )

cat("\nGroup summary:\n")

print(
  group_summary
)

# ------------------------------------------------------------
# 15. Save phenotype information
# ------------------------------------------------------------

write.csv(
  sample_info,
  "data/processed/sample_information.csv",
  row.names = FALSE
)

# Save group vector
saveRDS(
  group,
  "data/processed/sample_groups.rds"
)

# Save processed phenotype table
saveRDS(
  sample_info,
  "data/processed/sample_information.rds"
)

# ------------------------------------------------------------
# 16. Create simple group plot
# ------------------------------------------------------------

group_plot <- ggplot(
  group_summary,
  aes(
    x = Group,
    y = n
  )
) +
  
  geom_col() +
  
  theme_minimal() +
  
  labs(
    title = "GSE10072 Sample Groups",
    x = "Sample Group",
    y = "Number of Samples"
  )

print(group_plot)

# ------------------------------------------------------------
# 17. Save plot
# ------------------------------------------------------------

ggsave(
  "results/figures/sample_groups.png",
  group_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 18. Final status
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 3 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat("\nNormal samples:", sum(group == "Normal"), "\n")
cat("Tumor samples:", sum(group == "Tumor"), "\n")

