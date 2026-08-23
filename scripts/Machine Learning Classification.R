# ============================================================
# Lung Cancer Bioinformatics Project
# Step 7: Machine Learning Classification
# Dataset: GSE10072
# Classification: Tumor vs Normal
# Model: Random Forest
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(tidyverse)
library(caret)
library(randomForest)
library(pROC)
library(limma)

# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------

expr <- readRDS(
  "data/raw/GSE10072_expression.rds"
)

group <- readRDS(
  "data/processed/sample_groups.rds"
)

# ------------------------------------------------------------
# 3. Basic verification
# ------------------------------------------------------------

cat("============================================\n")
cat("MACHINE LEARNING CLASSIFICATION\n")
cat("============================================\n")

cat("\nExpression dimensions:\n")
cat("Probes:", nrow(expr), "\n")
cat("Samples:", ncol(expr), "\n")

cat("\nGroups:\n")
print(table(group))

# ------------------------------------------------------------
# 4. Make sure group is correctly ordered
# ------------------------------------------------------------

group <- factor(
  group,
  levels = c(
    "Normal",
    "Tumor"
  )
)

if (any(is.na(group))) {
  stop("ERROR: Missing group labels detected.")
}

# ============================================================
# 5. TRAIN / TEST SPLIT
# ============================================================

set.seed(123)

train_index <- createDataPartition(
  group,
  p = 0.80,
  list = FALSE
)

train_expr <- expr[
  ,
  train_index,
  drop = FALSE
]

test_expr <- expr[
  ,
  -train_index,
  drop = FALSE
]

train_group <- group[
  train_index
]

test_group <- group[
  -train_index
]

cat("\n============================================\n")
cat("TRAIN / TEST SPLIT\n")
cat("============================================\n")

cat(
  "\nTraining samples:",
  length(train_group),
  "\n"
)

cat(
  "Testing samples:",
  length(test_group),
  "\n"
)

cat("\nTraining groups:\n")
print(table(train_group))

cat("\nTesting groups:\n")
print(table(test_group))

# ============================================================
# 6. FEATURE SELECTION
# IMPORTANT:
# Feature selection is performed ONLY on training data.
# ============================================================

cat("\nPerforming training-only feature selection...\n")

train_design <- model.matrix(
  ~ train_group
)

colnames(train_design) <- c(
  "Intercept",
  "Tumor_vs_Normal"
)

# limma model

train_fit <- lmFit(
  train_expr,
  train_design
)

train_fit <- eBayes(
  train_fit
)

train_deg <- topTable(
  train_fit,
  coef = "Tumor_vs_Normal",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# ------------------------------------------------------------
# 7. Select top 50 features
# ------------------------------------------------------------

top_feature_number <- min(
  50,
  nrow(train_deg)
)

selected_features <- rownames(
  train_deg
)[
  1:top_feature_number
]

# Make sure selected probes exist

selected_features <- selected_features[
  selected_features %in% rownames(expr)
]

cat(
  "\nSelected features:",
  length(selected_features),
  "\n"
)

if (length(selected_features) < 2) {
  stop(
    "ERROR: Not enough features available for machine learning."
  )
}

# ============================================================
# 8. CREATE TRAINING FEATURE MATRIX
# ============================================================

train_x <- t(
  train_expr[
    selected_features,
    ,
    drop = FALSE
  ]
)

# Create testing feature matrix

test_x <- t(
  test_expr[
    selected_features,
    ,
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# 9. Make feature names safe for R
# ------------------------------------------------------------

safe_feature_names <- make.names(
  colnames(train_x),
  unique = TRUE
)

colnames(train_x) <- safe_feature_names
colnames(test_x) <- safe_feature_names

# ------------------------------------------------------------
# 10. Remove problematic features
# ------------------------------------------------------------

feature_sd <- apply(
  train_x,
  2,
  sd
)

valid_features <- names(
  feature_sd[
    is.finite(feature_sd) &
      feature_sd > 0
  ]
)

train_x <- train_x[
  ,
  valid_features,
  drop = FALSE
]

test_x <- test_x[
  ,
  valid_features,
  drop = FALSE
]

cat(
  "Final ML features:",
  ncol(train_x),
  "\n"
)

# ============================================================
# 11. CHECK FOR MISSING VALUES
# ============================================================

if (
  any(
    !is.finite(
      train_x
    )
  )
) {
  stop(
    "ERROR: Invalid values detected in training data."
  )
}

if (
  any(
    !is.finite(
      test_x
    )
  )
) {
  stop(
    "ERROR: Invalid values detected in testing data."
  )
}

cat("✓ Training and testing matrices verified.\n")

# ============================================================
# 12. RANDOM FOREST
# ============================================================

cat("\n============================================\n")
cat("RANDOM FOREST\n")
cat("============================================\n")

set.seed(123)

rf_model <- randomForest(
  x = train_x,
  y = train_group,
  ntree = 500,
  importance = TRUE
)

cat("\n✓ Random Forest trained successfully.\n")

print(rf_model)

# ============================================================
# 13. PREDICTION
# ============================================================

cat("\nGenerating predictions...\n")

rf_prediction <- predict(
  rf_model,
  newdata = test_x
)

rf_probability <- predict(
  rf_model,
  newdata = test_x,
  type = "prob"
)[,
  "Tumor"
]

# ============================================================
# 14. CONFUSION MATRIX
# ============================================================

cat("\n============================================\n")
cat("CONFUSION MATRIX\n")
cat("============================================\n")

conf_matrix <- confusionMatrix(
  rf_prediction,
  test_group,
  positive = "Tumor"
)

print(conf_matrix)

# ------------------------------------------------------------
# Performance metrics
# ------------------------------------------------------------

accuracy <- as.numeric(
  conf_matrix$overall[
    "Accuracy"
  ]
)

sensitivity <- as.numeric(
  conf_matrix$byClass[
    "Sensitivity"
  ]
)

specificity <- as.numeric(
  conf_matrix$byClass[
    "Specificity"
  ]
)

# ============================================================
# 15. ROC-AUC
# ============================================================

cat("\nCalculating ROC-AUC...\n")

roc_object <- roc(
  response = test_group,
  predictor = rf_probability,
  levels = c(
    "Normal",
    "Tumor"
  ),
  direction = "<",
  quiet = TRUE
)

auc_value <- as.numeric(
  auc(
    roc_object
  )
)

cat(
  "\nROC-AUC:",
  round(
    auc_value,
    4
  ),
  "\n"
)

# ============================================================
# 16. ROC PLOT
# ============================================================

roc_df <- data.frame(
  False_Positive_Rate =
    1 - roc_object$specificities,
  
  True_Positive_Rate =
    roc_object$sensitivities
)

roc_plot <- ggplot(
  roc_df,
  aes(
    x = False_Positive_Rate,
    y = True_Positive_Rate
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_abline(
    linetype = "dashed"
  ) +
  
  theme_minimal() +
  
  labs(
    title = "Random Forest ROC Curve",
    subtitle = paste0(
      "GSE10072 | AUC = ",
      round(
        auc_value,
        3
      )
    ),
    x = "False Positive Rate",
    y = "True Positive Rate"
  )

ggsave(
  "results/figures/09_random_forest_ROC.png",
  roc_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("✓ ROC curve saved.\n")

# ============================================================
# 17. FEATURE IMPORTANCE
# ============================================================

importance_matrix <- importance(
  rf_model
)

importance_df <- data.frame(
  Feature = rownames(
    importance_matrix
  ),
  
  MeanDecreaseGini =
    importance_matrix[
      ,
      "MeanDecreaseGini"
    ]
)

importance_df <- importance_df %>%
  
  arrange(
    desc(
      MeanDecreaseGini
    )
  )

# Save complete importance

write.csv(
  importance_df,
  "results/tables/random_forest_feature_importance.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Top 20 features
# ------------------------------------------------------------

top_importance <- importance_df %>%
  head(20)

# ============================================================
# 18. FEATURE IMPORTANCE PLOT
# ============================================================

importance_plot <- ggplot(
  top_importance,
  aes(
    x = reorder(
      Feature,
      MeanDecreaseGini
    ),
    y = MeanDecreaseGini
  )
) +
  
  geom_col() +
  
  coord_flip() +
  
  theme_minimal() +
  
  labs(
    title = "Top Random Forest Features",
    x = "Probe",
    y = "Mean Decrease Gini"
  )

ggsave(
  "results/figures/10_random_forest_feature_importance.png",
  importance_plot,
  width = 9,
  height = 7,
  dpi = 300
)

cat(
  "✓ Feature importance plot saved.\n"
)

# ============================================================
# 19. SAVE PREDICTIONS
# ============================================================

prediction_results <- data.frame(
  
  Sample = rownames(
    test_x
  ),
  
  Actual = test_group,
  
  Predicted = rf_prediction,
  
  Tumor_Probability = rf_probability
  
)

write.csv(
  prediction_results,
  "results/tables/random_forest_predictions.csv",
  row.names = FALSE
)

# ============================================================
# 20. SAVE PERFORMANCE
# ============================================================

performance <- data.frame(
  
  Metric = c(
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "ROC_AUC"
  ),
  
  Value = c(
    accuracy,
    sensitivity,
    specificity,
    auc_value
  )
  
)

write.csv(
  performance,
  "results/tables/random_forest_performance.csv",
  row.names = FALSE
)

# ============================================================
# 21. SAVE MODEL
# ============================================================

saveRDS(
  rf_model,
  "data/processed/random_forest_model.rds"
)

saveRDS(
  selected_features,
  "data/processed/selected_ML_features.rds"
)

# ============================================================
# FINAL STATUS
# ============================================================

cat("\n============================================\n")
cat("STEP 7 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat(
  "\nSelected features:",
  length(valid_features),
  "\n"
)

cat(
  "Accuracy:",
  round(
    accuracy,
    3
  ),
  "\n"
)

cat(
  "Sensitivity:",
  round(
    sensitivity,
    3
  ),
  "\n"
)

cat(
  "Specificity:",
  round(
    specificity,
    3
  ),
  "\n"
)

cat(
  "ROC-AUC:",
  round(
    auc_value,
    3
  ),
  "\n"
)

cat("\nGenerated figures:\n")
cat("09_random_forest_ROC.png\n")
cat("10_random_forest_feature_importance.png\n")
