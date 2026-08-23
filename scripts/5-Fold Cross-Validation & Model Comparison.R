# ============================================================
# Lung Cancer Bioinformatics Project
# STEP 8: 5-Fold Cross-Validation + Model Comparison
# Dataset: GSE10072
# Models: Random Forest + SVM
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACKAGES
# ============================================================

required_packages <- c(
  "limma",
  "randomForest",
  "e1071",
  "pROC"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
    
  }
  
}

library(limma)
library(randomForest)
library(e1071)
library(pROC)

# ============================================================
# 2. LOAD DATA
# ============================================================

expr <- readRDS(
  "data/raw/GSE10072_expression.rds"
)

group <- readRDS(
  "data/processed/sample_groups.rds"
)

group <- factor(
  group,
  levels = c(
    "Normal",
    "Tumor"
  )
)

cat("============================================\n")
cat("STEP 8: MODEL VALIDATION\n")
cat("============================================\n")

cat("\nTotal probes:", nrow(expr), "\n")
cat("Total samples:", ncol(expr), "\n")

cat("\nGroup distribution:\n")
print(table(group))

# ============================================================
# 3. CREATE STRATIFIED 5-FOLD SPLITS
# ============================================================

set.seed(123)

normal_indices <- sample(
  which(group == "Normal")
)

tumor_indices <- sample(
  which(group == "Tumor")
)

normal_folds <- split(
  normal_indices,
  rep(
    1:5,
    length.out = length(normal_indices)
  )
)

tumor_folds <- split(
  tumor_indices,
  rep(
    1:5,
    length.out = length(tumor_indices)
  )
)

folds <- vector(
  "list",
  5
)

for (i in 1:5) {
  
  folds[[i]] <- sort(
    c(
      normal_folds[[i]],
      tumor_folds[[i]]
    )
  )
  
}

cat("\n✓ 5 stratified folds created.\n")

# ============================================================
# 4. STORAGE
# ============================================================

rf_auc <- numeric(5)
svm_auc <- numeric(5)

rf_accuracy <- numeric(5)
svm_accuracy <- numeric(5)

# ============================================================
# 5. CROSS-VALIDATION
# ============================================================

for (fold in 1:5) {
  
  cat("\n")
  cat("============================================\n")
  cat("FOLD", fold, "OF 5\n")
  cat("============================================\n")
  
  # ----------------------------------------------------------
  # Test and training indices
  # ----------------------------------------------------------
  
  test_indices <- folds[[fold]]
  
  train_indices <- setdiff(
    seq_len(ncol(expr)),
    test_indices
  )
  
  train_expr <- expr[
    ,
    train_indices,
    drop = FALSE
  ]
  
  test_expr <- expr[
    ,
    test_indices,
    drop = FALSE
  ]
  
  train_group <- group[
    train_indices
  ]
  
  test_group <- group[
    test_indices
  ]
  
  cat(
    "Training samples:",
    length(train_group),
    "\n"
  )
  
  cat(
    "Testing samples:",
    length(test_group),
    "\n"
  )
  
  # ==========================================================
  # 6. FEATURE SELECTION
  # TRAINING DATA ONLY
  # ==========================================================
  
  cat("\nSelecting features...\n")
  
  design <- model.matrix(
    ~ train_group
  )
  
  colnames(design) <- c(
    "Intercept",
    "Tumor_vs_Normal"
  )
  
  fit <- lmFit(
    train_expr,
    design
  )
  
  fit <- eBayes(
    fit
  )
  
  deg <- topTable(
    fit,
    coef = "Tumor_vs_Normal",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  
  number_features <- min(
    50,
    nrow(deg)
  )
  
  selected_features <- rownames(
    deg
  )[1:number_features]
  
  selected_features <- selected_features[
    selected_features %in% rownames(expr)
  ]
  
  cat(
    "Selected features:",
    length(selected_features),
    "\n"
  )
  
  # ==========================================================
  # 7. CREATE TRAINING / TESTING MATRICES
  # ==========================================================
  
  train_x <- t(
    train_expr[
      selected_features,
      ,
      drop = FALSE
    ]
  )
  
  test_x <- t(
    test_expr[
      selected_features,
      ,
      drop = FALSE
    ]
  )
  
  # Safe feature names
  
  safe_names <- make.names(
    colnames(train_x),
    unique = TRUE
  )
  
  colnames(train_x) <- safe_names
  colnames(test_x) <- safe_names
  
  # ----------------------------------------------------------
  # Remove zero variance features
  # ----------------------------------------------------------
  
  feature_sd <- apply(
    train_x,
    2,
    sd
  )
  
  keep <- is.finite(feature_sd) &
    feature_sd > 0
  
  train_x <- train_x[
    ,
    keep,
    drop = FALSE
  ]
  
  test_x <- test_x[
    ,
    keep,
    drop = FALSE
  ]
  
  # ==========================================================
  # 8. RANDOM FOREST
  # ==========================================================
  
  cat("\nTraining Random Forest...\n")
  
  set.seed(
    100 + fold
  )
  
  rf_model <- randomForest(
    x = train_x,
    y = train_group,
    ntree = 300
  )
  
  rf_pred <- predict(
    rf_model,
    newdata = test_x
  )
  
  rf_prob <- predict(
    rf_model,
    newdata = test_x,
    type = "prob"
  )[, "Tumor"]
  
  # Accuracy
  
  rf_accuracy[fold] <- mean(
    rf_pred == test_group
  )
  
  # AUC
  
  rf_roc <- roc(
    response = test_group,
    predictor = rf_prob,
    levels = c(
      "Normal",
      "Tumor"
    ),
    direction = "<",
    quiet = TRUE
  )
  
  rf_auc[fold] <- as.numeric(
    auc(
      rf_roc
    )
  )
  
  cat(
    "RF Accuracy:",
    round(
      rf_accuracy[fold],
      3
    ),
    "\n"
  )
  
  cat(
    "RF AUC:",
    round(
      rf_auc[fold],
      3
    ),
    "\n"
  )
  
  # ==========================================================
  # 9. SVM
  # ==========================================================
  
  cat("\nTraining SVM...\n")
  
  # ----------------------------------------------------------
  # Scaling based ONLY on training data
  # ----------------------------------------------------------
  
  train_center <- apply(
    train_x,
    2,
    mean
  )
  
  train_scale <- apply(
    train_x,
    2,
    sd
  )
  
  train_scale[
    train_scale == 0
  ] <- 1
  
  train_x_scaled <- scale(
    train_x,
    center = train_center,
    scale = train_scale
  )
  
  test_x_scaled <- scale(
    test_x,
    center = train_center,
    scale = train_scale
  )
  
  # ----------------------------------------------------------
  # SVM model
  # ----------------------------------------------------------
  
  svm_model <- svm(
    x = train_x_scaled,
    y = train_group,
    kernel = "radial",
    probability = TRUE
  )
  
  # ----------------------------------------------------------
  # Predictions
  # ----------------------------------------------------------
  
  svm_pred <- predict(
    svm_model,
    test_x_scaled,
    probability = TRUE
  )
  
  svm_probabilities <- attr(
    svm_pred,
    "probabilities"
  )
  
  svm_prob <- svm_probabilities[
    ,
    "Tumor"
  ]
  
  # Accuracy
  
  svm_accuracy[fold] <- mean(
    svm_pred == test_group
  )
  
  # AUC
  
  svm_roc <- roc(
    response = test_group,
    predictor = svm_prob,
    levels = c(
      "Normal",
      "Tumor"
    ),
    direction = "<",
    quiet = TRUE
  )
  
  svm_auc[fold] <- as.numeric(
    auc(
      svm_roc
    )
  )
  
  cat(
    "SVM Accuracy:",
    round(
      svm_accuracy[fold],
      3
    ),
    "\n"
  )
  
  cat(
    "SVM AUC:",
    round(
      svm_auc[fold],
      3
    ),
    "\n"
  )
}

# ============================================================
# 10. SUMMARY
# ============================================================

cat("\n")
cat("============================================\n")
cat("5-FOLD CROSS-VALIDATION RESULTS\n")
cat("============================================\n")

rf_mean_auc <- mean(
  rf_auc
)

rf_sd_auc <- sd(
  rf_auc
)

rf_mean_accuracy <- mean(
  rf_accuracy
)

svm_mean_auc <- mean(
  svm_auc
)

svm_sd_auc <- sd(
  svm_auc
)

svm_mean_accuracy <- mean(
  svm_accuracy
)

cat(
  "\nRandom Forest Mean AUC:",
  round(
    rf_mean_auc,
    3
  ),
  "\n"
)

cat(
  "Random Forest SD:",
  round(
    rf_sd_auc,
    3
  ),
  "\n"
)

cat(
  "Random Forest Accuracy:",
  round(
    rf_mean_accuracy,
    3
  ),
  "\n"
)

cat(
  "\nSVM Mean AUC:",
  round(
    svm_mean_auc,
    3
  ),
  "\n"
)

cat(
  "SVM SD:",
  round(
    svm_sd_auc,
    3
  ),
  "\n"
)

cat(
  "SVM Accuracy:",
  round(
    svm_mean_accuracy,
    3
  ),
  "\n"
)

# ============================================================
# 11. MODEL COMPARISON TABLE
# ============================================================

model_comparison <- data.frame(
  
  Model = c(
    "Random Forest",
    "SVM"
  ),
  
  Mean_AUC = c(
    rf_mean_auc,
    svm_mean_auc
  ),
  
  SD_AUC = c(
    rf_sd_auc,
    svm_sd_auc
  ),
  
  Mean_Accuracy = c(
    rf_mean_accuracy,
    svm_mean_accuracy
  )
)

print(
  model_comparison
)

# ============================================================
# 12. BEST MODEL
# ============================================================

best_index <- which.max(
  model_comparison$Mean_AUC
)

best_model <- model_comparison$Model[
  best_index
]

best_auc <- model_comparison$Mean_AUC[
  best_index
]

cat(
  "\nBest model:",
  best_model,
  "\n"
)

cat(
  "Best Mean AUC:",
  round(
    best_auc,
    3
  ),
  "\n"
)

# ============================================================
# 13. SAVE RESULTS
# ============================================================

write.csv(
  model_comparison,
  "results/tables/cross_validated_model_comparison.csv",
  row.names = FALSE
)

fold_results <- data.frame(
  
  Fold = 1:5,
  
  RF_AUC = rf_auc,
  
  RF_Accuracy = rf_accuracy,
  
  SVM_AUC = svm_auc,
  
  SVM_Accuracy = svm_accuracy
)

write.csv(
  fold_results,
  "results/tables/fold_level_model_results.csv",
  row.names = FALSE
)

# ============================================================
# 14. MODEL COMPARISON GRAPH
# ============================================================

png(
  "results/figures/11_model_comparison.png",
  width = 2400,
  height = 1800,
  res = 300
)

barplot(
  model_comparison$Mean_AUC,
  
  names.arg = model_comparison$Model,
  
  ylim = c(
    0,
    1
  ),
  
  main = "5-Fold Cross-Validated Model Comparison",
  
  ylab = "Mean ROC-AUC",
  
  xlab = "Model"
)

abline(
  h = 0.5,
  lty = 2
)

dev.off()

# ============================================================
# 15. FINAL
# ============================================================

cat("\n")
cat("============================================\n")
cat("STEP 8 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat(
  "\nRandom Forest Mean AUC:",
  round(
    rf_mean_auc,
    3
  ),
  "\n"
)

cat(
  "SVM Mean AUC:",
  round(
    svm_mean_auc,
    3
  ),
  "\n"
)

cat(
  "Best model:",
  best_model,
  "\n"
)

cat(
  "Best Mean AUC:",
  round(
    best_auc,
    3
  ),
  "\n"
)
