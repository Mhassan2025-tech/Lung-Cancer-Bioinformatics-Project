# ============================================================
# Lung Cancer Bioinformatics Project
# Step 2: GEO Dataset Download & Inspection
# Dataset: GSE10072
# ============================================================

# ------------------------------------------------------------
# 1. Clean environment
# ------------------------------------------------------------

rm(list = ls())

# ------------------------------------------------------------
# 2. Load required packages
# ------------------------------------------------------------

library(GEOquery)
library(Biobase)

# ------------------------------------------------------------
# 3. Dataset information
# ------------------------------------------------------------

gse_id <- "GSE10072"

cat("============================================\n")
cat("GEO DATA DOWNLOAD\n")
cat("Dataset:", gse_id, "\n")
cat("============================================\n\n")

# ------------------------------------------------------------
# 4. Download GEO dataset
# ------------------------------------------------------------

cat("Downloading GEO dataset...\n")
cat("This may take a few minutes.\n\n")

gse <- getGEO(
  gse_id,
  GSEMatrix = TRUE,
  AnnotGPL = TRUE
)

cat("Download completed successfully.\n\n")

# ------------------------------------------------------------
# 5. Check GEO objects
# ------------------------------------------------------------

cat("Number of GEO objects returned:\n")
print(length(gse))

cat("\nGEO object names:\n")
print(names(gse))

# ------------------------------------------------------------
# 6. Select the first platform
# ------------------------------------------------------------

data_obj <- gse[[1]]

cat("\nObject class:\n")
print(class(data_obj))

# ------------------------------------------------------------
# 7. Extract expression matrix
# ------------------------------------------------------------

# GEOquery can return different object types depending
# on the dataset/Bioconductor version.

if (inherits(data_obj, "ExpressionSet")) {
  
  cat("\nDetected ExpressionSet object.\n")
  
  expr <- exprs(data_obj)
  pheno <- pData(data_obj)
  
} else if (inherits(data_obj, "SummarizedExperiment")) {
  
  cat("\nDetected SummarizedExperiment object.\n")
  
  expr <- assay(data_obj)
  pheno <- as.data.frame(colData(data_obj))
  
} else {
  
  stop(
    paste(
      "Unsupported GEO object type:",
      paste(class(data_obj), collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 8. Convert phenotype data to data frame
# ------------------------------------------------------------

pheno <- as.data.frame(pheno)

# ------------------------------------------------------------
# 9. Dataset dimensions
# ------------------------------------------------------------

cat("\n============================================\n")
cat("DATASET INFORMATION\n")
cat("============================================\n")

cat("\nExpression matrix dimensions:\n")
cat("Genes/Probes:", nrow(expr), "\n")
cat("Samples:", ncol(expr), "\n")

cat("\nPhenotype table dimensions:\n")
cat("Samples:", nrow(pheno), "\n")
cat("Variables:", ncol(pheno), "\n")

# ------------------------------------------------------------
# 10. Verify sample numbers
# ------------------------------------------------------------

if (ncol(expr) != nrow(pheno)) {
  
  stop(
    "ERROR: Number of expression samples does not match phenotype samples."
  )
  
} else {
  
  cat("\n✓ Expression and phenotype sample counts match.\n")
}

# ------------------------------------------------------------
# 11. Check sample names
# ------------------------------------------------------------

cat("\nFirst 10 expression sample names:\n")
print(colnames(expr)[1:min(10, ncol(expr))])

# ------------------------------------------------------------
# 12. Check phenotype columns
# ------------------------------------------------------------

cat("\n============================================\n")
cat("PHENOTYPE INFORMATION\n")
cat("============================================\n")

cat("\nAvailable phenotype columns:\n")

print(colnames(pheno))

# ------------------------------------------------------------
# 13. Display important phenotype information
# ------------------------------------------------------------

cat("\nFirst 6 rows of phenotype data:\n\n")

print(
  pheno[
    1:min(6, nrow(pheno)),
    1:min(10, ncol(pheno)),
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# 14. Check sample titles
# ------------------------------------------------------------

if ("title" %in% colnames(pheno)) {
  
  cat("\nSample titles:\n")
  print(pheno$title)
  
}

# ------------------------------------------------------------
# 15. Check source information
# ------------------------------------------------------------

if ("source_name_ch1" %in% colnames(pheno)) {
  
  cat("\nSource information:\n")
  print(pheno$source_name_ch1)
  
}

# ------------------------------------------------------------
# 16. Check characteristics columns
# ------------------------------------------------------------

characteristic_columns <- grep(
  "characteristics",
  colnames(pheno),
  ignore.case = TRUE,
  value = TRUE
)

if (length(characteristic_columns) > 0) {
  
  cat("\nCharacteristics columns found:\n")
  print(characteristic_columns)
  
  cat("\nFirst 6 rows of characteristics:\n")
  
  print(
    pheno[
      1:min(6, nrow(pheno)),
      characteristic_columns,
      drop = FALSE
    ]
  )
  
}

# ------------------------------------------------------------
# 17. Expression data quality check
# ------------------------------------------------------------

cat("\n============================================\n")
cat("EXPRESSION DATA QC\n")
cat("============================================\n")

cat("\nMissing expression values:\n")

missing_values <- sum(is.na(expr))

print(missing_values)

if (missing_values == 0) {
  
  cat("✓ No missing expression values detected.\n")
  
} else {
  
  cat("⚠ Missing values detected.\n")
  
}

# ------------------------------------------------------------
# 18. Expression value summary
# ------------------------------------------------------------

cat("\nExpression value summary:\n")

print(
  summary(
    as.vector(expr)
  )
)

# ------------------------------------------------------------
# 19. Check expression matrix type
# ------------------------------------------------------------

cat("\nExpression matrix class:\n")
print(class(expr))

# ------------------------------------------------------------
# 20. Save raw objects
# ------------------------------------------------------------

cat("\n============================================\n")
cat("SAVING DATA\n")
cat("============================================\n")

saveRDS(
  data_obj,
  "data/raw/GSE10072_data_object.rds"
)

saveRDS(
  expr,
  "data/raw/GSE10072_expression.rds"
)

saveRDS(
  pheno,
  "data/raw/GSE10072_phenotype.rds"
)

cat("\n✓ GEO object saved.\n")
cat("✓ Expression matrix saved.\n")
cat("✓ Phenotype data saved.\n")

# ------------------------------------------------------------
# 21. Final status
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 2 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat("\nDataset:", gse_id, "\n")
cat("Genes/Probes:", nrow(expr), "\n")
cat("Samples:", ncol(expr), "\n")
cat("Missing values:", missing_values, "\n")
