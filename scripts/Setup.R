install.packages(c(
  "tidyverse",
  "pheatmap",
  "RColorBrewer",
  "caret",
  "randomForest",
  "e1071",
  "pROC"
))

BiocManager::install(c(
  "GEOquery",
  "limma",
  "Biobase",
  "SummarizedExperiment",
  "AnnotationDbi",
  "hgu133a.db"
), ask = FALSE, update = FALSE)

library(GEOquery)
library(limma)
library(Biobase)
library(SummarizedExperiment)
library(AnnotationDbi)
library(hgu133a.db)

library(tidyverse)
library(pheatmap)
library(caret)
library(randomForest)
library(pROC)

cat("ENVIRONMENT READY")

folders <- c(
  "data/raw",
  "data/processed",
  "results/figures",
  "results/tables",
  "scripts",
  "report"
)

for (folder in folders) {
  dir.create(
    folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ============================================================
# Lung Cancer Bioinformatics Project
# Step 1: Project Setup
# ============================================================

# Clear environment
rm(list = ls())

# Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------

library(tidyverse)
library(GEOquery)
library(limma)
library(Biobase)
library(SummarizedExperiment)
library(pheatmap)
library(RColorBrewer)
library(caret)
library(randomForest)
library(e1071)
library(pROC)
library(AnnotationDbi)
library(hgu133a.db)

# ------------------------------------------------------------
# Project information
# ------------------------------------------------------------

project_name <- "Lung Cancer Gene Expression Analysis"

dataset_id <- "GSE10072"

cat("============================================\n")
cat("Project:", project_name, "\n")
cat("Dataset:", dataset_id, "\n")
cat("Working directory:", getwd(), "\n")
cat("============================================\n")

