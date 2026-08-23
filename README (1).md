# Lung Cancer Bioinformatics Analysis

## Overview
Computational analysis of lung cancer gene-expression data from **GSE10072**, integrating differential expression, machine learning, functional enrichment, and candidate biomarker identification.

## Objectives
- Perform quality control and PCA.
- Identify differentially expressed genes (DEGs).
- Annotate DEGs and generate heatmaps.
- Build Random Forest and SVM classifiers.
- Evaluate models with stratified 5-fold cross-validation.
- Perform GO Biological Process and KEGG enrichment.
- Integrate DEG significance with machine-learning feature importance.
- Identify computational candidate biomarkers.

## Dataset
**Dataset:** GSE10072  
**Source:** NCBI Gene Expression Omnibus (GEO)  
**Organism:** *Homo sapiens*

## Workflow
```text
GSE10072
   ↓
Expression Data → Quality Control → PCA
   ↓
Differential Expression → Annotation → Heatmap
   ↓
Random Forest + SVM → 5-Fold Cross-Validation
   ↓
GO + KEGG Enrichment
   ↓
DEG + ML Integration
   ↓
Candidate Biomarkers → Biological Interpretation
```

## Methods

### Quality Control
Expression distributions, density, and sample correlations were examined to assess data quality.

### PCA
Principal Component Analysis was used to visualize major sources of variation and separation between tumor and normal samples.

### Differential Expression
Differential expression was performed using **limma**.

Criteria:
- Adjusted P-value < 0.05
- Absolute log2 fold change >= 1

### Gene Annotation and Heatmap
Significant probes were mapped to human gene symbols, and top DEGs were visualized using a standardized expression heatmap.

### Machine Learning
Random Forest and Support Vector Machine (SVM) classifiers were evaluated using stratified 5-fold cross-validation.

### Functional Enrichment
Significant genes were analyzed using:
- Gene Ontology (GO) Biological Process enrichment
- KEGG pathway enrichment

### Integrated Biomarker Analysis
Differential-expression results were integrated with Random Forest feature importance to identify genes supported by both statistical and machine-learning evidence.

## Key Findings

High-ranking integrated candidate genes included:

**VWF, LDB2, CDH5, EDNRB, GRK5, FHL1, PECAM1, DUOX1, FOXF1, GPM6A**

Other candidates included **CAV1, CD36, CLEC3B, TCF21, TEK, FAM107A, AGER, SPP1, ABCA8, and FABP4**.

Enrichment analysis highlighted:
- Extracellular matrix organization
- Extracellular structure organization
- Angiogenesis
- Vasculature development
- Sprouting angiogenesis
- Integrin-related signaling
- ECM-receptor interaction
- Complement and coagulation cascades

Together, the findings suggest a strong signal involving **vascular/endothelial biology, angiogenesis, and extracellular-matrix remodeling**.

## Machine Learning Results

Random Forest achieved a mean 5-fold ROC-AUC of **1.000**.

SVM also showed strong performance, with fold-level ROC-AUC values ranging from approximately **0.905 to 1.000**.

The perfect Random Forest result should **not** be interpreted as clinical validation. Independent external validation is required to assess generalizability.

## Limitations
- Single GEO dataset.
- No independent external validation.
- Very high model performance may be dataset-specific.
- Candidate biomarkers are computational candidates, not clinically validated biomarkers.
- Enrichment identifies associations and does not establish causality.
- Experimental and independent-cohort validation is required.

## Reproducibility

The project was developed in **R/RStudio**.

Major packages include:
`limma`, `randomForest`, `e1071`, `pROC`, `pheatmap`, `clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `hgu133a.db`, `ggplot2`, and `tidyverse`.

Project outputs are organized into:

```text
data/
scripts/
results/
├── figures/
└── tables/
```

## Project Structure

```text
Lung-Cancer-Bioinformatics/
├── data/
│   ├── raw/
│   └── processed/
├── scripts/
├── results/
│   ├── figures/
│   └── tables/
├── Lung-Cancer-Bioinformatics.Rproj
└── README.md
```

## AI Assistance Disclosure

AI tools were used as a coding and troubleshooting assistant for R code development, debugging, workflow planning, bioinformatics concept explanations, and visualization guidance.

The analyses were executed in the local R environment, and the generated outputs were reviewed and interpreted as part of the project workflow.

## Conclusion

This project demonstrates an integrated workflow for analyzing lung cancer gene-expression data using statistical analysis, machine learning, and functional enrichment.

The combined results highlighted **vascular/endothelial biology, angiogenesis, and extracellular-matrix remodeling** as important biological themes.

The identified genes should be considered **computational candidate biomarkers for further investigation**, rather than clinically validated biomarkers.

## Author

**Muhammad Hassan**  
BS Biotechnology  
Bioinformatics | Data Analysis | Machine Learning

## Disclaimer

This project is intended for educational and research purposes. The identified genes, pathways, and machine-learning models should not be interpreted as clinically validated diagnostic or prognostic tools.
