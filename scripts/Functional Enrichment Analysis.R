# ============================================================
# Lung Cancer Bioinformatics Project
# STEP 9: Functional Enrichment Analysis
# Dataset: GSE10072
# Analysis: GO + KEGG
# ============================================================

rm(list = ls())

# ============================================================
# 1. INSTALL / LOAD REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "ggplot2"
)

# Install BiocManager if necessary

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
  
}

# Install missing Bioconductor packages

for (pkg in required_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
    
  }
  
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# ============================================================
# 2. LOAD DEG RESULTS
# ============================================================

deg <- read.csv(
  "results/tables/annotated_differential_expression_results.csv",
  stringsAsFactors = FALSE
)

# ============================================================
# 3. CHECK DATA
# ============================================================

cat("============================================\n")
cat("STEP 9: FUNCTIONAL ENRICHMENT ANALYSIS\n")
cat("============================================\n")

cat(
  "\nTotal annotated probes:",
  nrow(deg),
  "\n"
)

# ------------------------------------------------------------
# Check Gene Symbol column
# ------------------------------------------------------------

if (!"Gene_Symbol" %in% colnames(deg)) {
  
  stop(
    "ERROR: Gene_Symbol column not found."
  )
  
}

# ============================================================
# 4. SELECT SIGNIFICANT GENES
# ============================================================

significant_deg <- deg[
  
  deg$adj.P.Val < 0.05 &
    abs(deg$logFC) >= 1 &
    !is.na(deg$Gene_Symbol) &
    deg$Gene_Symbol != "",
  
]

# Remove duplicate symbols

gene_symbols <- unique(
  significant_deg$Gene_Symbol
)

cat(
  "\nSignificant unique genes:",
  length(gene_symbols),
  "\n"
)

if (length(gene_symbols) < 10) {
  
  stop(
    "ERROR: Too few significant genes for enrichment analysis."
  )
  
}

# ============================================================
# 5. CONVERT SYMBOLS → ENTREZ IDs
# ============================================================

cat("\nConverting gene symbols to Entrez IDs...\n")

gene_conversion <- bitr(
  
  gene_symbols,
  
  fromType = "SYMBOL",
  
  toType = "ENTREZID",
  
  OrgDb = org.Hs.eg.db
  
)

# Remove duplicates

gene_conversion <- gene_conversion[
  !duplicated(
    gene_conversion$ENTREZID
  ),
]

entrez_genes <- unique(
  gene_conversion$ENTREZID
)

cat(
  "Successfully converted:",
  length(entrez_genes),
  "genes\n"
)

# ============================================================
# 6. CREATE BACKGROUND GENE LIST
# ============================================================

cat("\nCreating background gene set...\n")

background_symbols <- unique(
  deg$Gene_Symbol[
    !is.na(deg$Gene_Symbol) &
      deg$Gene_Symbol != ""
  ]
)

background_conversion <- bitr(
  
  background_symbols,
  
  fromType = "SYMBOL",
  
  toType = "ENTREZID",
  
  OrgDb = org.Hs.eg.db
  
)

background_entrez <- unique(
  background_conversion$ENTREZID
)

cat(
  "Background genes:",
  length(background_entrez),
  "\n"
)

# ============================================================
# 7. GENE ONTOLOGY — BIOLOGICAL PROCESS
# ============================================================

cat("\nRunning GO Biological Process enrichment...\n")

go_bp <- enrichGO(
  
  gene = entrez_genes,
  
  universe = background_entrez,
  
  OrgDb = org.Hs.eg.db,
  
  keyType = "ENTREZID",
  
  ont = "BP",
  
  pAdjustMethod = "BH",
  
  pvalueCutoff = 0.05,
  
  qvalueCutoff = 0.20,
  
  readable = TRUE
  
)

cat(
  "GO terms found:",
  nrow(as.data.frame(go_bp)),
  "\n"
)

# ============================================================
# 8. SAVE GO RESULTS
# ============================================================

go_results <- as.data.frame(
  go_bp
)

write.csv(
  
  go_results,
  
  "results/tables/GO_Biological_Process_enrichment.csv",
  
  row.names = FALSE
  
)

# ============================================================
# 9. GO DOT PLOT
# ============================================================

if (
  nrow(go_results) > 0
) {
  
  go_plot <- dotplot(
    
    go_bp,
    
    showCategory = 15,
    
    font.size = 10
    
  ) +
    
    ggtitle(
      "GO Biological Process Enrichment"
    )
  
  ggsave(
    
    "results/figures/12_GO_BP_dotplot.png",
    
    go_plot,
    
    width = 10,
    
    height = 8,
    
    dpi = 300
    
  )
  
  cat(
    "✓ GO dotplot saved.\n"
  )
  
} else {
  
  cat(
    "No significant GO terms detected.\n"
  )
  
}

# ============================================================
# 10. KEGG PATHWAY ENRICHMENT
# ============================================================

cat("\nRunning KEGG pathway enrichment...\n")

kegg_results <- enrichKEGG(
  
  gene = entrez_genes,
  
  universe = background_entrez,
  
  organism = "hsa",
  
  pvalueCutoff = 0.05,
  
  pAdjustMethod = "BH",
  
  qvalueCutoff = 0.20
  
)

kegg_df <- as.data.frame(
  kegg_results
)

cat(
  "KEGG pathways found:",
  nrow(kegg_df),
  "\n"
)

# ============================================================
# 11. SAVE KEGG RESULTS
# ============================================================

write.csv(
  
  kegg_df,
  
  "results/tables/KEGG_pathway_enrichment.csv",
  
  row.names = FALSE
  
)

# ============================================================
# 12. KEGG DOT PLOT
# ============================================================

if (
  nrow(kegg_df) > 0
) {
  
  kegg_plot <- dotplot(
    
    kegg_results,
    
    showCategory = 15,
    
    font.size = 10
    
  ) +
    
    ggtitle(
      "KEGG Pathway Enrichment"
    )
  
  ggsave(
    
    "results/figures/13_KEGG_dotplot.png",
    
    kegg_plot,
    
    width = 10,
    
    height = 8,
    
    dpi = 300
    
  )
  
  cat(
    "✓ KEGG dotplot saved.\n"
  )
  
} else {
  
  cat(
    "No significant KEGG pathways detected.\n"
  )
  
}

# ============================================================
# 13. TOP GO TERMS
# ============================================================

if (
  nrow(go_results) > 0
) {
  
  top_go <- go_results[
    
    order(
      go_results$p.adjust
    ),
    
  ]
  
  top_go <- head(
    top_go,
    20
  )
  
  write.csv(
    
    top_go,
    
    "results/tables/top_20_GO_terms.csv",
    
    row.names = FALSE
    
  )
  
}

# ============================================================
# 14. TOP KEGG PATHWAYS
# ============================================================

if (
  nrow(kegg_df) > 0
) {
  
  top_kegg <- kegg_df[
    
    order(
      kegg_df$p.adjust
    ),
    
  ]
  
  top_kegg <- head(
    top_kegg,
    20
  )
  
  write.csv(
    
    top_kegg,
    
    "results/tables/top_20_KEGG_pathways.csv",
    
    row.names = FALSE
    
  )
  
}

# ============================================================
# 15. SAVE GENE LISTS
# ============================================================

write.csv(
  
  data.frame(
    Gene_Symbol = gene_symbols
  ),
  
  "results/tables/significant_gene_symbols.csv",
  
  row.names = FALSE
  
)

write.csv(
  
  gene_conversion,
  
  "results/tables/gene_symbol_to_entrez_conversion.csv",
  
  row.names = FALSE
  
)

# ============================================================
# 16. FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================\n")
cat("STEP 9 COMPLETED SUCCESSFULLY ✓\n")
cat("============================================\n")

cat(
  "\nSignificant genes:",
  length(gene_symbols),
  "\n"
)

cat(
  "Entrez IDs:",
  length(entrez_genes),
  "\n"
)

cat(
  "GO Biological Processes:",
  nrow(go_results),
  "\n"
)

cat(
  "KEGG pathways:",
  nrow(kegg_df),
  "\n"
)

cat("\nGenerated files:\n")

cat(
  "12_GO_BP_dotplot.png\n"
)

cat(
  "13_KEGG_dotplot.png\n"
)

cat(
  "GO_Biological_Process_enrichment.csv\n"
)

cat(
  "KEGG_pathway_enrichment.csv\n"
)
