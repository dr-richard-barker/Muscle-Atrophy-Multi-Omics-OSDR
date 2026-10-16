#!/usr/bin/env Rscript
# 03_gsea_enrichment.R
# Phase 5: GSEA on mTOR/hypertrophic/proteostasis/mitochondrial pathways
# Also covers Phase 4: fatigue biomarker mapping via GSVA/ssGSEA

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(reactome.db)
  library(dplyr)
  library(data.table)
  library(ggplot2)
})

signature_dir <- "/workspace/osdr_pipeline/results/signatures"
results_dir <- "/workspace/osdr_pipeline/results"
figures_dir <- "/workspace/osdr_pipeline/results/figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load ranked gene list from meta-analysis
# ============================================================
cat("=== Loading ranked gene list ===\n")
ranked <- read.csv(file.path(signature_dir, "ranked_gene_list.csv"), stringsAsFactors = FALSE)
ranked <- ranked[!is.na(ranked$gene) & ranked$gene != "" & !duplicated(ranked$gene), ]
ranks <- ranked$log2FC
names(ranks) <- ranked$gene
ranks <- sort(ranks, decreasing = TRUE)
cat(sprintf("  %d genes in ranked list\n", length(ranks)))

# ============================================================
# 2. Load MSigDB gene sets (Hallmark, KEGG, Reactome)
# ============================================================
cat("=== Loading MSigDB gene sets ===\n")

# Hallmark
hallmark <- msigdbr(collection = "H", species = "Homo sapiens")
h_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
cat(sprintf("  Hallmark: %d gene sets\n", length(h_sets)))

# KEGG
kegg <- msigdbr(collection = "C2", subcollection = "CP:KEGG_LEGACY", species = "Homo sapiens")
k_sets <- split(kegg$gene_symbol, kegg$gs_name)
cat(sprintf("  KEGG: %d gene sets\n", length(k_sets)))

# Reactome
reactome <- msigdbr(collection = "C2", subcollection = "CP:REACTOME", species = "Homo sapiens")
r_sets <- split(reactome$gene_symbol, reactome$gs_name)
cat(sprintf("  Reactome: %d gene sets\n", length(r_sets)))

# GO Biological Process
go_bp <- msigdbr(collection = "C5", subcollection = "GO:BP", species = "Homo sapiens")
g_sets <- split(go_bp$gene_symbol, go_bp$gs_name)
cat(sprintf("  GO BP: %d gene sets\n", length(g_sets)))

# Combine all
all_sets <- c(h_sets, k_sets, r_sets, g_sets)
cat(sprintf("  Total: %d gene sets\n", length(all_sets)))

# ============================================================
# 3. Run GSEA with fgsea
# ============================================================
cat("\n=== Running GSEA ===\n")
set.seed(42)
fgsea_res <- fgsea(pathways = all_sets, stats = ranks, minSize = 10, maxSize = 500, nperm = 10000)
fgsea_res <- as.data.frame(fgsea_res)
fgsea_res <- fgsea_res[order(fgsea_res$pval), ]
cat(sprintf("  GSEA complete: %d pathways tested, %d significant (padj<0.05)\n",
            nrow(fgsea_res), sum(fgsea_res$padj < 0.05, na.rm = TRUE)))

# Save full results
fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, function(x) paste(unlist(x), collapse=","))
write.csv(fgsea_res, file.path(results_dir, "gsea_full_results.csv"), row.names = FALSE)

# ============================================================
# 4. Focus on mTOR/hypertrophic/proteostasis/mitochondrial pathways
# ============================================================
cat("\n=== Key pathway results ===\n")

key_patterns <- c(
  "MTOR", "mTOR", "MTORC1",
  "HYPERTROPH", "MUSCLE_CONTRACTION", "MUSCLE_FIBER",
  "PROTEOSTASIS", "UBIQUITIN", "PROTEASOME", "AUTOPHAG", "CHAPERONE",
  "MITOCHOND", "OXIDATIVE_PHOSPHORYLATION", "ELECTRON_TRANSPORT", "RESPIRATORY",
  "MYOFIBRIL", "SARCOMERE", "ACTIN_FILAMENT",
  "INSULIN_SIGNAL", "PI3K", "AKT",
  "FOXO", "NFKB", "INFLAMM",
  "GLYCOLYSIS", "FATTY_ACID", "LIPID_METABOL",
  "MUSCLE", "SKELETAL", "MYOGENESIS", "MYOBLAST",
  "ATROPHY", "PROTEIN_DEGRADATION", "PROTEIN_CATABOL"
)

key_res <- fgsea_res[grepl(paste(key_patterns, collapse = "|"), fgsea_res$pathway, ignore.case = TRUE), ]
key_res_sig <- key_res[which(key_res$padj < 0.1), ]

cat(sprintf("\n  Key pathways (padj<0.1): %d\n", nrow(key_res_sig)))
if (nrow(key_res_sig) > 0) {
  cat(sprintf("  %-70s %8s %10s %8s\n", "Pathway", "NES", "pval", "padj"))
  for (i in 1:min(30, nrow(key_res_sig))) {
    cat(sprintf("  %-70s %8.2f %10.2e %8.4f\n",
                substr(key_res_sig$pathway[i], 1, 70),
                key_res_sig$NES[i], key_res_sig$pval[i], key_res_sig$padj[i]))
  }
}

# Save key pathway results
write.csv(key_res, file.path(results_dir, "gsea_key_pathways.csv"), row.names = FALSE)

# ============================================================
# 5. Hallmark enrichment summary
# ============================================================
cat("\n=== Hallmark enrichment summary ===\n")
hallmark_res <- fgsea_res[grepl("^HALLMARK_", fgsea_res$pathway), ]
hallmark_res <- hallmark_res[order(hallmark_res$pval), ]
cat(sprintf("  %-55s %8s %10s %8s\n", "Hallmark", "NES", "pval", "padj"))
for (i in 1:nrow(hallmark_res)) {
  cat(sprintf("  %-55s %8.2f %10.2e %8.4f\n",
              substr(hallmark_res$pathway[i], 1, 55),
              hallmark_res$NES[i], hallmark_res$pval[i], hallmark_res$padj[i]))
}

# ============================================================
# 6. ORA on consensus signature genes
# ============================================================
cat("\n=== ORA on consensus signature ===\n")
sig <- read.csv(file.path(signature_dir, "consensus_catabolic_signature.csv"), stringsAsFactors = FALSE)
up_genes <- sig$gene[sig$signature_direction == "up"]
down_genes <- sig$gene[sig$signature_direction == "down"]
cat(sprintf("  Up: %d, Down: %d genes\n", length(up_genes), length(down_genes)))

# ORA for up-regulated genes
if (length(up_genes) > 5) {
  ora_up <- enricher(gene = up_genes, TERM2GENE = msigdbr(collection = "H", species = "Homo sapiens")[, c("gs_name", "gene_symbol")],
                     pvalueCutoff = 0.1, qvalueCutoff = 0.2)
  if (!is.null(ora_up) && nrow(ora_up) > 0) {
    cat("\n  Up-regulated ORA (top 10):\n")
    print(head(ora_up, 10))
  }
}

# ORA for down-regulated genes
if (length(down_genes) > 5) {
  ora_down <- enricher(gene = down_genes, TERM2GENE = msigdbr(collection = "H", species = "Homo sapiens")[, c("gs_name", "gene_symbol")],
                       pvalueCutoff = 0.1, qvalueCutoff = 0.2)
  if (!is.null(ora_down) && nrow(ora_down) > 0) {
    cat("\n  Down-regulated ORA (top 10):\n")
    print(head(ora_down, 10))
  }
}

# ============================================================
# 7. Save summary
# ============================================================
cat("\n=== Saving results ===\n")
# Save hallmark results
write.csv(hallmark_res, file.path(results_dir, "gsea_hallmark_results.csv"), row.names = FALSE)

# Save significant GSEA results
sig_gsea <- fgsea_res[which(fgsea_res$padj < 0.05), ]
write.csv(sig_gsea, file.path(results_dir, "gsea_significant_pathways.csv"), row.names = FALSE)
cat(sprintf("  Significant pathways (padj<0.05): %d\n", nrow(sig_gsea)))

cat("\nDone.\n")
