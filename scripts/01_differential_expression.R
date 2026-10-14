#!/usr/bin/env Rscript
# 01_differential_expression.R (v2)
# Phase 2: Extract DE results from pre-computed files
# Uses OSDR pre-computed DE files + EDGE test results + DESeq2 results

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(data.table)
})

results_dir <- "/workspace/osdr_pipeline/results/tables"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: save DE results in standard format
save_de <- function(gene, log2FC, pvalue, padj, organism, accession, contrast, filename) {
  d <- data.frame(
    gene = gene, log2FC = log2FC, pvalue = pvalue, padj = padj,
    direction = ifelse(log2FC > 0, "up", "down"),
    organism = organism, accession = accession, contrast = contrast,
    stringsAsFactors = FALSE
  )
  d <- d[!is.na(d$pvalue), ]
  d <- d[order(d$padj), ]
  write.csv(d, file.path(results_dir, filename), row.names = FALSE)
  cat(sprintf("  %s: %d genes, %d significant (padj<0.05)\n",
              filename, nrow(d), sum(d$padj < 0.05, na.rm = TRUE)))
}

# Helper: EDGE test "Fold change" columns are signed, not ratios: >= 1 for
# higher expression, <= -1 for lower, never between -1 and 1. log2() of the
# raw column turns every negative (down-regulated) value into NaN, which is
# what the OSD-530 and GSE213808 tables used to contain. Convert the
# magnitude and keep the sign, and stop if a value breaks the convention.
# This does not check which group is the numerator; confirm that from the
# group means in the same file.
edge_fc_to_log2 <- function(fc) {
  fc <- as.numeric(fc)
  bad <- !is.na(fc) & abs(fc) < 1
  if (any(bad)) {
    stop(sprintf("%d EDGE fold changes lie strictly between -1 and 1; not the signed convention",
                 sum(bad)))
  }
  sign(fc) * log2(abs(fc))
}

# ============================================================
# OSD-401: Mouse RNA-seq, SpaceX-11/RR-5
# Pre-computed DE file with Space Flight vs Ground Control contrasts
# ============================================================
cat("=== OSD-401 ===\n")
d401 <- fread("/workspace/osdr_pipeline/data/osdr_raw/OSD-401/GLDS-401_rna_seq_differential_expression.csv")
# Use Space Flight vs Ground Control (pooled tissues)
# Column: Log2fc_.Space.Flight...Gastrocnemius.v.Ground.Control...Gastrocnemius.
# Also: Log2fc_.Space.Flight...Quadriceps.femoris.v.Ground.Control...Quadriceps.femoris.
# Use the gastrocnemius comparison as primary (most relevant for atrophy)
fc_col <- grep("Log2fc.*Space.Flight.*Gastrocnemius.*Ground.Control.*Gastrocnemius", colnames(d401), value = TRUE)
pval_col <- grep("^P.value.*Space.Flight.*Gastrocnemius.*Ground.Control.*Gastrocnemius", colnames(d401), value = TRUE)
padj_col <- grep("Adj.p.value.*Space.Flight.*Gastrocnemius.*Ground.Control.*Gastrocnemius", colnames(d401), value = TRUE)

if (length(fc_col) > 0 && length(pval_col) > 0 && length(padj_col) > 0) {
  save_de(d401$SYMBOL, d401[[fc_col[1]]], d401[[pval_col[1]]], d401[[padj_col[1]]],
          "Mus musculus", "OSD-401", "SpaceFlight_vs_GroundControl",
          "OSD-401_DE_SpaceFlight_vs_GroundControl.csv")
} else {
  cat("  WARNING: Could not find Space Flight vs Ground Control columns\n")
  cat("  Available Log2fc columns:\n")
  print(grep("Log2fc", colnames(d401), value = TRUE))
}

# ============================================================
# OSD-684: Human tissue chip RNA-seq, SpaceX-21
# ============================================================
cat("\n=== OSD-684 ===\n")
d684 <- fread("/workspace/osdr_pipeline/data/osdr_raw/OSD-684/GLDS-615_rna_seq_differential_expression_GLbulkRNAseq.csv")
cat("  Cols (first 20):", paste(head(colnames(d684), 20), collapse=", "), "\n")
# Check for Space Flight vs Ground Control
fc_cols <- grep("Log2fc|log2FC|log2fc", colnames(d684), ignore.case = TRUE, value = TRUE)
cat("  Log2FC columns:", paste(fc_cols, collapse="\n    "), "\n")

if (length(fc_cols) > 0) {
  # Use first Space Flight vs Ground Control contrast
  sf_col <- grep("Space.Flight.*Ground.Control|flight.*control", fc_cols, ignore.case = TRUE, value = TRUE)
  if (length(sf_col) > 0) {
    fc_col <- sf_col[1]
    pval_col <- sub("Log2fc", "P.value", fc_col, ignore.case = TRUE)
    padj_col <- sub("Log2fc", "Adj.p.value", fc_col, ignore.case = TRUE)
    # Check if p-value column exists
    if (!(pval_col %in% colnames(d684))) {
      pval_col <- grep(sub("Log2fc", "P.value", fc_col, ignore.case = TRUE), colnames(d684), value = TRUE)[1]
    }
    if (!(padj_col %in% colnames(d684))) {
      padj_col <- grep(sub("Log2fc", "Adj.p.value", fc_col, ignore.case = TRUE), colnames(d684), value = TRUE)[1]
    }
    # Find gene symbol column
    sym_col <- if ("SYMBOL" %in% colnames(d684)) "SYMBOL" else if ("symbol" %in% colnames(d684)) "symbol" else "ENSEMBL"
    save_de(d684[[sym_col]], d684[[fc_col]], d684[[pval_col]], d684[[padj_col]],
            "Homo sapiens", "OSD-684", "SpaceFlight_vs_GroundControl",
            "OSD-684_DE_SpaceFlight_vs_GroundControl.csv")
  } else {
    # Use first available contrast
    fc_col <- fc_cols[1]
    pval_col <- grep("P.value", colnames(d684), value = TRUE)[1]
    padj_col <- grep("Adj.p.value", colnames(d684), value = TRUE)[1]
    sym_col <- if ("SYMBOL" %in% colnames(d684)) "SYMBOL" else "ENSEMBL"
    save_de(d684[[sym_col]], d684[[fc_col]], d684[[pval_col]], d684[[padj_col]],
            "Homo sapiens", "OSD-684", "SpaceFlight_vs_GroundControl",
            "OSD-684_DE_SpaceFlight_vs_GroundControl.csv")
  }
}

# ============================================================
# OSD-21: Mouse microarray, STS-135
# ============================================================
cat("\n=== OSD-21 ===\n")
d21 <- fread("/workspace/osdr_pipeline/data/osdr_raw/OSD-21/GLDS-21_array_differential_expression.csv")
fc_cols21 <- grep("Log2fc|log2FC", colnames(d21), ignore.case = TRUE, value = TRUE)
cat("  Log2FC columns:", paste(head(fc_cols21, 5), collapse="\n    "), "\n")
# Find Space Flight vs Ground Control
sf_col21 <- grep("Space.Flight.*Ground.Control|FLT.*CTRL|FLT.*GC", fc_cols21, ignore.case = TRUE, value = TRUE)
if (length(sf_col21) > 0) {
  fc_col <- sf_col21[1]
  pval_col <- grep("P.value", colnames(d21), value = TRUE)[grep("Space.Flight.*Ground.Control|FLT.*CTRL", grep("P.value", colnames(d21), value = TRUE), ignore.case = TRUE)][1]
  padj_col <- grep("Adj.p.value", colnames(d21), value = TRUE)[grep("Space.Flight.*Ground.Control|FLT.*CTRL", grep("Adj.p.value", colnames(d21), value = TRUE), ignore.case = TRUE)][1]
  sym_col <- if ("SYMBOL" %in% colnames(d21)) "SYMBOL" else if ("Gene.symbol" %in% colnames(d21)) "Gene.symbol" else colnames(d21)[1]
  save_de(d21[[sym_col]], d21[[fc_col]], d21[[pval_col]], d21[[padj_col]],
          "Mus musculus", "OSD-21", "SpaceFlight_vs_GroundControl",
          "OSD-21_DE_SpaceFlight_vs_GroundControl.csv")
} else {
  # Use first contrast
  fc_col <- fc_cols21[1]
  pval_col <- grep("^P.value", colnames(d21), value = TRUE)[1]
  padj_col <- grep("Adj.p.value", colnames(d21), value = TRUE)[1]
  sym_col <- if ("SYMBOL" %in% colnames(d21)) "SYMBOL" else if ("Gene.symbol" %in% colnames(d21)) "Gene.symbol" else colnames(d21)[1]
  save_de(d21[[sym_col]], d21[[fc_col]], d21[[pval_col]], d21[[padj_col]],
          "Mus musculus", "OSD-21", "SpaceFlight_vs_GroundControl",
          "OSD-21_DE_SpaceFlight_vs_GroundControl.csv")
}

# ============================================================
# OSD-111: Mouse microarray, Bion-M1 (EDL/Soleus)
# ============================================================
cat("\n=== OSD-111 ===\n")
d111 <- fread("/workspace/osdr_pipeline/data/osdr_raw/OSD-111/GLDS-111_array_differential_expression.csv")
fc_cols111 <- grep("Log2fc|log2FC", colnames(d111), ignore.case = TRUE, value = TRUE)
cat("  Log2FC columns:", paste(head(fc_cols111, 5), collapse="\n    "), "\n")
sf_col111 <- grep("Space.Flight.*Ground.Control|FLT.*GC|FLT.*CTRL", fc_cols111, ignore.case = TRUE, value = TRUE)
if (length(sf_col111) > 0) {
  fc_col <- sf_col111[1]
  pval_col <- grep("P.value", colnames(d111), value = TRUE)[grep("Space.Flight.*Ground.Control|FLT.*GC|FLT.*CTRL", grep("P.value", colnames(d111), value = TRUE), ignore.case = TRUE)][1]
  padj_col <- grep("Adj.p.value", colnames(d111), value = TRUE)[grep("Space.Flight.*Ground.Control|FLT.*GC|FLT.*CTRL", grep("Adj.p.value", colnames(d111), value = TRUE), ignore.case = TRUE)][1]
  sym_col <- if ("SYMBOL" %in% colnames(d111)) "SYMBOL" else if ("Gene.symbol" %in% colnames(d111)) "Gene.symbol" else colnames(d111)[1]
  save_de(d111[[sym_col]], d111[[fc_col]], d111[[pval_col]], d111[[padj_col]],
          "Mus musculus", "OSD-111", "SpaceFlight_vs_GroundControl",
          "OSD-111_DE_SpaceFlight_vs_GroundControl.csv")
} else if (length(fc_cols111) > 0) {
  fc_col <- fc_cols111[1]
  pval_col <- grep("^P.value", colnames(d111), value = TRUE)[1]
  padj_col <- grep("Adj.p.value", colnames(d111), value = TRUE)[1]
  sym_col <- if ("SYMBOL" %in% colnames(d111)) "SYMBOL" else if ("Gene.symbol" %in% colnames(d111)) "Gene.symbol" else colnames(d111)[1]
  save_de(d111[[sym_col]], d111[[fc_col]], d111[[pval_col]], d111[[padj_col]],
          "Mus musculus", "OSD-111", "SpaceFlight_vs_GroundControl",
          "OSD-111_DE_SpaceFlight_vs_GroundControl.csv")
}

# ============================================================
# OSD-135: Mouse microarray, Bion-M1 (Longissimus dorsi/Tongue)
# ============================================================
cat("\n=== OSD-135 ===\n")
d135 <- fread("/workspace/osdr_pipeline/data/osdr_raw/OSD-135/GLDS-135_array_differential_expression.csv")
fc_cols135 <- grep("Log2fc|log2FC", colnames(d135), ignore.case = TRUE, value = TRUE)
cat("  Log2FC columns:", paste(head(fc_cols135, 5), collapse="\n    "), "\n")
sf_col135 <- grep("Space.Flight.*Ground.Control|FLT.*GC|FLT.*CTRL", fc_cols135, ignore.case = TRUE, value = TRUE)
if (length(sf_col135) > 0) {
  fc_col <- sf_col135[1]
  pval_col <- grep("P.value", colnames(d135), value = TRUE)[grep("Space.Flight.*Ground.Control|FLT.*GC|FLT.*CTRL", grep("P.value", colnames(d135), value = TRUE), ignore.case = TRUE)][1]
  padj_col <- grep("Adj.p.value", colnames(d135), value = TRUE)[grep("Space.Flight.*Ground.Control|FLT.*GC|FLT.*CTRL", grep("Adj.p.value", colnames(d135), value = TRUE), ignore.case = TRUE)][1]
  sym_col <- if ("SYMBOL" %in% colnames(d135)) "SYMBOL" else if ("Gene.symbol" %in% colnames(d135)) "Gene.symbol" else colnames(d135)[1]
  save_de(d135[[sym_col]], d135[[fc_col]], d135[[pval_col]], d135[[padj_col]],
          "Mus musculus", "OSD-135", "SpaceFlight_vs_GroundControl",
          "OSD-135_DE_SpaceFlight_vs_GroundControl.csv")
} else if (length(fc_cols135) > 0) {
  fc_col <- fc_cols135[1]
  pval_col <- grep("^P.value", colnames(d135), value = TRUE)[1]
  padj_col <- grep("Adj.p.value", colnames(d135), value = TRUE)[1]
  sym_col <- if ("SYMBOL" %in% colnames(d135)) "SYMBOL" else if ("Gene.symbol" %in% colnames(d135)) "Gene.symbol" else colnames(d135)[1]
  save_de(d135[[sym_col]], d135[[fc_col]], d135[[pval_col]], d135[[padj_col]],
          "Mus musculus", "OSD-135", "SpaceFlight_vs_GroundControl",
          "OSD-135_DE_SpaceFlight_vs_GroundControl.csv")
}

# ============================================================
# OSD-370: Human bedrest microarray
# ============================================================
cat("\n=== OSD-370 ===\n")
d370 <- fread("/workspace/osdr_pipeline/data/osdr_raw/OSD-370/GLDS-370_array_differential_expression_GLmicroarray.csv")
fc_cols370 <- grep("Log2fc|log2FC", colnames(d370), ignore.case = TRUE, value = TRUE)
cat("  Log2FC columns:", paste(head(fc_cols370, 3), collapse="\n    "), "\n")
# Use During Bed Rest vs Before Bed Rest
br_col <- grep("During.Bed.Rest.*Before.Bed.Rest|Bed.Rest.*baseline.*before", fc_cols370, ignore.case = TRUE, value = TRUE)
if (length(br_col) > 0) {
  fc_col <- br_col[1]
  pval_col <- grep("P.value", colnames(d370), value = TRUE)[grep("During.Bed.Rest.*Before.Bed.Rest|Bed.Rest.*baseline.*before", grep("P.value", colnames(d370), value = TRUE), ignore.case = TRUE)][1]
  padj_col <- grep("Adj.p.value", colnames(d370), value = TRUE)[grep("During.Bed.Rest.*Before.Bed.Rest|Bed.Rest.*baseline.*before", grep("Adj.p.value", colnames(d370), value = TRUE), ignore.case = TRUE)][1]
  sym_col <- if ("SYMBOL" %in% colnames(d370)) "SYMBOL" else if ("Gene.symbol" %in% colnames(d370)) "Gene.symbol" else colnames(d370)[1]
  save_de(d370[[sym_col]], d370[[fc_col]], d370[[pval_col]], d370[[padj_col]],
          "Homo sapiens", "OSD-370", "DuringBedRest_vs_BeforeBedRest",
          "OSD-370_DE_DuringBedRest_vs_BeforeBedRest.csv")
} else if (length(fc_cols370) > 0) {
  fc_col <- fc_cols370[1]
  pval_col <- grep("^P.value", colnames(d370), value = TRUE)[1]
  padj_col <- grep("Adj.p.value", colnames(d370), value = TRUE)[1]
  sym_col <- if ("SYMBOL" %in% colnames(d370)) "SYMBOL" else if ("Gene.symbol" %in% colnames(d370)) "Gene.symbol" else colnames(d370)[1]
  save_de(d370[[sym_col]], d370[[fc_col]], d370[[pval_col]], d370[[padj_col]],
          "Homo sapiens", "OSD-370", "DuringBedRest_vs_BeforeBedRest",
          "OSD-370_DE_DuringBedRest_vs_BeforeBedRest.csv")
}

# ============================================================
# OSD-571: Inspiration4 proteomics
# Pre-computed limma results in xlsx
# ============================================================
cat("\n=== OSD-571 ===\n")
f571 <- "/workspace/osdr_pipeline/data/osdr_raw/OSD-571/GLDS-563_proteomics_Plasma_Proteomics_Processed_Data.xlsx"
sheets571 <- excel_sheets(f571)
cat("  Sheets:", paste(head(sheets571, 5), collapse=", "), "\n")
# Read first sheet - skip metadata rows
df571 <- read_excel(f571, sheet = 1, skip = 6, n_max = 5)
cat("  Cols:", paste(colnames(df571), collapse=", "), "\n")
df571_full <- read_excel(f571, sheet = 1, skip = 6)
if ("Gene" %in% colnames(df571_full) && "logFC" %in% colnames(df571_full)) {
  save_de(df571_full$Gene, df571_full$logFC, df571_full$P.Value, df571_full$adj.P.Val,
          "Homo sapiens", "OSD-571", "PostFlight_vs_PreFlight",
          "OSD-571_DE_PostFlight_vs_PreFlight.csv")
}

# ============================================================
# OSD-530: cfRNA, JAXA CFE
# 9group pairwise file with EDGE test results
# ============================================================
cat("\n=== OSD-530 ===\n")
f530 <- "/workspace/osdr_pipeline/data/osdr_raw/OSD-530/GLDS-530_rna-seq_TGB_050_1_2_64samples_9group_totalcount_all0removed_scalingnormalized_pairwise_analysis_included.xlsx"
df530 <- read_excel(f530, sheet = 1)
cols530 <- colnames(df530)
# Aggregate across Pre vs Flight1-4 using Fisher's method
flight_comparisons <- c("Pre vs Flight1", "Pre vs Flight2", "Pre vs Flight3", "Pre vs Flight4")
de_list_530 <- list()
for (comp in flight_comparisons) {
  pval_col <- cols530[grepl(comp, cols530) & grepl(" - P-value$", cols530) & !grepl("FDR", cols530)]
  fc_col <- cols530[grepl(comp, cols530) & grepl("Fold change$", cols530)]
  fdr_col <- cols530[grepl(comp, cols530) & grepl("FDR p-value correction$", cols530)]
  if (length(pval_col) > 0 && length(fc_col) > 0) {
    de_list_530[[comp]] <- data.frame(
      gene = df530[[1]], log2FC = edge_fc_to_log2(df530[[fc_col[1]]]),
      pvalue = as.numeric(df530[[pval_col[1]]]), padj = as.numeric(df530[[fdr_col[1]]]),
      stringsAsFactors = FALSE
    )
    de_list_530[[comp]] <- de_list_530[[comp]][!is.na(de_list_530[[comp]]$pvalue), ]
  }
}
# Fisher's method
if (length(de_list_530) > 0) {
  all_genes <- unique(unlist(lapply(de_list_530, function(x) x$gene)))
  combined <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  for (comp in names(de_list_530)) {
    d <- de_list_530[[comp]]
    colnames(d)[2:4] <- c(paste0("lFC_", comp), paste0("p_", comp), paste0("padj_", comp))
    combined <- merge(combined, d, by = "gene", all.x = TRUE)
  }
  p_cols <- grep("^p_", colnames(combined), value = TRUE)
  fc_cols <- grep("^lFC_", colnames(combined), value = TRUE)
  combined$mean_log2FC <- rowMeans(combined[, fc_cols], na.rm = TRUE)
  p_mat <- as.matrix(combined[, p_cols])
  p_mat[p_mat == 0] <- 1e-300
  p_mat[is.na(p_mat)] <- 1
  fisher_stat <- -2 * rowSums(log(p_mat), na.rm = TRUE)
  n_studies <- rowSums(!is.na(as.matrix(combined[, p_cols])))
  combined$pvalue <- pchisq(fisher_stat, df = 2 * n_studies, lower.tail = FALSE)
  combined$padj <- p.adjust(combined$pvalue, method = "BH")
  save_de(combined$gene, combined$mean_log2FC, combined$pvalue, combined$padj,
          "Homo sapiens", "OSD-530", "Flight_vs_Pre",
          "OSD-530_DE_Flight_vs_Pre.csv")
}

# ============================================================
# GSE213808: JAXA mouse cfRNA
# EDGE test: MG vs GC (space flight vs ground control)
# ============================================================
cat("\n=== GSE213808 ===\n")
f_geo <- "/workspace/osdr_pipeline/data/geo_raw/GSE213808_TGB_022_3group_N6_totalcount_quantilenormalized.xlsx"
df_geo <- read_excel(f_geo, sheet = 1)
cols_geo <- colnames(df_geo)
# MG vs GC
comp <- "MG vs GC"
pval_col <- cols_geo[grepl(comp, cols_geo) & grepl(" - P-value$", cols_geo) & !grepl("FDR", cols_geo)]
fc_col <- cols_geo[grepl(comp, cols_geo) & grepl("Fold change$", cols_geo)]
fdr_col <- cols_geo[grepl(comp, cols_geo) & grepl("FDR p-value correction$", cols_geo)]
if (length(pval_col) > 0) {
  save_de(df_geo[[1]], edge_fc_to_log2(df_geo[[fc_col[1]]]),
          as.numeric(df_geo[[pval_col[1]]]), as.numeric(df_geo[[fdr_col[1]]]),
          "Mus musculus", "GSE213808", "SpaceFlight_vs_GroundControl",
          "GSE213808_DE_SpaceFlight_vs_GroundControl.csv")
}

# ============================================================
# OSD-569: Inspiration4 direct RNA-seq
# Each sheet has DESeq2 results in last columns
# ============================================================
cat("\n=== OSD-569 ===\n")
f569 <- "/workspace/osdr_pipeline/data/osdr_raw/OSD-569/GLDS-561_long-readRNAseq_Direct_RNA_seq_Gene_Expression_Processed.xlsx"
sheets569 <- excel_sheets(f569)
de_list_569 <- list()
for (sheet_name in sheets569) {
  raw <- read_excel(f569, sheet = sheet_name, n_max = 10, col_names = FALSE)
  row9 <- as.character(unlist(raw[9, ]))
  log2fc_pos <- which(row9 == "log2FC")[1]
  pval_pos <- which(row9 == "p-value")[1]
  padj_pos <- which(row9 == "adjusted p-value")[1]
  if (!is.na(log2fc_pos)) {
    df <- read_excel(f569, sheet = sheet_name, skip = 9, col_names = TRUE)
    de569 <- data.frame(
      gene = sub("\\.\\d+$", "", df[[1]]),
      log2FC = as.numeric(df[[log2fc_pos]]),
      pvalue = as.numeric(df[[pval_pos]]),
      padj = as.numeric(df[[padj_pos]]),
      stringsAsFactors = FALSE
    )
    de569 <- de569[!is.na(de569$pvalue), ]
    de_list_569[[sheet_name]] <- de569
  }
}
# Fisher's method across sheets
if (length(de_list_569) > 0) {
  all_genes <- unique(unlist(lapply(de_list_569, function(x) x$gene)))
  combined569 <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
  for (sn in names(de_list_569)) {
    d <- de_list_569[[sn]]
    colnames(d)[2:4] <- c(paste0("lFC_", sn), paste0("p_", sn), paste0("padj_", sn))
    combined569 <- merge(combined569, d, by = "gene", all.x = TRUE)
  }
  p_cols <- grep("^p_", colnames(combined569), value = TRUE)
  fc_cols <- grep("^lFC_", colnames(combined569), value = TRUE)
  combined569$mean_log2FC <- rowMeans(combined569[, fc_cols], na.rm = TRUE)
  p_mat <- as.matrix(combined569[, p_cols])
  p_mat[p_mat == 0] <- 1e-300
  p_mat[is.na(p_mat)] <- 1
  fisher_stat <- -2 * rowSums(log(p_mat), na.rm = TRUE)
  n_studies <- rowSums(!is.na(as.matrix(combined569[, p_cols])))
  combined569$pvalue <- pchisq(fisher_stat, df = 2 * n_studies, lower.tail = FALSE)
  combined569$padj <- p.adjust(combined569$pvalue, method = "BH")
  save_de(combined569$gene, combined569$mean_log2FC, combined569$pvalue, combined569$padj,
          "Homo sapiens", "OSD-569", "PostFlight_vs_PreFlight",
          "OSD-569_DE_PostFlight_vs_PreFlight.csv")
}

# ============================================================
# Summary
# ============================================================
cat("\n=== DE Results Summary ===\n")
de_files <- list.files(results_dir, pattern = "_DE_.*\\.csv$", full.names = TRUE)
for (f in de_files) {
  d <- read.csv(f)
  cat(sprintf("  %-55s %6d genes, %5d sig (padj<0.05)\n",
              basename(f), nrow(d), sum(d$padj < 0.05, na.rm = TRUE)))
}
cat("\nDone.\n")
