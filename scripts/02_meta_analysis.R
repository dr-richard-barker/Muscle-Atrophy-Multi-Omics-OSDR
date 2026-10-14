#!/usr/bin/env Rscript
# 02_meta_analysis.R
# Phase 3: Cross-species meta-analysis & consensus catabolic signature
# Uses homologene R package for ortholog mapping

suppressPackageStartupMessages({
  library(metafor)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(homologene)
  library(dplyr)
  library(data.table)
})

results_dir <- "/workspace/osdr_pipeline/results/tables"
signature_dir <- "/workspace/osdr_pipeline/results/signatures"
dir.create(signature_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Load all DE results
# ============================================================
cat("=== Loading DE results ===\n")
de_files <- list.files(results_dir, pattern = "_DE_.*\\.csv$", full.names = TRUE)
de_list <- list()
for (f in de_files) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  study_id <- paste(d$accession[1], d$contrast[1], sep = "_")
  de_list[[study_id]] <- d
  cat(sprintf("  %s: %d genes, org=%s, ID example: %s\n",
              study_id, nrow(d), d$organism[1], d$gene[1]))
}

# ============================================================
# 2. Build ortholog mapping using homologene
# ============================================================
cat("\n=== Building ortholog table via homologene ===\n")
data(homologeneData, package = "homologene")
mouse_h <- homologeneData[homologeneData$Taxonomy == 10090, c("HID", "Gene.ID", "Gene.Symbol")]
human_h <- homologeneData[homologeneData$Taxonomy == 9606, c("HID", "Gene.ID", "Gene.Symbol")]
colnames(mouse_h) <- c("HID", "mouse_entrez", "mouse_symbol")
colnames(human_h) <- c("HID", "human_entrez", "human_symbol")
ortho_table <- merge(mouse_h, human_h, by = "HID")
ortho_table <- ortho_table[, c("mouse_symbol", "mouse_entrez", "human_symbol", "human_entrez")]
ortho_table <- ortho_table[!duplicated(ortho_table$mouse_symbol), ]
cat(sprintf("  Ortholog pairs: %d mouse -> human\n", nrow(ortho_table)))

# ============================================================
# 3. Map each study's gene IDs to human HGNC symbols
# ============================================================
cat("\n=== Mapping study gene IDs to human symbols ===\n")

map_study_to_human <- function(de_df, ortho_table) {
  org <- de_df$organism[1]
  genes <- as.character(de_df$gene)
  acc <- de_df$accession[1]
  ctr <- de_df$contrast[1]

  if (org == "Homo sapiens") {
    # Determine ID type by checking majority of genes
    n_ens <- sum(grepl("^ENSG\\d+", genes))
    n_entrez <- sum(grepl("^\\d+$", genes))
    if (n_ens > length(genes) * 0.5) {
      genes_clean <- sub("\\.\\d+$", "", genes)
      # Filter to only valid Ensembl IDs
      valid_ens <- unique(genes_clean[genes_clean %in% keys(org.Hs.eg.db, keytype = "ENSEMBL")])
      if (length(valid_ens) > 0) {
        mapping <- AnnotationDbi::select(org.Hs.eg.db, keys = valid_ens,
                                          columns = "SYMBOL", keytype = "ENSEMBL")
        mapping <- mapping[!is.na(mapping$SYMBOL) & mapping$SYMBOL != "", ]
        mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
        de_df$hgnc_symbol <- mapping$SYMBOL[match(genes_clean, mapping$ENSEMBL)]
      } else {
        de_df$hgnc_symbol <- NA
      }
    } else if (n_entrez > length(genes) * 0.5) {
      valid_entrez <- unique(genes[genes %in% keys(org.Hs.eg.db, keytype = "ENTREZID")])
      if (length(valid_entrez) > 0) {
        mapping <- AnnotationDbi::select(org.Hs.eg.db, keys = valid_entrez,
                                          columns = "SYMBOL", keytype = "ENTREZID")
        mapping <- mapping[!is.na(mapping$SYMBOL) & mapping$SYMBOL != "", ]
        mapping <- mapping[!duplicated(mapping$ENTREZID), ]
        de_df$hgnc_symbol <- mapping$SYMBOL[match(genes, mapping$ENTREZID)]
      } else {
        de_df$hgnc_symbol <- NA
      }
    } else {
      de_df$hgnc_symbol <- genes
    }
  } else if (org == "Mus musculus") {
    n_ens <- sum(grepl("^ENSMUSG\\d+", genes))
    n_entrez <- sum(grepl("^\\d+$", genes))
    if (n_ens > length(genes) * 0.5) {
      valid_ens <- unique(genes[genes %in% keys(org.Mm.eg.db, keytype = "ENSEMBL")])
      if (length(valid_ens) > 0) {
        mapping <- AnnotationDbi::select(org.Mm.eg.db, keys = valid_ens,
                                          columns = c("ENTREZID", "SYMBOL"), keytype = "ENSEMBL")
        mapping <- mapping[!is.na(mapping$ENTREZID), ]
        mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
        de_df$mouse_entrez <- mapping$ENTREZID[match(genes, mapping$ENSEMBL)]
        de_df$mouse_symbol <- mapping$SYMBOL[match(genes, mapping$ENSEMBL)]
        de_df$hgnc_symbol <- ortho_table$human_symbol[match(de_df$mouse_entrez, ortho_table$mouse_entrez)]
        na_idx <- is.na(de_df$hgnc_symbol)
        if (any(na_idx)) {
          de_df$hgnc_symbol[na_idx] <- ortho_table$human_symbol[match(de_df$mouse_symbol[na_idx], ortho_table$mouse_symbol)]
        }
      } else {
        de_df$hgnc_symbol <- NA
      }
    } else if (n_entrez > length(genes) * 0.5) {
      de_df$hgnc_symbol <- ortho_table$human_symbol[match(genes, ortho_table$mouse_entrez)]
      na_idx <- is.na(de_df$hgnc_symbol)
      if (any(na_idx)) {
        valid_entrez <- unique(genes[na_idx][genes[na_idx] %in% keys(org.Mm.eg.db, keytype = "ENTREZID")])
        if (length(valid_entrez) > 0) {
          mapping <- AnnotationDbi::select(org.Mm.eg.db, keys = valid_entrez,
                                            columns = "SYMBOL", keytype = "ENTREZID")
          mapping <- mapping[!is.na(mapping$SYMBOL) & mapping$SYMBOL != "", ]
          mapping <- mapping[!duplicated(mapping$ENTREZID), ]
          mouse_syms <- mapping$SYMBOL[match(genes[na_idx], mapping$ENTREZID)]
          de_df$hgnc_symbol[na_idx] <- ortho_table$human_symbol[match(mouse_syms, ortho_table$mouse_symbol)]
        }
      }
    } else {
      # Mouse gene symbols -> ortholog
      de_df$hgnc_symbol <- ortho_table$human_symbol[match(genes, ortho_table$mouse_symbol)]
    }
  }

  n_before <- nrow(de_df)
  de_df <- de_df[!is.na(de_df$hgnc_symbol) & de_df$hgnc_symbol != "", ]
  de_df$accession <- acc
  de_df$contrast <- ctr
  cat(sprintf("  %s_%s: %d -> %d genes mapped (%.1f%%)\n",
              acc, ctr, n_before, nrow(de_df), 100 * nrow(de_df) / n_before))
  return(de_df)
}

de_list_mapped <- lapply(de_list, map_study_to_human, ortho_table = ortho_table)

# ============================================================
# 4. Prepare effect sizes
# ============================================================
cat("\n=== Preparing effect sizes ===\n")

prepare_effect_sizes <- function(de_df) {
  d <- data.frame(
    gene = de_df$hgnc_symbol,
    log2FC = as.numeric(de_df$log2FC),
    pvalue = as.numeric(de_df$pvalue),
    padj = as.numeric(de_df$padj),
    accession = de_df$accession,
    contrast = de_df$contrast,
    organism = de_df$organism,
    stringsAsFactors = FALSE
  )
  d$pvalue <- pmax(pmin(d$pvalue, 1 - 1e-10), 1e-10)
  d$z <- qnorm(d$pvalue / 2, lower.tail = FALSE)
  d$se <- abs(d$log2FC) / d$z
  d$se[d$se == 0 | !is.finite(d$se)] <- 10
  d$se <- pmin(d$se, 100)
  d <- d[order(d$padj), ]
  d <- d[!duplicated(d$gene), ]
  return(d)
}

de_effects <- lapply(de_list_mapped, prepare_effect_sizes)

# ============================================================
# 5. Random-effects meta-analysis (vectorized)
# ============================================================
cat("\n=== Running random-effects meta-analysis ===\n")

effect_dt <- rbindlist(lapply(names(de_effects), function(nm) {
  d <- de_effects[[nm]]
  d$study <- nm
  d
}))

all_genes <- unique(effect_dt$gene)
cat(sprintf("  Union of genes: %d across %d studies\n", length(all_genes), length(de_effects)))

# Vectorized random-effects meta-analysis using the DerSimonian-Laird estimator
# For each gene: yi = log2FC, vi = se^2
# Fixed effect: w_i = 1/vi, theta_FE = sum(w_i*yi)/sum(w_i)
# Q = sum(w_i*(yi-theta_FE)^2), df = k-1
# tau^2 = max(0, (Q - df) / (sum(w_i) - sum(w_i^2)/sum(w_i)))
# Random effects: w_i_RE = 1/(vi + tau^2), theta_RE = sum(w_i_RE*yi)/sum(w_i_RE)
# SE_RE = sqrt(1/sum(w_i_RE)), z = theta_RE/SE_RE, p = 2*pnorm(-|z|)

cat("  Computing vectorized REML-equivalent meta-analysis...\n")

# Split by gene and compute
setkey(effect_dt, gene)

# Precompute per-gene statistics using data.table
meta_results <- effect_dt[, {
  # Remove NA values
  valid <- !is.na(log2FC) & !is.na(se) & is.finite(log2FC) & is.finite(se) & se > 0
  yi <- log2FC[valid]
  vi <- se[valid]^2
  pv <- pvalue[valid]
  k <- length(yi)
  
  if (k == 0) {
    list(combined_log2FC = NA_real_, combined_se = NA_real_,
         combined_pvalue = 1, n_studies = 0L, I2 = NA_real_)
  } else if (k == 1) {
    list(combined_log2FC = yi[1], combined_se = sqrt(vi[1]),
         combined_pvalue = pv[1], n_studies = 1L, I2 = NA_real_)
  } else {
    # Fixed effect
    w <- 1 / vi
    theta_fe <- sum(w * yi) / sum(w)
    Q <- sum(w * (yi - theta_fe)^2)
    df <- k - 1
    
    # DerSimonian-Laird tau^2
    tau2 <- 0
    if (df > 0 && Q > df && is.finite(Q)) {
      denom <- sum(w) - sum(w^2) / sum(w)
      if (denom > 0) tau2 <- (Q - df) / denom
    }
    
    # Random effects
    w_re <- 1 / (vi + tau2)
    theta_re <- sum(w_re * yi) / sum(w_re)
    se_re <- sqrt(1 / sum(w_re))
    z <- theta_re / se_re
    pval <- 2 * pnorm(-abs(z))
    
    # I^2
    I2_val <- if (df > 0 && Q > df && is.finite(Q)) max(0, 100 * (Q - df) / Q) else 0
    
    list(combined_log2FC = theta_re, combined_se = se_re,
         combined_pvalue = pval, n_studies = k, I2 = I2_val)
  }
}, by = gene]

meta_results$combined_padj <- p.adjust(meta_results$combined_pvalue, method = "BH")
meta_results$direction <- ifelse(meta_results$combined_log2FC > 0, "up", "down")
meta_results <- meta_results[order(combined_padj), ]

# ============================================================
# 6. Define consensus catabolic signature
# ============================================================
cat("\n=== Defining consensus catabolic signature ===\n")

sig_up <- meta_results[combined_log2FC > 0.25 & combined_padj < 0.05, ]
sig_down <- meta_results[combined_log2FC < -0.25 & combined_padj < 0.05, ]

cat(sprintf("  Up-regulated: %d genes\n", nrow(sig_up)))
cat(sprintf("  Down-regulated: %d genes\n", nrow(sig_down)))

write.csv(as.data.frame(meta_results), file.path(signature_dir, "meta_analysis_full.csv"), row.names = FALSE)

signature <- rbind(
  cbind(as.data.frame(sig_up), signature_direction = "up"),
  cbind(as.data.frame(sig_down), signature_direction = "down")
)
write.csv(signature, file.path(signature_dir, "consensus_catabolic_signature.csv"), row.names = FALSE)

ranked_list <- meta_results$combined_log2FC
names(ranked_list) <- meta_results$gene
ranked_list <- sort(ranked_list, decreasing = TRUE)
saveRDS(ranked_list, file.path(signature_dir, "ranked_gene_list.rds"))
write.csv(data.frame(gene = names(ranked_list), log2FC = ranked_list),
          file.path(signature_dir, "ranked_gene_list.csv"), row.names = FALSE)

saveRDS(list(up = sig_up$gene, down = sig_down$gene, full = as.data.frame(meta_results)),
        file.path(signature_dir, "catabolic_signature_gene_sets.rds"))

# ============================================================
# 7. Validate against known atrophy markers
# ============================================================
cat("\n=== Validating against known atrophy markers ===\n")

known_markers <- c(
  "TRIM63", "FBXO32", "MTOR", "RPS6KB1", "RPS6KB2", "EIF4EBP1", "AKT1", "AKT2",
  "FOXO1", "FOXO3", "ULK1", "ATG7", "BECN1", "SQSTM1", "MAP1LC3B",
  "PPARGC1A", "ESRRA", "NDUFS1", "COX5B", "CS",
  "MYH7", "MYH2", "ACTN3", "TTN", "ACTB",
  "IL6", "TNF", "NFKB1", "STAT3", "JAK2",
  "IGF1", "IGF1R", "PIK3CA", "PIK3R1",
  "MYC", "HIF1A", "GSK3B", "TSC1", "TSC2", "RHEB", "DDIT4"
)

known_in_data <- known_markers[known_markers %in% meta_results$gene]
cat(sprintf("\n  %d/%d known markers found\n", length(known_in_data), length(known_markers)))
cat(sprintf("  %-12s %10s %12s %8s %8s\n", "Gene", "log2FC", "padj", "n_studies", "direction"))
for (g in known_in_data) {
  row <- meta_results[gene == g, ]
  cat(sprintf("  %-12s %10.4f %12.2e %8d %8s\n",
              g, row$combined_log2FC, row$combined_padj, row$n_studies, row$direction))
}

cat("\n  Key marker validation:\n")
key_markers <- c("TRIM63", "FBXO32", "MTOR", "RPS6KB1", "FOXO3", "PPARGC1A", "FOXO1", "ULK1")
for (g in key_markers) {
  if (g %in% meta_results$gene) {
    row <- meta_results[gene == g, ]
    cat(sprintf("    %s: log2FC=%.4f, padj=%.2e, n_studies=%d, direction=%s\n",
                g, row$combined_log2FC, row$combined_padj, row$n_studies, row$direction))
  } else {
    cat(sprintf("    %s: NOT FOUND\n", g))
  }
}

cat(sprintf("\n=== Meta-analysis complete ===\n"))
cat(sprintf("  Total genes: %d\n", nrow(meta_results)))
cat(sprintf("  Significant up (|lFC|>0.25, padj<0.05): %d\n", nrow(sig_up)))
cat(sprintf("  Significant down (|lFC|>0.25, padj<0.05): %d\n", nrow(sig_down)))
cat(sprintf("  Results saved to: %s\n", signature_dir))
