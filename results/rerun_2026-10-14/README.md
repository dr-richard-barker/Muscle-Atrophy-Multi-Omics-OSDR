# Re-run after the EDGE fold-change fix (2026-10-14)

`scripts/01_differential_expression.R` built the OSD-530 and GSE213808 tables by taking
`log2()` of EDGE "Fold change" columns. Those columns are signed (≥ 1 up, ≤ −1 down), so
every down-regulated gene became NaN, and the GSE213808 values also pointed the wrong
way. This folder holds the two corrected DE tables and the meta-analysis
(`scripts/02_meta_analysis.R`) re-run on them.

**Nothing downstream has been re-run yet.** `tables/`, the GSEA (03), LINCS (05),
supplement and nutrition outputs, the figures, the manuscript, the report PDF and the
site all still describe the earlier analysis.

## Inputs

| File | Source | Bytes | SHA-256 |
|---|---|---|---|
| `GLDS-530_rna-seq_TGB_050_1_2_64samples_9group_totalcount_all0removed_scalingnormalized_pairwise_analysis_included.xlsx` | NASA OSDR, OSD-530 | 11,499,522 | `d4e7747790711871d428b31ebd2f58b998332cd4e7a45f0316f95165b1fde933` |
| `GSE213808_TGB_022_3group_N6_totalcount_quantilenormalized.xlsx` | NCBI GEO, GSE213808 | 14,454,079 | `2c2e9f7b20882c9d1df74ee05dca973b27cb07ecd4cb0627d3f9e69fbbecb454` |

The other eight DE tables are unchanged and come from `zenodo_package.zip`
(`zenodo_package/results/tables/`).

## What the source files show

- No EDGE fold change lies strictly between −1 and 1 in any comparison of either file.
- An "A vs B" fold change is **B relative to A**. Compared with the group means in the
  same files, its sign agrees with B/A for 79–97% of expressed genes in all seven
  comparisons (OSD-530 Pre vs Flight1–4; GSE213808 MG vs GC, AG vs GC, AG vs MG), and
  for all 16 genes with FDR < 0.05 in those comparisons. So OSD-530 "Pre vs Flight" is
  Flight/Pre, as the pipeline assumed, and GSE213808 "MG vs GC" is GC/MG, which the
  pipeline had labelled flight vs ground without negating.

## Corrected DE tables (`de_tables/`)

| Table | Rows | Negative log2FC | Positive | Zero | NaN | padj < 0.05 |
|---|---|---|---|---|---|---|
| GSE213808 (before) | 49,585 | 0 | 4,954 | 38,173 | 6,458 | 0 |
| GSE213808 (now) | 49,585 | 4,954 | 6,458 | 38,173 | 0 | 0 |
| OSD-530 (before) | 26,845 | 0 | 15,524 | 2,221 | 9,100 | 63 (27 up, 36 NaN) |
| OSD-530 (now) | 26,845 | 17,097 | 7,527 | 2,221 | 0 | 63 (24 up, 39 down) |

P-values and adjusted p-values are identical to before; only log2FC changed. In
GSE213808 every previously finite value is now its negation. The 13 mitochondrially
encoded genes in OSD-530 are unchanged (mean log2FC +1.372): they were positive in all
four comparisons, so the bug did not affect them.

## Meta-analysis

Run with R 4.6.1, readxl 1.5.0.1, dplyr 1.2.1, data.table 1.18.6.1, metafor 5.2.1,
AnnotationDbi 1.74.0, org.Hs.eg.db 3.23.1, org.Mm.eg.db 3.23.0, homologene
1.4.68.19.3.27. The original run's package versions were not recorded.

These packages do not reproduce `tables/` even from the original DE tables, so the two
runs here are compared with each other, not with `tables/`:

| | Genes | Signature | Up | Down | Single-study genes in signature |
|---|---|---|---|---|---|
| `tables/` (original run) | 32,770 | 141 | 77 | 64 | 54 |
| `meta_control/` (original DE tables, these packages) | 33,003 | 132 | 68 | 64 | 44 |
| `meta_fixed/` (corrected DE tables, these packages) | 33,003 | 140 | 60 | 80 | 51 |

- `meta_control` vs `tables/`: 127 signature genes shared. Part of the difference is gene
  symbols: the current annotation names mitochondrially encoded genes `ND1`, `ND2`, `ND4`
  rather than `MT-ND1` and so on.
- `meta_fixed` vs `meta_control`: 116 signature genes kept, 16 dropped, 24 added; combined
  log2FC changes for 22,260 of 33,003 genes and changes sign for 6,365.

## Reproducing

1. Put the two input files above in one folder.
2. In a copy of `scripts/01_differential_expression.R`, set `results_dir`, `f530` and
   `f_geo` to local paths, and run its helper functions plus the OSD-530 and GSE213808
   sections (the other sections need inputs not included here).
3. Put those two tables and the other eight DE tables from `zenodo_package.zip` in one
   folder, set `results_dir` and `signature_dir` in a copy of
   `scripts/02_meta_analysis.R`, and run it.
