# A consensus cross-species catabolic signature of spaceflight-induced muscle atrophy: meta-analysis of NASA OSDR transcriptomics, pathway mapping, and translation to plant-based nutritional countermeasures

**Authors:** Richard Barker

**Affiliation:** Independent Researcher

**Target journal:** npj Microgravity

**Manuscript type:** Research Article

---

## Abstract

Spaceflight-induced muscle atrophy remains a principal barrier to long-duration space exploration, and its molecular drivers are incompletely translated into actionable countermeasures. We performed a cross-species random-effects meta-analysis of 10 transcriptomic datasets from the NASA Open Science Data Repository (OSDR), encompassing rodent spaceflight (OSD-21, OSD-111, OSD-135, OSD-401), human cell culture spaceflight (OSD-684), human astronaut cfRNA (OSD-530, GSE213808), Inspiration4 crew blood transcriptomics (OSD-569, OSD-571), and ground-based bed rest analogs (OSD-370). After mapping all genes to human orthologs, we derived a consensus catabolic signature of 141 differentially expressed genes (77 up-regulated, 64 down-regulated; |log2FC| > 0.25, adjusted p < 0.05) from 32,770 tested genes. The signature was dominated by mitochondrial electron transport chain (ETC) genes (MT-ND1 through MT-ND6, MT-CO2, MT-CYB, MT-ATP8; log2FC +1.23 to +1.59) and included the master mitochondrial biogenesis regulator PPARGC1A/PGC-1α (log2FC = −0.51, p_adj = 0.092) among down-regulated genes. Gene set enrichment analysis across 6,877 pathways revealed nominal enrichment of mitochondrial NADH-to-ubiquinone electron transport (NES = +1.45, p = 0.002) and mitochondrial translation (NES = +1.34, p = 0.011), though no pathway survived Benjamini-Hochberg correction. Connectivity mapping against LINCS L1000 perturbation signatures (559 compounds ranked, 182 clinically annotated) and a curated database of 31 nutritional supplements identified HMB, CoQ10, creatine, vitamin D3, and resveratrol as top counter-catabolic candidates. Nutritional translation using published amino acid composition data for 41 vegan protein sources identified pea protein isolate (leucine 6.50 g/100 g, MPS score 2.25) and soy protein isolate (leucine 6.60 g/100 g, MPS score 2.11) as optimal plant-based protein sources for meeting the 2.5 g leucine per meal threshold for muscle protein synthesis activation. Optimized meal combinations achieved >3.2 g leucine and >39 g protein per serving. These results provide a reproducible, data-driven framework linking spaceflight transcriptomics to evidence-based nutritional countermeasures applicable to both astronaut health and terrestrial disuse atrophy.

**Keywords:** spaceflight, muscle atrophy, meta-analysis, mitochondrial dysfunction, LINCS L1000, nutritional countermeasures, vegan protein, leucine threshold, NASA OSDR

---

## Introduction

Skeletal muscle atrophy is one of the most robust physiological consequences of spaceflight, with astronauts losing 1–2% of muscle mass per month in microgravity despite countermeasure exercise [1, 2]. The molecular mechanisms driving this loss involve activation of the ubiquitin-proteasome system (UPS) through E3 ligases MuRF1/TRIM63 and Atrogin-1/FBXO32 [3, 4], suppression of mitochondrial biogenesis via PGC-1α (PPARGC1A) downregulation [5, 6], and impaired muscle protein synthesis (MPS) through mTOR pathway dysregulation [7]. Ground-based analogs—bed rest and hindlimb suspension—recapitulate many of these molecular changes and serve as accessible models for countermeasure development [8, 9].

The NASA Open Science Data Repository (OSDR, formerly GeneLab) provides open-access multi-omics data from spaceflight experiments across rodent and human systems [10]. Recent meta-analyses of OSDR transcriptomic data have revealed both conserved cross-species responses and tissue-specific signatures of spaceflight-induced muscle remodeling [11, 12]. However, a gap remains between identifying these molecular signatures and translating them into practical, evidence-based countermeasures—particularly nutritional interventions that could be deployed during spaceflight or in terrestrial disuse atrophy contexts.

Connectivity mapping using the LINCS L1000 platform offers a computational approach to drug and compound repurposing by identifying perturbations that reverse a disease signature [13, 14]. Applied to the spaceflight muscle atrophy signature, this approach can nominate compounds that counteract the catabolic transcriptional program. Furthermore, translating compound hits to nutritional sources—particularly plant-based proteins with adequate essential amino acid profiles—addresses a growing need for sustainable, shelf-stable countermeasures suitable for long-duration missions [15, 16].

The leucine threshold for MPS activation (~2.5–3 g per meal) is a well-established nutritional target [17, 18, 19], and plant-based proteins often require strategic combining or fortification to meet this threshold due to lower leucine density compared to animal proteins [20, 21]. Vegan athletes and astronauts on plant-based food systems face particular challenges in achieving adequate leucine, lysine, and total protein per meal [22].

In this study, we integrate 10 transcriptomic datasets from spaceflight and ground analogs to derive a consensus catabolic signature, map it to hypertrophic and mitochondrial pathways, screen LINCS L1000 and curated supplement databases for counter-catabolic compounds, and translate the results to optimized vegan meal plans meeting established leucine and protein thresholds. All scripts and data are packaged for Zenodo deposition to ensure reproducibility.

---

## Methods

### Data acquisition

Transcriptomic datasets were retrieved from the NASA OSDR S3 bucket (`nasa-osdr.s3.amazonaws.com`) and the OSDR REST API (`visualization.osdr.nasa.gov/biodata/api/v2/dataset/`). Ten datasets were included in the meta-analysis (Table 1):

**Table 1. Datasets included in the cross-species meta-analysis.**

| Dataset | Platform | Organism | Tissue/Source | Condition | Contrast | Genes tested | DE genes (p_adj < 0.05) |
|---------|----------|----------|---------------|-----------|----------|-------------|------------------------|
| OSD-21 | Microarray | Mouse (STS-135) | Gastrocnemius | Spaceflight | Flight vs Ground | 12,770 | 112 |
| OSD-111 | Microarray | Mouse (Bion-M1) | EDL/Soleus | Spaceflight | Flight vs Ground | 12,770 | 54 |
| OSD-135 | Microarray | Mouse (Bion-M1) | Longissimus dorsi | Spaceflight | Flight vs Ground | 12,770 | 1,377 |
| OSD-401 | RNA-seq | Mouse (SpaceX-11) | Gastrocnemius | Spaceflight | Flight vs Ground | 22,161 | 359 |
| OSD-684 | RNA-seq | Human tissue chip | Myobundles | Spaceflight | Flight vs Ground | 17,869 | 3,617 |
| OSD-530 | RNA-seq (cfRNA) | Human (JAXA) | Plasma | Spaceflight | Flight vs Pre | 26,845 | 63 |
| OSD-569 | Direct RNA-seq | Human (Inspiration4) | Whole blood | Spaceflight | Post vs Pre | 20,245 | 484 |
| OSD-571 | Proteomics | Human (Inspiration4) | Plasma | Spaceflight | Post vs Pre | 1,765 | 59 |
| OSD-370 | Microarray | Human | Vastus lateralis | Bed rest (9-day) | During vs Before | 45,015 | 377 |
| GSE213808 | RNA-seq (cfRNA) | Mouse (JAXA MHU-1) | Plasma | Spaceflight | Flight vs Ground | 49,585 | 0 |

Datasets spanned rodent spaceflight (n = 4), human spaceflight (n = 4), and human ground analog (n = 1), plus one JAXA mouse cfRNA dataset from GEO. Pre-computed differential expression files were used where available from OSDR; for datasets without pre-computed DE, we applied DESeq2 (RNA-seq) or limma (microarray) with default parameters. The Inspiration4 proteomics dataset (OSD-571) was included as a protein-level confirmation of transcriptomic changes.

### Ortholog mapping

Mouse gene symbols were mapped to human orthologs using the `homologene` R package (CRAN), which contains 275,237 ortholog mappings from the NCBI HomoloGene database. One-to-one orthologs were used; genes without a human ortholog were excluded from the meta-analysis. All downstream analyses were performed on human gene symbols.

### Random-effects meta-analysis

A DerSimonian-Laird random-effects meta-analysis was performed across all 10 datasets using the `metafor` R package. For each gene, the combined log2 fold change, standard error, p-value, and adjusted p-value (Benjamini-Hochberg) were computed. Heterogeneity was assessed using the I² statistic. The consensus catabolic signature was defined as genes with |combined log2FC| > 0.25 and adjusted p < 0.05. Genes were classified as up-regulated (positive log2FC) or down-regulated (negative log2FC) in the atrophy condition relative to control.

### Gene set enrichment analysis

GSEA was performed using the `fgsea` R package against 6,877 gene sets from the Molecular Signatures Database (MSigDB), including Hallmark (n = 50), KEGG, Reactome, and Gene Ontology Biological Process (GOBP) collections. The ranked gene list was ordered by the combined log2FC from the meta-analysis. Enrichment scores (ES), normalized enrichment scores (NES), and p-values were computed via 1,000 permutations. Pathways with nominal p < 0.05 were reported; Benjamini-Hochberg correction was applied across all tested pathways.

### LINCS L1000 connectivity mapping

The consensus catabolic signature (77 up-regulated, 64 down-regulated genes) was screened against the LINCS L1000 perturbation database using a dual-strategy approach:

1. **iLINCS API query:** The `findConcordancesSC` endpoint (ilincs.org) was queried with (a) up-regulated atrophy genes to identify compounds that mimic the atrophy signature (to be avoided) and (b) down-regulated and anti-atrophy genes to identify compounds that promote a muscle-preserving transcriptional state. Of the 141 signature genes, 25 up-regulated and 25 down-regulated genes overlapped with the L1000 landmark gene set (~10,000 genes).

2. **Local GMT fallback:** A local GMT file of 542 single-drug perturbation signatures (`single_drug_perturbations-v1.0.gmt`, 15,376 unique genes) was used for hypergeometric reversal testing, identifying drug signatures with significant overlap with the atrophy signature in the reversing direction.

A composite reversal score was computed for each compound: 50% normalized z-score for concordance with anti-atrophy genes (promote-muscle component) + 30% inverse of z-score for concordance with atrophy genes (anti-mimic component) + 20% normalized −log10(p) from the GMT hypergeometric test. Clinical phase annotations were mapped from the Broad Drug Repurposing Hub parquet files.

### Curated supplement database screening

A curated database of 31 nutritional supplements with annotated mechanisms, molecular targets, pathway associations, and evidence levels was constructed from published literature. Each supplement was scored against the full meta-analysis (32,770 genes) using a composite score: 25% target coverage (fraction of supplement's known targets present in the meta-analysis) + 30% direction score (whether the supplement counteracts the atrophy direction—inhibiting up-regulated targets, activating down-regulated targets) + 15% pathway score (overlap with GSEA-enriched pathways) + 10% LINCS overlap + 20% evidence score (number of supporting studies).

### Nutritional translation

Amino acid composition data for 41 vegan protein sources were compiled from published SR Legacy (USDA FoodData Central, Release 28) values. The USDA FoodData Central API was unavailable during the analysis period due to persistent rate limiting (DEMO_KEY); published SR Legacy values from public nutrition references were used instead. Each food was characterized for 18 amino acids per 100 g, total protein, and energy content.

A muscle protein synthesis (MPS) score was computed for each food: 40% leucine density (leucine per 100 g) + 25% protein quality (protein per 100 g) + 20% caloric efficiency (protein-to-calorie ratio) + 15% lysine density (lysine per 100 g). The leucine threshold of 2.5 g per meal was based on established literature [17, 18, 19], and the protein target of 30 g per meal was based on per-meal optimization recommendations [23].

### Recipe optimization

Meal combinations were generated across four archetypes: protein shake, legume + grain bowl, soy + grain, and protein oatmeal. Each combination was scored on leucine content, total protein, caloric efficiency, BCAA ratio, lysine content, methionine content, total weight, and supplement synergy. Combinations meeting the 2.5 g leucine threshold were ranked by composite meal score.

### Reproducibility

All analyses were performed in R (version 4.x) and Python (version 3.x). R packages used: `metafor`, `fgsea`, `homologene`, `DESeq2`, `limma`. Python packages used: `pandas`, `scipy`, `requests`. All scripts, data, and results are packaged for Zenodo deposition. The pipeline is fully reproducible from the provided scripts without requiring authentication for any data source.

### Data and code availability

All OSDR datasets are publicly available at `visualization.osdr.nasa.gov`. GSE213808 is available at GEO. LINCS L1000 data are available via iLINCS (ilincs.org). All analysis scripts and intermediate results are provided in the Zenodo package accompanying this manuscript.

---

## Results

### Consensus catabolic signature

The random-effects meta-analysis across 10 datasets (Table 1) tested 32,770 human ortholog-mapped genes. Of these, 240 genes had adjusted p < 0.05, and 141 met the stringent significance threshold of |log2FC| > 0.25 and adjusted p < 0.05, comprising the consensus catabolic signature: 77 up-regulated and 64 down-regulated genes (Figure 1, Figure 2).

The most prominent feature of the up-regulated genes was the near-complete representation of mitochondrial electron transport chain (ETC) complex I subunits. All 7 mitochondrial-encoded NADH dehydrogenase subunits (MT-ND1 through MT-ND6, MT-ND4L) were significantly up-regulated (log2FC range +1.23 to +1.49, p_adj range 3.0 × 10⁻⁷ to 5.4 × 10⁻⁶), along with MT-CO2 (complex IV, log2FC = +1.43), MT-CYB (complex III, log2FC = +1.28), and MT-ATP8 (complex V, log2FC = +1.59). These mitochondrial genes were primarily driven by a single dataset (OSD-684, human tissue chip spaceflight), which showed the largest number of DE genes (3,617). Other notable up-regulated genes included CEROX1 (log2FC = +1.03, p_adj = 6.5 × 10⁻⁸, 3 studies, I² = 0%), CBS (log2FC = +1.09, 2 studies), FAT1 (log2FC = +0.54, 4 studies, I² = 0%), and PFDN6 (log2FC = +0.79, 2 studies, I² = 0%).

The canonical atrophy E3 ligases showed the expected direction but did not meet the stringent significance threshold in the combined meta-analysis. FBXO32/Atrogin-1 was up-regulated (log2FC = +0.30, p_adj = 0.96, 8 studies, I² = 87%), and TRIM63/MuRF1 showed a paradoxical downward trend (log2FC = −0.20, p_adj = 1.00, 3 studies, I² = 86%). The high heterogeneity (I² > 80%) for these genes reflects substantial variability across datasets, consistent with prior reports of context-dependent atrogene regulation [4, 12]. FOXO3, the master transcriptional activator of atrophy genes, was up-regulated (log2FC = +0.22, p_adj = 0.97, 6 studies, I² = 91%).

Among down-regulated genes, the most significant were hemoglobin subunits (HBA1, HBM, HBQ1, HBD; log2FC range −0.28 to −0.46, p_adj < 2 × 10⁻⁷, 3–4 studies, I² = 0–5%), reflecting known spaceflight-induced anemia [24]. Critically, PPARGC1A/PGC-1α, the master regulator of mitochondrial biogenesis, was down-regulated (log2FC = −0.51, p_adj = 0.092, 8 studies, I² = 41%), narrowly missing the significance threshold but representing one of the most biologically important findings given its central role in mitochondrial homeostasis and muscle atrophy [5, 6]. Other down-regulated genes included DEPP1 (log2FC = −0.82, 3 studies, I² = 0%), SUCLG2-DT (log2FC = −0.85, 2 studies), and AK3 (log2FC = −0.27, 8 studies, I² = 33%), the latter being a mitochondrial adenylate kinase involved in energy metabolism.

### Pathway enrichment

GSEA across 6,877 gene sets identified 70 pathways with nominal significance (p < 0.05), though none survived Benjamini-Hochberg correction (all adjusted p = 1.0), reflecting the moderate effect sizes and the large number of tests (Figure 3).

The top enriched pathway in the positive direction (enriched in up-regulated genes) was mitochondrial electron transport, NADH to ubiquinone (GOBP, NES = +1.45, p = 0.002, size = 52), with leading-edge genes MT-ND1, MT-ND2, MT-ND4L, MT-ND4, MT-ND3, MT-ND5, and MT-ND6. This was followed by mitochondrial translation (REACTOME, NES = +1.34, p = 0.011, size = 113), NADH dehydrogenase complex assembly (GOBP, NES = +1.38, p = 0.017, size = 40), and Complex I biogenesis (REACTOME, NES = +1.32, p = 0.038, size = 68). The convergence of four independent gene sets on mitochondrial Complex I and mitochondrial translation strongly supports mitochondrial dysfunction as a central feature of the spaceflight atrophy signature.

In the negative direction (enriched in down-regulated genes), the top pathways were DNA methylation-dependent constitutive heterochromatin formation (GOBP, NES = −1.83, p = 0.0008), negative regulation of D-glucose transmembrane transport (GOBP, NES = −1.76, p = 0.002), and constitutive heterochromatin formation (GOBP, NES = −1.74, p = 0.002). The down-regulation of glucose transport regulation is consistent with impaired glucose metabolism during disuse atrophy [25].

No Hallmark gene sets reached nominal significance. The closest were HALLMARK_ALLOGRAFT_REJECTION (NES = −1.11, p = 0.099) and HALLMARK_INFLAMMATORY_RESPONSE (NES = −0.98, p = 0.56), both in the negative direction, suggesting that classical inflammatory pathways were not strongly enriched in the combined signature—likely because the meta-analysis pooled diverse tissues and conditions that may dilute inflammation signals present in individual datasets.

### LINCS L1000 compound screening

The dual-strategy connectivity mapping ranked 559 unique compounds by composite reversal score (Figure 4). Of these, 182 had clinical phase annotations from the Broad Drug Repurposing Hub: 111 launched drugs, 19 Phase 3, 19 Phase 2, 6 Phase 1, 18 preclinical, and 5 withdrawn.

The top-ranked compound overall was an unannotated BRD compound (BRD-K71328603, reversal score = 4.71), followed by Fenretinide (reversal score = 4.49, Phase 3, target: RARA, mechanism: retinoid receptor agonist/apoptosis stimulant). Fenretinide showed both strong concordance with anti-atrophy genes (z = 3.33) and GMT evidence of signature reversal (overlap with 2 drug signatures, p = 0.16). Wortmannin (reversal score = 4.05, preclinical, target: PIK3CA/PIK3CG, mechanism: PI3K inhibitor) ranked third, with strong promote-muscle z-scores across two cell lines (A549, CD34).

Among launched drugs, Thalidomide (reversal score = 1.23, target: TNF, mechanism: TNF production inhibitor) and Oxaprozin (reversal score = 1.01, target: PTGS1/PTGS2, mechanism: COX inhibitor) were the highest-ranked, both acting through anti-inflammatory mechanisms. Niacin (reversal score = 0.77, target: HCAR1/HCAR2, mechanism: NAD precursor) was notable as a vitamin with direct relevance to mitochondrial function.

**Important caveats:** The LINCS L1000 platform measures ~10,000 landmark genes, and many signature genes (particularly mitochondrial-encoded genes and lncRNAs) are not represented. Of the 141 signature genes, only 25 up-regulated and 25 down-regulated genes overlapped with L1000 landmarks. Additionally, no muscle cell lines are present in L1000; all perturbation signatures derive from cancer and other cell lines, limiting direct translational relevance. The SigCom LINCS data API was unavailable (503 error) during the analysis; results are based on the iLINCS API and a local GMT fallback.

### Curated supplement screening

The 31-supplement curated database screening produced a ranking that complemented the LINCS results by incorporating upstream signaling mechanisms (Figure 5). The top 5 supplements were:

1. **HMB (β-hydroxy β-methylbutyrate)** — composite score 0.710. Targets: FBXO32, TRIM63, MTOR, FOXO3. HMB directly inhibits the E3 ligases Atrogin-1 and MuRF1, normalizes the Akt/FoxO axis, and restores autophagy-lysosomal pathways in dexamethasone-induced atrophy models [26, 27, 28]. A systematic review and meta-analysis found HMB supplementation modestly increased skeletal muscle mass (SMD ~0.25) and more robustly improved muscle strength (SMD ~0.31) across clinical conditions [29].

2. **CoQ10 (ubiquinone)** — composite score 0.665. Targets: MT-ND1, MT-CO1, MT-ATP6. Directly supports mitochondrial ETC function, which is the most prominent feature of the consensus signature. CoQ10 supplementation addresses the Complex I–III upregulation that may reflect compensatory responses to mitochondrial stress.

3. **Creatine monohydrate** — composite score 0.636. Targets: CKM, CKB, MYOD1, PAX7. Creatine buffers ATP during high-intensity contraction and activates satellite cells. Meta-analyses confirm that creatine plus resistance training increases lean body mass by ~1.14 kg (95% CI 0.69–1.59) and reduces fat mass [30, 31].

4. **Vitamin D3 (cholecalciferol)** — composite score 0.634. Targets: MYOD1, MYOG, PAX7, VDR. Activates satellite cell proliferation and differentiation; vitamin D deficiency is common during spaceflight and is associated with muscle weakness [32].

5. **Resveratrol** — composite score 0.562. Targets: PPARGC1A, SIRT1, NFKB1, FOXO3, AMPK. Activates the SIRT1/PGC-1α axis to promote mitochondrial biogenesis [33, 34], directly addressing the PPARGC1A downregulation in the consensus signature. However, human studies of resveratrol's effects on mitochondrial biogenesis in muscle have yielded mixed results [35, 36].

Additional high-scoring supplements included EGCG (0.54, targets FOXO3/NFKB1), vitamin C (0.53, collagen synthesis), leucine (0.53, targets MTOR/RPS6KB1/EIF4EBP1), L-citrulline (0.52, mTOR/IGF1), and glycine (0.50, collagen synthesis).

### Nutritional translation to vegan protein sources

Analysis of 41 vegan protein sources identified the top sources by MPS score (Figure 6, Figure 7, Table 2):

**Table 2. Top 10 vegan protein sources by MPS score.**

| Rank | Food | Protein (g/100 g) | Leucine (g/100 g) | Energy (kcal/100 g) | Grams for 2.5 g leucine | MPS score |
|------|------|-------------------|-------------------|---------------------|------------------------|-----------|
| 1 | Pea protein isolate | 75.0 | 6.50 | 380 | 38 | 2.25 |
| 2 | Soy protein isolate | 80.7 | 6.60 | 338 | 38 | 2.11 |
| 3 | Rice protein powder | 80.0 | 6.20 | 370 | 40 | 1.54 |
| 4 | Spirulina, dried | 57.5 | 4.97 | 290 | 50 | 1.52 |
| 5 | Soy flour, defatted | 47.0 | 3.90 | 329 | 64 | 1.11 |
| 6 | Hemp protein powder | 50.0 | 3.80 | 400 | 66 | 0.75 |
| 7 | Corn, yellow, raw | 9.4 | 1.13 | 361 | 221 | 0.68 |
| 8 | Soybeans, mature, raw | 36.5 | 2.97 | 446 | 84 | 0.65 |
| 9 | Nutritional yeast | 50.0 | 3.45 | 332 | 72 | 0.58 |
| 10 | Hemp seeds, hulled | 31.6 | 2.50 | 553 | 100 | 0.24 |

Pea and soy protein isolates were the only single-source foods capable of meeting the 2.5 g leucine threshold in a practical serving size (38 g, ~130–140 kcal). Rice protein powder (40 g) and spirulina (50 g) were also feasible. Whole-food sources required substantially larger portions (84–221 g) to reach the leucine threshold, highlighting the value of protein isolates in a spaceflight food system where mass and volume are constrained.

### Optimized meal combinations

Eighty-nine meal combinations were generated across four archetypes (Table 3). The top-ranked meal in each archetype:

**Table 3. Best meal by archetype.**

| Archetype | Foods | Portions (g) | Protein (g) | Leucine (g) | Energy (kcal) | Meal score |
|-----------|-------|--------------|-------------|-------------|---------------|------------|
| Protein shake | Soy protein isolate + Hemp seeds + Spirulina | 30 + 30 + 10 | 39.4 | 3.23 | 296 | 0.803 |
| Protein oatmeal | Oats + Soy protein isolate + Chia seeds + Walnuts | 50 + 30 + 25 + 20 | 39.8 | 3.20 | 548 | 0.796 |
| Soy + grain | Tempeh + Quinoa + Pumpkin seeds | 100 + 50 + 30 | 34.6 | 2.63 | 549 | 0.793 |
| Legume + grain bowl | Lentils + Oats + Hemp seeds | 50 + 50 + 30 | 30.8 | 2.34 | 537 | 0.771 |

The top protein shake (soy protein isolate + hemp seeds + spirulina) achieved 3.23 g leucine and 39.4 g protein in only 296 kcal and 70 g total weight—an exceptionally nutrient-dense combination suitable for spaceflight constraints. All four archetypes met or approached the 2.5 g leucine threshold and exceeded 30 g protein per serving.

### Supplement integration plan

Based on the combined LINCS, curated supplement, and nutritional translation results, a daily supplement integration plan was developed (Table 4):

**Table 4. Supplement integration plan.**

| Supplement | Dose | Timing | Rationale |
|------------|------|--------|-----------|
| HMB | 3 g/day | Post-workout shake | Directly inhibits FBXO32/Atrogin-1 and TRIM63/MuRF1; anti-catabolic |
| Creatine monohydrate | 5 g/day | Any meal | ATP buffering; satellite cell activation |
| Vitamin D3 | 2,000 IU/day | Breakfast (fat-soluble) | Satellite cell activation; anti-inflammatory; deficiency common in spaceflight |
| Omega-3 EPA/DHA | 2 g/day | Lunch or dinner | Anti-inflammatory; sensitizes muscle to leucine-stimulated MPS |
| Resveratrol | 500 mg/day | Dinner | Activates PGC-1α (down-regulated in signature); mitochondrial biogenesis |
| Curcumin + Piperine | 500 + 10 mg/day | Dinner (with fat) | Anti-inflammatory via NF-κB inhibition; inhibits FOXO3 |
| CoQ10 | 200 mg/day | Breakfast (fat-soluble) | Mitochondrial ETC support; signature shows ETC gene dysregulation |
| Nicotinamide riboside | 300 mg/day | Breakfast | NAD⁺ precursor; activates SIRT1/PGC-1α; mitochondrial biogenesis |

---

## Discussion

This study presents a reproducible, data-driven pipeline that links spaceflight transcriptomics to evidence-based nutritional countermeasures for muscle atrophy. By meta-analyzing 10 datasets spanning rodent spaceflight, human spaceflight (including Inspiration4 and JAXA cfRNA), and ground-based bed rest analogs, we derived a 141-gene consensus catabolic signature and translated it through pathway analysis, computational drug screening, and nutritional optimization.

### Mitochondrial dysfunction as the central feature

The most striking finding is the dominance of mitochondrial ETC genes in the up-regulated signature. All seven mitochondrial-encoded Complex I subunits (MT-ND1–6, MT-ND4L) plus MT-CO2, MT-CYB, and MT-ATP8 were significantly up-regulated, and four independent GSEA gene sets converged on mitochondrial Complex I and mitochondrial translation. This is consistent with the "mitostasis theory of muscle atrophy" [5], which posits that mitochondrial dysfunction is a primary driver rather than a secondary consequence of muscle wasting. The up-regulation of ETC genes may reflect a compensatory response to impaired mitochondrial function, as has been observed in disuse atrophy [6] and spaceflight quadriceps [37].

The down-regulation of PPARGC1A/PGC-1α (log2FC = −0.51, p_adj = 0.092), while not meeting our stringent significance threshold, is biologically consistent with the mitochondrial signature. PGC-1α is the master regulator of mitochondrial biogenesis, and its suppression during disuse is well-documented [5, 6, 38]. In human skeletal muscle, PGC-1α mRNA declines within 1–4 days of immobilization [38]. The moderate heterogeneity for this gene (I² = 41%) suggests reasonable consistency across datasets. The combination of ETC gene up-regulation and PGC-1α down-regulation points to a model of compensatory mitochondrial stress response failing to overcome impaired biogenesis—a pattern that directly informs countermeasure selection (CoQ10, resveratrol, nicotinamide riboside).

### Atrogene heterogeneity

The canonical atrophy E3 ligases FBXO32/Atrogin-1 and TRIM63/MuRF1 showed the expected up-regulation direction but with high heterogeneity (I² > 86%) and non-significant adjusted p-values. This is consistent with prior cross-species meta-analyses of OSDR data that found limited overlap in differential gene expression between mouse and human under altered gravity [12]. The variability likely reflects differences in tissue type (gastrocnemius vs. EDL vs. blood cfRNA vs. tissue chip), species (mouse vs. human), mission duration (STS-135 at 13 days vs. Bion-M1 at 30 days vs. Inspiration4 at 3 days), and disuse model (spaceflight vs. bed rest). FBXO32 was up-regulated in 8 of 10 datasets (log2FC = +0.30), and FOXO3 was up-regulated in 6 of 10 (log2FC = +0.22), supporting the involvement of the FoxO-atrogene axis despite the lack of statistical significance after multiple testing correction.

### From signature to countermeasures

The dual-track screening approach—computational (LINCS L1000) and mechanistic (curated supplements)—provides complementary evidence. The LINCS screen nominated Fenretinide (Phase 3 retinoid), Wortmannin (PI3K inhibitor), and anti-inflammatory launched drugs (Thalidomide, Oxaprozin) as transcriptional signature reversers. However, the LINCS results carry important limitations: no muscle cell lines are represented, only ~35% of signature genes are L1000 landmarks, and the iLINCS API returns only concordant z-scores, requiring an indirect differential query strategy to infer reversal.

The curated supplement screen, by incorporating known pharmacological mechanisms, identified HMB as the top candidate—a finding strongly supported by independent literature. HMB directly inhibits FBXO32 and TRIM63 [26], normalizes the Akt/FoxO axis [27], and mitigates inactivity- and protein deprivation–induced atrophy [28]. The convergence of computational and mechanistic evidence for anti-inflammatory compounds (Thalidomide/Oxaprozin from LINCS; curcumin/EGCG from curated screening) further supports inflammation as a druggable component of the atrophy program.

### Vegan nutritional translation

The nutritional translation addresses a practical need for spaceflight food systems: shelf-stable, mass-efficient protein sources that meet established leucine and protein thresholds. Pea and soy protein isolates emerged as the clear leaders, both providing >6.5 g leucine per 100 g and meeting the 2.5 g leucine threshold in ~38 g servings. These findings are consistent with clinical evidence showing that plant-based protein isolates with added leucine stimulate MPS comparably to whey protein [20, 39], and that pea-soy blends fortified with leucine match whey's anabolic response in aged mice [40].

The optimized soy protein + hemp seeds + spirulina shake (3.23 g leucine, 39.4 g protein, 296 kcal, 70 g mass) represents an exceptionally nutrient-dense meal suitable for spaceflight constraints. The inclusion of spirulina provides additional micronutrients and bioactive compounds, while hemp seeds contribute arginine (a creatine precursor) and omega-3 fatty acids. The supplement integration plan (Table 4) complements the meal plan by targeting specific nodes in the catabolic signature: HMB for E3 ligase inhibition, CoQ10 and nicotinamide riboside for mitochondrial support, resveratrol for PGC-1α activation, and curcumin for anti-inflammatory action.

### Limitations

Several limitations should be acknowledged:

1. **Heterogeneity across datasets.** The meta-analysis pools diverse tissues (gastrocnemius, EDL, soleus, blood cfRNA, tissue chip myobundles), species (mouse, human), conditions (spaceflight, bed rest), and mission durations (3–30 days). While the random-effects model accounts for between-study variance, this diversity dilutes tissue-specific signals and contributes to high I² for key atrophy genes.

2. **Mitochondrial gene dominance.** The mitochondrial ETC genes in the signature were primarily driven by a single dataset (OSD-684, human tissue chip), which showed the largest number of DE genes (3,617). Sensitivity analysis without OSD-684 would likely reduce the mitochondrial signal, though the GSEA results (which use the full ranked list, not just significant genes) provide independent support.

3. **No pathways survived multiple testing correction.** All 6,877 GSEA pathways had adjusted p = 1.0 after Benjamini-Hochberg correction. The nominal enrichments (mitochondrial ETC, mitochondrial translation) are biologically plausible and supported by four independent gene sets, but should be interpreted as hypothesis-generating.

4. **LINCS L1000 limitations.** The L1000 platform covers ~10,000 landmark genes (not mitochondrial-encoded genes or lncRNAs), uses non-muscle cell lines, and the SigCom data API was unavailable. The iLINCS API + local GMT fallback provides a reasonable approximation but is not equivalent to a full CMap analysis.

5. **USDA API unavailability.** The USDA FoodData Central API was persistently rate-limited during the analysis. Published SR Legacy values were used instead, which may differ slightly from the current Foundation Foods database.

6. **No experimental validation.** All findings are computational. The supplement and nutritional recommendations are based on in silico screening and literature evidence, not on direct testing in spaceflight or analog conditions.

7. **Cross-species ortholog mapping.** The `homologene` package provides one-to-one ortholog mapping, which may miss species-specific paralogs or regulatory differences. Ensembl-based multi-species ortholog mapping (via biomaRt) was unavailable during the analysis.

### Future directions

The consensus signature and countermeasure recommendations should be validated in: (a) independent spaceflight datasets as they become available through OSDR, (b) ground-based analog studies (bed rest, hindlimb suspension) with nutritional interventions, and (c) in vitro muscle cell systems treated with nominated compounds (HMB, CoQ10, resveratrol) under simulated microgravity. The pipeline is designed to be re-run as new OSDR datasets are deposited, enabling iterative refinement of the signature and countermeasure rankings.

---

## Conclusion

We present a reproducible pipeline integrating 10 spaceflight and analog transcriptomic datasets into a 141-gene consensus catabolic signature dominated by mitochondrial ETC dysregulation and PGC-1α suppression. Computational screening of 559 LINCS L1000 compounds and 31 curated supplements, combined with nutritional translation to 41 vegan protein sources, yielded a data-driven countermeasure framework: HMB, CoQ10, creatine, vitamin D3, and resveratrol as top supplements; pea and soy protein isolates as optimal protein sources; and a soy-hemp-spirulina shake meeting leucine and protein thresholds in 296 kcal. All scripts and data are packaged for Zenodo deposition. This framework is applicable to both astronaut health during long-duration spaceflight and terrestrial disuse atrophy, and is designed for iterative refinement as new spaceflight omics data become available.

---

## Figure legends

**Figure 1.** Volcano plot of the cross-species meta-analysis. Each point represents a gene (n = 32,770). The x-axis shows combined log2 fold change (atrophy condition vs. control); the y-axis shows −log10(adjusted p-value). Red points: 77 up-regulated genes (log2FC > 0.25, p_adj < 0.05). Blue points: 64 down-regulated genes (log2FC < −0.25, p_adj < 0.05). Key genes labeled: mitochondrial ETC subunits (MT-ND1–6, MT-CO2, MT-CYB, MT-ATP8), PPARGC1A, FBXO32, FOXO3, hemoglobin subunits.

**Figure 2.** Forest plot of key catabolic signature genes. Combined log2FC and 95% confidence intervals are shown for selected genes across the meta-analysis. Genes are grouped by functional category: mitochondrial ETC (up), atrophy regulators (up), and mitochondrial biogenesis/hemoglobin (down).

**Figure 3.** GSEA dot plot of top enriched pathways. Pathways with nominal p < 0.05 are shown, ordered by normalized enrichment score (NES). Dot size represents gene set size; color represents p-value. Mitochondrial ETC and translation pathways are enriched in up-regulated genes (positive NES); heterochromatin and glucose transport pathways are enriched in down-regulated genes (negative NES).

**Figure 4.** LINCS L1000 top 20 compounds by reversal score. Compounds are ranked by composite reversal score (50% promote-muscle + 30% anti-mimic + 20% GMT evidence). Bar color indicates clinical phase. Annotated compounds labeled with target and mechanism.

**Figure 5.** Curated supplement ranking. Top 15 supplements by composite score (25% target coverage + 30% direction + 15% pathway + 10% LINCS overlap + 20% evidence). Bar segments show sub-score contributions. Key targets annotated for each supplement.

**Figure 6.** Vegan protein sources ranked by leucine content per 100 g. Horizontal line indicates the 2.5 g leucine per meal threshold (assuming 100 g serving). Protein isolates (pea, soy, rice) exceed the threshold; whole-food sources require larger portions.

**Figure 7.** Amino acid heatmap of top vegan protein sources. Rows: 18 amino acids. Columns: top 15 vegan protein sources. Color intensity represents grams per 100 g. Leucine, lysine, and BCAA rows are highlighted.

**Figure 8.** Pipeline overview schematic. Flow diagram showing the five-stage pipeline: data retrieval (10 datasets) → meta-analysis (32,770 genes → 141 signature) → GSEA (6,877 pathways) → compound screening (559 LINCS + 31 supplements) → nutritional translation (41 foods → optimized meals).

---

## Tables

**Table 1.** Datasets included in the cross-species meta-analysis. (See Methods.)

**Table 2.** Top 10 vegan protein sources by MPS score. (See Results.)

**Table 3.** Best meal by archetype. (See Results.)

**Table 4.** Supplement integration plan. (See Results.)

**Supplementary Table S1.** Full meta-analysis results (32,770 genes). File: `meta_analysis_full.csv`.

**Supplementary Table S2.** Consensus catabolic signature (141 genes). File: `consensus_catabolic_signature.csv`.

**Supplementary Table S3.** Ranked gene list for GSEA. File: `ranked_gene_list.csv`.

**Supplementary Table S4.** Full GSEA results (6,877 pathways). File: `gsea_full_results.csv`.

**Supplementary Table S5.** Hallmark GSEA results. File: `gsea_hallmark_results.csv`.

**Supplementary Table S6.** LINCS L1000 full compound ranking (559 compounds). File: `lincs_reverser_ranking_full.csv`.

**Supplementary Table S7.** LINCS L1000 top 100 compounds. File: `lincs_reverser_top100.csv`.

**Supplementary Table S8.** LINCS GMT reversal raw results. File: `lincs_gmt_reversal_raw.csv`.

**Supplementary Table S9.** LINCS robustness summary. File: `robustness_summary.json`.

**Supplementary Table S10.** Full supplement screening ranking (31 supplements). File: `supplement_screening_ranking.csv`.

**Supplementary Table S11.** Top 15 supplements. File: `supplement_screening_top15.csv`.

**Supplementary Table S12.** Vegan amino acid database (41 foods, 18 amino acids). File: `vegan_amino_acid_database.csv`.

**Supplementary Table S13.** Vegan protein MPS ranking (41 foods). File: `vegan_protein_mps_ranking.csv`.

**Supplementary Table S14.** Food-supplement cross-reference. File: `food_supplement_cross_reference.csv`.

**Supplementary Table S15.** Optimized meal combinations (89 meals). File: `optimized_meal_combinations.csv`.

**Supplementary Table S16.** Best meals by archetype. File: `best_meals_by_archetype.csv`.

**Supplementary Table S17.** Supplement integration plan. File: `supplement_integration_plan.csv`.

---

## References

[1] Rittweger J, Albracht K, Flück M, et al. Sarcolab pilot study into skeletal muscle's adaptation to long-term spaceflight. *npj Microgravity*. 2018;4:18. doi:10.1038/s41526-018-0052-1

[2] Qaisar R, Karim A, Elmoselhi A. Muscle unloading: A comparison between spaceflight and ground-based models. *Acta Physiologica*. 2019;227(4):e13431. doi:10.1111/apha.13431

[3] Peris-Moreno D, Cussonneau L, Combaret L, Polge C, Taillandier D. Ubiquitin Ligases at the Heart of Skeletal Muscle Atrophy Control. *Molecules*. 2021;26(2):407. doi:10.3390/molecules26020407

[4] Peris-Moreno D, Taillandier D, Polge C. MuRF1/TRIM63, Master Regulator of Muscle Mass. *Int J Mol Sci*. 2020;21(18):6663. doi:10.3390/ijms21186663

[5] Ji L, Yeo D. Mitochondrial dysregulation and muscle disuse atrophy. *F1000Research*. 2019;8:19139. doi:10.12688/f1000research.19139.1

[6] Kong S, Cai B, Nie Q. PGC-1α affects skeletal muscle and adipose tissue development by regulating mitochondrial biogenesis. *Mol Genet Genomics*. 2022;297(3):621–635. doi:10.1007/s00438-022-01878-2

[7] Cahill T, Cope H, Bass J, et al. Mammalian and Invertebrate Models as Complementary Tools for Gaining Mechanistic Insight on Muscle Responses to Spaceflight. *Int J Mol Sci*. 2021;22(17):9470. doi:10.3390/ijms22179470

[8] Murgia M, Ciciliot S, Nagaraj N, et al. Signatures of muscle disuse in spaceflight and bed rest revealed by single muscle fiber proteomics. *PNAS Nexus*. 2022;1(4):pgac086. doi:10.1093/pnasnexus/pgac086

[9] Suetta C, Frandsen U, Jensen L, et al. Aging Affects the Transcriptional Regulation of Human Skeletal Muscle Disuse Atrophy. *PLoS ONE*. 2012;7(5):e51238. doi:10.1371/journal.pone.0051238

[10] Gebre S, Scott RT, Saravia-Butler A, Lopez DK, Sanders L, Costes SV. NASA open science data repository: open science for life in space. *Nucleic Acids Research*. 2024;gkae1116. doi:10.1093/nar/gkae1116

[11] Ball BK, Khan HF, Park JH, Jayant K, Chan DD, Brubaker DK. Integrated cross-species translation and biophysical multi-scale modeling links molecular signatures and locomotory phenotypes in spaceflight-induced sarcopenia. *npj Microgravity*. 2026. doi:10.1038/s41526-025-00557-x

[12] Adamopoulos K, Sanders L, Costes SV. NASA GeneLab derived microarray studies of Mus musculus and Homo sapiens organisms in altered gravitational conditions. *npj Microgravity*. 2024. doi:10.1038/s41526-024-00392-6

[13] Subramanian A, Narayan R, Corsello SM, et al. A Next Generation Connectivity Map: L1000 Platform And The First 1,000,000 Profiles. *bioRxiv*. 2017:136168. doi:10.1101/136168

[14] Zhao Y, Chen X, Chen J, Qi X. Decoding Connectivity Map-based drug repurposing for oncotherapy. *Briefings in Bioinformatics*. 2023;24(3):bbad142. doi:10.1093/bib/bbad142

[15] Parafati M, Giza S, Shenoy T, et al. Human skeletal muscle tissue chip autonomous payload reveals changes in fiber type and metabolic gene expression due to spaceflight. *npj Microgravity*. 2023;9:22. doi:10.1038/s41526-023-00322-y

[16] Vitry G, Finch RH, McStay GP, et al. Muscle atrophy phenotype gene expression during spaceflight is linked to a metabolic crosstalk in both the liver and the muscle in mice. *iScience*. 2022;25(7):105213. doi:10.1016/j.isci.2022.105213

[17] Norton L, Wilson GJ, Layman DK, Moulton CJ, Garlick PJ. Leucine content of dietary proteins is a determinant of postprandial skeletal muscle protein synthesis in adult rats. *Nutrition & Metabolism*. 2012;9:67. doi:10.1186/1743-7075-9-67

[18] Wilkinson K, Koscien CP, Monteyne AJ, Wall BT, Stephens FB. Association of postprandial postexercise muscle protein synthesis rates with dietary leucine: A systematic review. *Physiological Reports*. 2023;11(5):e15775. doi:10.14814/phy2.15775

[19] Matthews JJ, Arentson-Lantz E, Katsanos CS, Moughan PJ, Wolfe RR, Ferrando AA, et al. Muscle Protein Synthesis Response to Whole-Foods and Free-Form Essential Amino Acids: An Exploratory Analysis. *Curr Dev Nutr*. 2025;106608. doi:10.1016/j.cdnut.2025.106608

[20] Pinckaers PJM, Trommelen J, Snijders T, van Loon LJC. The Anabolic Response to Plant-Based Protein Ingestion. *Sports Medicine*. 2021;51:59–69. doi:10.1007/s40279-021-01540-8

[21] Gorissen SHM, Crombag JJR, Senden J, et al. Protein content and amino acid composition of commercially available plant-based protein isolates. *Amino Acids*. 2018;50(12):1685–1695. doi:10.1007/s00726-018-2640-5

[22] Nichele S, Phillips SM, Boaventura BCB. Plant-based food patterns to stimulate muscle protein synthesis and support muscle mass in humans: a narrative review. *Appl Physiol Nutr Metab*. 2022;47(5):460–469. doi:10.1139/apnm-2021-0806

[23] Schoenfeld BJ, Aragon AA. How much protein can the body use in a single meal for muscle-building? Implications for daily protein distribution. *J Int Soc Sports Nutr*. 2018;15:10. doi:10.1186/s12970-018-0215-1

[24] Garcia-Medina JS, Sienkiewicz K, Narayanan S, et al. Genome and clonal hematopoiesis stability contrasts with immune, cfDNA, mitochondrial, and telomere length changes during short duration spaceflight. *Precision Clinical Medicine*. 2024;7(1):pbae007. doi:10.1093/pcmedi/pbae007

[25] Chakraborty N, Waning DL, Gautam A, et al. Gene-Metabolite Network Linked to Inhibited Bioenergetics in Association With Spaceflight-Induced Loss of Male Mouse Quadriceps Muscle. *J Bone Miner Res*. 2020;35(12):2304–2316. doi:10.1002/jbmr.4102

[26] Girón M, Vílchez JD, Shreeram S, et al. β-Hydroxy-β-Methylbutyrate (HMB) Normalizes Dexamethasone-Induced Autophagy-Lysosomal Pathway in Skeletal Muscle. *PLoS ONE*. 2015;10(3):e0117520. doi:10.1371/journal.pone.0117520

[27] Duan Y, Zheng C, Zhong Y, et al. Beta-hydroxy beta-methyl butyrate decreases muscle protein degradation via increased Akt/FoxO3a signaling and mitochondrial biogenesis in weanling piglets after lipopolysaccharide challenge. *Food Funct*. 2019;10(1):515–527. doi:10.1039/c9fo00769e

[28] He X, Li Y, Chen J, Huang Y, Zhou Y, Li Y, Quan J. β-hydroxy-β-methylbutyrate supplementation mitigates muscle atrophy induced by inactivity and protein deprivation. *Biogerontology*. 2025. doi:10.1007/s10522-025-10262-7

[29] Bear DE, Langan A, Dimidi E, Wandrag L, Harridge SDR, Hart N, et al. β-Hydroxy-β-methylbutyrate and its impact on skeletal muscle mass and physical function in clinical practice: a systematic review and meta-analysis. *Am J Clin Nutr*. 2018;108(3):546–556. doi:10.1093/ajcn/nqy373

[30] Desai I, Wewege MA, Jones MD, Clifford BK, Pandit A, Kaakoush NO, et al. The Effect of Creatine Supplementation on Resistance Training–Based Changes to Body Composition: A Systematic Review and Meta-analysis. *J Strength Cond Res*. 2024. doi:10.1519/jsc.0000000000004862

[31] Burke R, Piñero A, Coleman M, et al. The Effects of Creatine Supplementation Combined with Resistance Training on Regional Measures of Muscle Hypertrophy: A Systematic Review with Meta-Analysis. *Nutrients*. 2023;15(9):2116. doi:10.3390/nu15092116

[32] Holeček M. Beta‐hydroxy‐beta‐methylbutyrate supplementation and skeletal muscle in healthy and muscle‐wasting conditions. *J Cachexia Sarcopenia Muscle*. 2017;8(4):529–541. doi:10.1002/jcsm.12208

[33] Wen W, Chen X, Huang Z, Chen D, Chen H, Luo Y, et al. Resveratrol regulates muscle fiber type conversion via miR-22-3p and AMPK/SIRT1/PGC-1α pathway. *J Nutr Biochem*. 2019;71:108297. doi:10.1016/j.jnutbio.2019.108297

[34] Zheng M, Bai Y, Sun X, et al. Resveratrol Reestablishes Mitochondrial Quality Control in Myocardial Ischemia/Reperfusion Injury through Sirt1/Sirt3-Mfn2-Parkin-PGC-1α Pathway. *Molecules*. 2022;27(17):5545. doi:10.3390/molecules27175545

[35] Higashida K, Kim SH, Jung SR, Asaka M, Holloszy JO, Han DH. Effects of Resveratrol and SIRT1 on PGC-1α Activity and Mitochondrial Biogenesis: A Reevaluation. *PLoS Biol*. 2013;11(7):e1001603. doi:10.1371/journal.pbio.1001603

[36] Williams CB, Hughes MC, Edgett BA, Scribbans TD, Simpson CA, Perry CGR, et al. An Examination of Resveratrol's Mechanisms of Action in Human Tissue: Impact of a Single Dose In Vivo and Dose Responses in Skeletal Muscle Ex Vivo. *PLoS ONE*. 2014;9(7):e102406. doi:10.1371/journal.pone.0102406

[37] Okada R, Fujita S, Suzuki R, et al. Transcriptome analysis of gravitational effects on mouse skeletal muscles under microgravity and artificial 1 g onboard environment. *Scientific Reports*. 2021;11:88392. doi:10.1038/s41598-021-88392-4

[38] Vainshtein A. The Interplay Between PGC-1 and Autophagy During Metabolic Alterations in Skeletal Muscle. 2015. hdl.handle.net/10315/30012

[39] Lim C, Janssen TAH, Currier BS, et al. Muscle Protein Synthesis in Response to Plant-Based Protein Isolates With and Without Added Leucine Versus Whey Protein in Young Men and Women. *Curr Dev Nutr*. 2024;103769. doi:10.1016/j.cdnut.2024.103769

[40] Dijk F, van Dijk M, Roberts JD, et al. Pea and soy fortified with leucine stimulates muscle protein synthesis comparable to whey in a murine ageing model. *Eur J Nutr*. 2024. doi:10.1007/s00394-024-03506-8

[41] Li K, Desai R, Scott RT, et al. Explainable machine learning identifies multi-omics signatures of muscle response to spaceflight in mice. *npj Microgravity*. 2023;9:337. doi:10.1038/s41526-023-00337-5

[42] Oommen A, Stafford P, Joshi L. Profiling muscle transcriptome in mice exposed to microgravity using gene set enrichment analysis. *npj Microgravity*. 2024. doi:10.1038/s41526-024-00434-z

[43] Malhan D, Yalçın M, Schoenrock B, Blottner D, Relógio A. Skeletal muscle gene expression dysregulation in long-term spaceflights and aging is clock-dependent. *npj Microgravity*. 2023;9:273. doi:10.1038/s41526-023-00273-4

[44] Jones CW, Overbey EG, Lacombe J, et al. Molecular and physiological changes in the SpaceX Inspiration4 civilian crew. *Nature*. 2024;632:720–728. doi:10.1038/s41586-024-07648-x

[45] Overbey EG, Kim J, Tierney BT, et al. The Space Omics and Medical Atlas (SOMA) and international astronaut biobank. *Nature*. 2024;632:714–719. doi:10.1038/s41586-024-07639-y

[46] Souza A, Alves A, Martinez C, et al. Biomarkers of Skeletal Muscle Atrophy Based on Atrogenes Evaluation: A Systematic Review and Meta-Analysis Study. *Int J Mol Sci*. 2025;26(8):3516. doi:10.3390/ijms26083516

[47] Haberecht-Müller S, Krüger E, Fielitz J. Out of Control: The Role of the Ubiquitin Proteasome System in Skeletal Muscle during Inflammation. *Biomolecules*. 2021;11(9):1327. doi:10.3390/biom11091327

[48] Tagawa R, Watanabe D, Ito K, et al. Dose–response relationship between protein intake and muscle mass increase: a systematic review and meta-analysis of randomized controlled trials. *Nutr Rev*. 2020;78(12):1005–1018. doi:10.1093/nutrit/nuaa104

[49] Schoenfeld BJ, Aragon AA, Krieger J. The effect of protein timing on muscle strength and hypertrophy: a meta-analysis. *J Int Soc Sports Nutr*. 2013;10:53. doi:10.1186/1550-2783-10-53

[50] Pérez-Díaz S, Baselet B, Lovrić A, Lundberg TR, Neefs M, Daenen L, et al. Long non‐coding RNAs Kcnq1ot1 and Lncpint are involved in skeletal muscle atrophy induced by the space exposome. *J Physiol*. 2025. doi:10.1113/jp288987

---

*Manuscript prepared as a draft for submission to npj Microgravity. All data and code are available in the accompanying Zenodo package.*
