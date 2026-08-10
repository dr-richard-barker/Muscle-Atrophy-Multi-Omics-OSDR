# Plan: Supplement Validation + Non-Vegan Protein Benchmark

## Summary

Two tasks extending the completed muscle atrophy pipeline:

1. **Computational validation of HMB, CoQ10, and resveratrol** as counter-catabolic candidates, using three complementary in-silico approaches: (a) LINCS L1000 perturbation signature reversal analysis against the consensus catabolic signature, (b) molecular docking of each compound against its key protein targets using AutoDock Vina via HPC, and (c) literature evidence synthesis from published in vitro atrophy studies. A detailed wet-lab protocol for future experimental validation is also produced.

2. **Non-vegan protein benchmark**: Extend the amino acid database with 16+ animal protein sources, compute the same MPS scoring framework, and produce a head-to-head comparison table and figure benchmarking vegan vs. animal protein sources for leucine density, protein quality, caloric efficiency, and leucine threshold feasibility.

## Task 1: Computational Validation of HMB, CoQ10, Resveratrol

### 1A. LINCS L1000 Signature Reversal Analysis

**What**: For each of the 3 target compounds (plus creatine and nicotinamide riboside as positive controls), extract all available LINCS L1000 perturbation signatures from the local GMT file and compute a formal reversal score against the consensus catabolic signature.

**Data sources**:
- Local GMT: `/mnt/datalake/LINCS1000/RNAseq_transcriptomics_genesets/single_drug_perturbations-v1.0.gmt` (542 drug signatures)
- Consensus signature: `/workspace/osdr_pipeline/results/signatures/consensus_catabolic_signature.csv` (77 up, 64 down)
- Broad Drug Repurposing Hub: SMILES and metadata for compound identification

**Method**: For each compound's up/down GMT signature pair:
- Compute hypergeometric overlap p-value between drug-up genes and atrophy-down genes (reversal direction)
- Compute hypergeometric overlap p-value between drug-dn genes and atrophy-up genes (reversal direction)
- Compute a composite reversal score: −log10(p) weighted by direction concordance
- Also compute a "mimicry" score (drug-up vs atrophy-up) as a negative control
- Aggregate across multiple signatures per compound (resveratrol has 13 signature pairs, creatine 1, NR 4)
- Compare against the full compound ranking from the original LINCS screen (559 compounds) to contextualize

**Available signatures**:
- Resveratrol: 13 up/dn pairs (26 GMT entries)
- Creatine: 1 up/dn pair
- Nicotinamide riboside: 4 up/dn pairs
- HMB: 0 entries in LINCS GMT (not in L1000) — will note as limitation
- CoQ10: 0 direct entries (coenzyme-Q10 exists in Broad hub but no LINCS perturbation signatures) — will note as limitation

**Output**: `supplement_validation_lincs_reversal.csv` — per-compound, per-signature reversal scores with hypergeometric p-values, overlap gene lists, and aggregate scores.

### 1B. Molecular Docking Against Key Protein Targets

**What**: Dock each compound (HMB, CoQ10, resveratrol, plus creatine as control) against its primary protein targets using AutoDock Vina via the HPC large-scale-virtual-screening tool.

**Compound SMILES** (from Broad Drug Repurposing Hub):
- Resveratrol: `Oc1ccc(\C=C\c2cc(O)cc(O)c2)cc1` (PubChem CID 445154)
- Creatine: `CN(CC(O)=O)C(N)=N` (PubChem CID 586)
- HMB (beta-hydroxy-beta-methylbutyrate): `CC(C)(O)CC(O)=O` (PubChem CID 129010087)
- CoQ10 (coenzyme-Q10): long SMILES from Broad hub (PubChem CID 5281915)

**Protein targets and PDB structures**:
| Compound | Target | PDB ID | Rationale |
|----------|--------|--------|-----------|
| Resveratrol | SIRT1 | 4I5I | Co-crystal with resveratrol; validates docking protocol |
| Resveratrol | FOXO3 | 2CO7 | DNA-binding domain; resveratrol inhibits FOXO3 transcriptional activity |
| HMB | MTOR (kinase domain) | 4JSN | HMB activates mTOR; docking predicts binding site |
| CoQ10 | Complex I (MT-ND1) | AlphaFold predicted | No PDB for mitochondrial-encoded; use AlphaFold structure |
| Creatine | CKM (creatine kinase) | 1CRK | Native substrate-enzyme structure; validates protocol |

**Method**:
1. Retrieve PDB structures (download from RCSB or predict via AlphaFold HPC for FBXO32/MT-ND1)
2. Prepare receptor: remove water, add hydrogens, compute charges
3. Prepare ligands from SMILES using RDKit (3D conformer generation)
4. Dock with AutoDock Vina (exhaustiveness=32, box centered on known active site or whole protein for AlphaFold models)
5. Report binding affinity (kcal/mol), RMSD, and key interacting residues
6. For SIRT1-resveratrol: validate against known co-crystal pose (4I5I) as positive control

**HPC execution**: Submit as `hpc_run_tool` with AutoDock Vina. Each docking job is ~5-15 min. Total: ~5 docking jobs (4 compounds × 1-2 targets each, plus 1 validation redock).

**Output**: `supplement_docking_results.csv` — per-compound, per-target binding affinities, RMSD, interacting residues. Docking pose PDBQT files saved to results.

### 1C. Literature Evidence Synthesis

**What**: Systematic literature search for published in vitro evidence of HMB, CoQ10, and resveratrol in muscle atrophy models (C2C12 myotubes, primary myotubes, dexamethasone-induced atrophy, serum starvation, clinostat/RPM microgravity simulation).

**Method**: Use `LiteratureSearch` to find papers for each compound × atrophy model combination. Extract: cell model, atrophy induction method, compound dose, key readouts (myotube diameter, atrogene expression, MPS signaling), and effect direction. Produce a structured evidence table.

**Search queries**:
- "HMB C2C12 myotube atrophy dexamethasone"
- "CoQ10 coenzyme Q10 muscle atrophy in vitro myotube"
- "Resveratrol C2C12 myotube atrophy dexamethasone microgravity"
- "Resveratrol muscle atrophy clinostat simulated microgravity"
- "HMB simulated microgravity muscle"

**Output**: `supplement_validation_literature_evidence.csv` — per-study evidence with model, dose, readout, effect size, and citation.

### 1D. Wet-Lab Protocol Document

**What**: A detailed experimental protocol for future wet-lab validation, specifying cell line, atrophy induction, microgravity simulation, compound treatment, readouts, controls, and statistical analysis plan.

**Content**:
- Cell model: C2C12 mouse myoblasts → differentiated myotubes (7 days)
- Atrophy induction: dexamethasone (10 μM, 48h) or serum starvation (0.5% FBS, 48h)
- Microgravity simulation: 3D clinostat or random positioning machine (RPM), 48-72h
- Compound treatment: HMB (5 mM), CoQ10 (10 μM), resveratrol (50 μM) — doses from literature
- Controls: vehicle (DMSO), atrophy-only, compound-only, positive control (IGF-1, 100 ng/mL)
- Readouts: myotube diameter (immunofluorescence microscopy), atrogene qPCR (FBXO32, TRIM63, FOXO3), MPS signaling (p-S6K1, p-4E-BP1 Western blot), mitochondrial membrane potential (JC-1)
- Experimental design: 2×2×2 factorial (atrophy × compound × microgravity), n=6 per condition, 3 biological replicates
- Statistics: two-way ANOVA with Tukey post-hoc, p<0.05, effect sizes reported

**Output**: `in_vitro_validation_protocol.md` — standalone protocol document.

### 1E. Validation Summary Figure

**What**: A multi-panel figure summarizing the computational validation results.

**Panels**:
- Panel A: LINCS reversal scores for the 5 compounds (bar chart with hypergeometric p-values)
- Panel B: Docking binding affinities (heatmap: compound × target)
- Panel C: Literature evidence summary (evidence level × compound matrix)
- Panel D: Integrated validation score (composite of LINCS + docking + literature)

**Output**: `fig9_supplement_validation.svg` + `.png` (300 DPI)

## Task 2: Non-Vegan Protein Benchmark

### 2A. Animal Protein Amino Acid Database

**What**: Build a comprehensive amino acid database for 16+ animal protein sources using published USDA SR Legacy values (same source as the vegan database for consistency).

**Sources to include** (comprehensive panel):
1. Whey protein isolate
2. Casein protein
3. Egg white (dried/powdered)
4. Egg, whole (raw)
5. Chicken breast (raw)
6. Lean beef (top sirloin, raw)
7. Salmon (Atlantic, raw)
8. Tuna (light, canned in water)
9. Greek yogurt (nonfat)
10. Milk (whole)
11. Cottage cheese (lowfat 2%)
12. Collagen peptides (hydrolyzed)
13. Beef protein isolate
14. Pork tenderloin (raw)
15. Turkey breast (raw)
16. Cod (Atlantic, raw)
17. Shrimp (raw)

**Data**: 18 amino acids per 100g, total protein, energy (kcal), from published SR Legacy values. Same column structure as the vegan database.

### 2B. MPS Scoring and Comparison

**What**: Apply the identical MPS scoring framework (40% leucine density + 25% protein quality + 20% caloric efficiency + 15% lysine density) to the animal protein sources, then produce a unified ranking combining vegan + animal sources.

**Method**:
- Compute all derived metrics (BCAA, leucine per 100g protein, grams for leucine threshold, kcal for leucine threshold) for animal sources
- Z-score normalize across the COMBINED dataset (vegan + animal) so scores are directly comparable
- Compute MPS scores on the combined scale
- Rank all sources together

### 2C. Comparative Benchmark Outputs

**What**: Tables and figures comparing vegan vs. animal protein sources.

**Outputs**:
1. `animal_amino_acid_database.csv` — 17 animal sources with 18 amino acids
2. `combined_protein_mps_ranking.csv` — unified ranking of all 58 sources (41 vegan + 17 animal)
3. `vegan_vs_animal_benchmark.csv` — head-to-head comparison table (top 10 vegan vs. top 10 animal by MPS score, with per-metric comparison)
4. `fig10_vegan_vs_animal_comparison.svg` + `.png` — multi-panel figure:
   - Panel A: Side-by-side bar chart of top 10 vegan vs. top 10 animal by MPS score
   - Panel B: Leucine per 100g scatter plot (vegan vs. animal, colored by category)
   - Panel C: Protein quality vs. caloric efficiency scatter (all 58 sources)
   - Panel D: Grams needed for 2.5g leucine threshold (horizontal bar chart, vegan vs. animal)

## Compute & Resource Estimate

| Component | Data Volume | RAM | Runtime | Target |
|-----------|-------------|-----|---------|--------|
| LINCS reversal analysis (Python) | <1 MB | 4 GB | ~5 min | worker-0 |
| Molecular docking (5 jobs via HPC) | <100 MB | 8 GB | ~15 min/job, parallel | HPC (AutoDock Vina) |
| Literature search (5 queries) | <1 MB | 4 GB | ~5 min | worker-0 |
| Animal protein database + scoring | <1 MB | 4 GB | ~10 min | worker-0 |
| Figures (2 figures, SVG+PNG) | <10 MB | 4 GB | ~10 min | worker-0 |
| **Total** | **<120 MB** | **8 GB** | **~30 min** (docking parallel) | **worker-0 + HPC** |

No additional machines needed. Docking jobs submit to HPC in parallel. All other work on worker-0.

## Key Assumptions

1. LINCS L1000 GMT signatures for resveratrol, creatine, and nicotinamide riboside are valid proxies for in vitro transcriptional effects (caveat: non-muscle cell lines)
2. HMB and CoQ10 lack LINCS perturbation signatures — their validation relies on docking and literature evidence only
3. PDB structures for SIRT1 (4I5I), FOXO3 (2CO7), MTOR (4JSN), CKM (1CRK) are available from RCSB
4. AlphaFold HPC can predict structures for MT-ND1 and FBXO32 if no PDB exists
5. USDA SR Legacy amino acid values for animal proteins are available from the same published sources used for the vegan database
6. The MPS scoring framework (designed for vegan sources) is applicable to animal sources without modification — the z-score normalization across the combined dataset ensures fair comparison
7. Docking binding affinities are computational predictions, not experimental measurements — they support plausibility but do not prove efficacy

## Testing & Acceptance Criteria

- LINCS reversal analysis produces per-compound scores with hypergeometric p-values for all available signatures
- Docking produces binding affinity (kcal/mol) for each compound-target pair, with SIRT1-resveratrol redock RMSD < 2.0 Å as protocol validation
- Literature search returns ≥5 relevant in vitro studies per compound
- Animal protein database contains ≥17 sources with complete amino acid profiles (all 18 AAs)
- Combined MPS ranking includes all 58 sources with directly comparable z-scores
- Both figures pass media output check (non-blank, readable, properly labeled)
- All outputs saved to `/mnt/results/` and added to the Zenodo package
