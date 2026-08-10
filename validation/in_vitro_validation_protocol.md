# In Vitro Validation Protocol: HMB, CoQ10, and Resveratrol in a Simulated Microgravity Muscle Atrophy Model

## Objective

To experimentally validate the counter-catabolic potential of beta-hydroxy-beta-methylbutyrate (HMB), coenzyme Q10 (CoQ10), and resveratrol in an in vitro muscle atrophy model combining glucocorticoid-induced catabolism with simulated microgravity, using C2C12 differentiated myotubes as the primary cell model.

## Rationale

Computational analysis (LINCS L1000 signature reversal, molecular docking, and literature synthesis) identified these three compounds as high-priority candidates for counteracting spaceflight-induced muscle atrophy. Resveratrol shows the strongest literature evidence (10 studies, 5 in microgravity/disuse models) and the best predicted binding to SIRT1 (-8.6 kcal/mol). HMB has robust evidence in dexamethasone-induced atrophy models (4 strong studies) and has been tested in human bed-rest (microgravity analog). CoQ10 supports mitochondrial function but has less direct atrophy evidence. This protocol provides the experimental framework to validate these computational predictions.

## Cell Model

### Cell Line
- **C2C12 mouse skeletal myoblasts** (ATCC CRL-1772)
- Passage number: 5-15 (avoid >20 passages)
- Mycoplasma testing: confirm negative before experiments

### Culture Conditions
- Growth medium: DMEM (high glucose, 4.5 g/L) + 10% FBS + 1% penicillin-streptomycin
- Differentiation medium: DMEM + 2% horse serum + 1% penicillin-streptomycin
- Incubation: 37C, 5% CO2, humidified atmosphere

### Differentiation Protocol
1. Seed C2C12 myoblasts at 2 x 10^5 cells/well in 6-well plates (or 35 mm dishes)
2. Culture in growth medium until ~90% confluence (typically 48h)
3. Switch to differentiation medium to induce fusion
4. Maintain in differentiation medium for 7 days, changing medium every 48h
5. Confirm myotube formation: >70% of nuclei within multi-nucleated myotubes (>=3 nuclei per tube)
6. Use mature myotubes (day 7 post-differentiation) for all experiments

## Atrophy Induction

### Method 1: Dexamethasone (Primary)
- **Agent**: Dexamethasone (Sigma D4902)
- **Dose**: 10 uM (dissolved in DMSO, final DMSO concentration <=0.1%)
- **Duration**: 48 hours
- **Rationale**: Well-established glucocorticoid-induced atrophy model; activates FOXO transcription factors, upregulates atrogin-1/MAFbx and MuRF1, inhibits mTORC1 signaling. Consistent with the consensus catabolic signature from our meta-analysis (FOXO3 log2FC=+0.221, FBXO32 log2FC=+0.303, TRIM63 log2FC=-0.196).

### Method 2: Serum Starvation (Alternative)
- **Condition**: 0.5% FBS in DMEM (vs. 2% horse serum in differentiation medium)
- **Duration**: 48 hours
- **Rationale**: Models nutrient deprivation component of microgravity-induced atrophy

## Simulated Microgravity

### Apparatus
- **3D Clinostat** (random positioning machine) OR **Random Positioning Machine (RPM)**
- Clinostat rotation: 2 rotations per minute (rpm) around two orthogonal axes
- Gravity vector averaging: effectively <10^-3 g at the cell center

### Microgravity Duration
- 48-72 hours concurrent with atrophy induction and compound treatment
- Control: identical plates maintained in normal gravity (1g) incubator

### Plate Compatibility
- Use 35 mm Petri dishes (not standard well plates) to ensure uniform fluid distribution during rotation
- Fill dishes completely with medium (no air gap) to prevent fluid shear artifacts
- Seal with parafilm to prevent leakage

## Compound Treatment

### Compounds and Doses

| Compound | Stock Solution | Working Concentration | Solvent | Rationale |
|----------|---------------|----------------------|---------|-----------|
| HMB (calcium salt) | 100 mM in water | 5 mM | Water | Based on Aversa et al. 2012 [1], Baptista et al. 2013 [2]; activates p70S6K/mTOR pathway |
| CoQ10 | 10 mM in DMSO | 10 uM | DMSO (<=0.1%) | Based on Wagner et al. 2012 [3]; mitochondrial support concentration |
| Resveratrol | 50 mM in DMSO | 50 uM | DMSO (<=0.1%) | Based on Alamdari et al. 2011 [4], Wang et al. 2014 [5]; SIRT1 activation dose |

### Treatment Groups

| Group | Atrophy Induction | Microgravity | Compound | n |
|-------|------------------|-------------|----------|---|
| 1. Vehicle control (1g) | - | - | Vehicle | 6 |
| 2. Vehicle control (ug) | - | + | Vehicle | 6 |
| 3. Atrophy only (1g) | Dexamethasone | - | Vehicle | 6 |
| 4. Atrophy only (ug) | Dexamethasone | + | Vehicle | 6 |
| 5. HMB only (1g) | - | - | HMB 5 mM | 6 |
| 6. HMB only (ug) | - | + | HMB 5 mM | 6 |
| 7. CoQ10 only (1g) | - | - | CoQ10 10 uM | 6 |
| 8. CoQ10 only (ug) | - | + | CoQ10 10 uM | 6 |
| 9. Resveratrol only (1g) | - | - | RSV 50 uM | 6 |
| 10. Resveratrol only (ug) | - | + | RSV 50 uM | 6 |
| 11. Atrophy + HMB (1g) | Dexamethasone | - | HMB 5 mM | 6 |
| 12. Atrophy + HMB (ug) | Dexamethasone | + | HMB 5 mM | 6 |
| 13. Atrophy + CoQ10 (1g) | Dexamethasone | - | CoQ10 10 uM | 6 |
| 14. Atrophy + CoQ10 (ug) | Dexamethasone | + | CoQ10 10 uM | 6 |
| 15. Atrophy + RSV (1g) | Dexamethasone | - | RSV 50 uM | 6 |
| 16. Atrophy + RSV (ug) | Dexamethasone | + | RSV 50 uM | 6 |
| 17. Positive control IGF-1 (1g) | Dexamethasone | - | IGF-1 100 ng/mL | 6 |
| 18. Positive control IGF-1 (ug) | Dexamethasone | + | IGF-1 100 ng/mL | 6 |

**Total: 18 conditions x 6 replicates = 108 dishes per biological replicate**
**Biological replicates: 3 independent experiments**
**Total dishes: 324**

### Treatment Timeline
1. Day 0: Seed C2C12 myoblasts
2. Day 2: Switch to differentiation medium (confluence reached)
3. Day 9: Mature myotubes formed (7 days post-differentiation)
4. Day 9: Begin treatment (add compounds + dexamethasone simultaneously)
5. Day 9: Place ug groups on clinostat/RPM
6. Day 11: End treatment (48h), harvest for readouts

## Readouts

### Primary Readouts

#### 1. Myotube Diameter (Morphology)
- **Method**: Immunofluorescence microscopy
- **Staining**: Anti-myosin heavy chain (MF20 antibody, DSHB) + DAPI
- **Imaging**: 5 random fields per dish, 20x objective
- **Measurement**: Measure diameter of >=50 myotubes per condition using ImageJ
- **Primary endpoint**: Mean myotube diameter (um)
- **Expected atrophy effect**: 20-30% reduction in diameter (dexamethasone) [1,4]

#### 2. Atrogene Expression (qPCR)
- **Targets**: Fbxo32 (Atrogin-1/MAFbx), Trim63 (MuRF1), Foxo3
- **Reference genes**: Gapdh, Rpl13a (geometric mean)
- **Method**: RT-qPCR (SYBR Green or TaqMan)
- **RNA extraction**: TRIzol or column-based (RNeasy Mini Kit)
- **cDNA**: SuperScript IV, oligo(dT) primers
- **qPCR**: 3 technical replicates per sample
- **Analysis**: delta-delta Ct method, fold change relative to vehicle control (1g)
- **Expected atrophy effect**: 3-10 fold upregulation of Fbxo32 and Trim63 [1,4,5]

### Secondary Readouts

#### 3. MPS Signaling (Western Blot)
- **Targets**: 
  - p-S6K1 (Thr389) / total S6K1 (mTORC1 activity)
  - p-4E-BP1 (Thr37/46) / total 4E-BP1 (mTORC1 activity)
  - p-AMPK (Thr172) / total AMPK (energy sensing)
  - SIRT1 protein level
  - p-FoxO1 (Ser256) / total FoxO1 (FOXO activity)
- **Method**: Western blot with chemiluminescence
- **Loading control**: beta-actin or vinculin
- **Quantification**: Densitometry (ImageJ), phospho/total ratio normalized to loading control
- **Expected HMB effect**: Increased p-S6K1/p-4E-BP1 (mTORC1 activation) [1,6]
- **Expected RSV effect**: Increased SIRT1, increased p-AMPK, decreased p-FoxO1 [4,5]

#### 4. Mitochondrial Membrane Potential (JC-1)
- **Method**: JC-1 fluorescent dye (5 ug/mL, 30 min, 37C)
- **Imaging**: Fluorescence microscopy or plate reader
- **Readout**: Red/green fluorescence ratio (healthy mitochondria = high red/green)
- **Expected CoQ10 effect**: Preserved red/green ratio under atrophy conditions [3]
- **Expected atrophy effect**: Decreased red/green ratio (depolarization) [7]

#### 5. Protein Synthesis Rate (Optional)
- **Method**: Surface sensing of translation (SUnSET) method with puromycin
- **Puromycin**: 1 ug/mL for 30 min before harvest
- **Detection**: Anti-puromycin Western blot (clone 12D10)
- **Readout**: Global protein synthesis rate

## Controls

### Negative Controls
- Vehicle-only (DMSO <=0.1% or water) at 1g and ug
- Untreated myotubes at 1g and ug

### Positive Controls
- **IGF-1 (100 ng/mL)**: Known anti-atrophic factor; activates PI3K/Akt/mTOR pathway
- Expected to partially prevent dexamethasone-induced atrophy

### Process Controls
- DMSO-only group (solvent control for CoQ10 and resveratrol)
- Water-only group (solvent control for HMB)

## Experimental Design

### Statistical Analysis Plan

#### Primary Analysis
- **Design**: 3-way factorial (atrophy x compound x gravity) with 6 replicates per condition
- **Analysis**: Three-way ANOVA with Tukey's HSD post-hoc test
- **Primary comparison**: Atrophy+compound (ug) vs. Atrophy-only (ug)
- **Significance threshold**: p < 0.05 (two-tailed)
- **Effect sizes**: Report Cohen's d for all pairwise comparisons
- **Multiple testing**: Tukey HSD controls family-wise error rate

#### Secondary Analysis
- **Dose-response**: If single dose shows effect, follow with 3-dose range finding
- **Interaction effects**: Test atrophy x gravity interaction (synergy between catabolism and microgravity)
- **Correlation**: Correlate myotube diameter with atrogene expression and signaling markers

#### Power Analysis
- Based on published C2C12 atrophy studies [1,4,5]:
  - Expected effect size: 20-30% diameter preservation
  - SD: ~10% of mean diameter
  - Effect size (Cohen's d): ~2.0 (large)
  - With n=6 per group: power >0.95 at alpha=0.05
  - With n=4 per group: power >0.80 at alpha=0.05
- n=6 provides adequate power with margin for technical failures

### Randomization and Blinding
- Randomize dish placement in incubator and on clinostat
- Blind the microscopist during image acquisition and diameter measurement
- Blind qPCR analysis (code samples before cDNA synthesis)

## Expected Outcomes

### Based on Computational Predictions and Literature

| Compound | Predicted Effect on Atrophy | Key Mechanism | Evidence Strength |
|----------|---------------------------|---------------|-------------------|
| Resveratrol | Strong protection | SIRT1 activation -> FoxO1/3 inhibition -> reduced atrogene expression | Strong (10 studies, 5 in microgravity models) |
| HMB | Moderate protection | mTORC1/p70S6K activation -> increased protein synthesis; PLD2 pathway | Strong (8 studies, 1 in bed rest) |
| CoQ10 | Modest protection | Mitochondrial membrane potential preservation; reduced oxidative stress | Moderate (6 studies, 1 in spaceflight) |

### Docking Predictions (AutoDock Vina + Boltz-2)
- Resveratrol-SIRT1: -8.6 kcal/mol (strong predicted binding, consistent with known SIRT1 activator)
- HMB-MTOR: -4.5 kcal/mol (moderate, consistent with indirect mTOR activation)
- CoQ10-SIRT1: Boltz-2 confidence 0.806, ligand_iptm 0.866 (high confidence ligand placement)
- CoQ10: Too large for AutoDock Vina (50+ rotatable bonds); Boltz-2 structure prediction used instead

## Timeline

| Day | Activity |
|-----|----------|
| -7 | Thaw C2C12 cells, expand in growth medium |
| 0 | Seed cells into 35 mm dishes (324 dishes) |
| 2 | Switch to differentiation medium |
| 4 | Change differentiation medium |
| 6 | Change differentiation medium |
| 8 | Change differentiation medium; confirm differentiation |
| 9 | Begin treatment: add compounds + dexamethasone; place ug groups on clinostat |
| 11 | End treatment (48h); harvest for readouts |
| 11-13 | Immunofluorescence imaging, RNA extraction, Western blot |
| 14-21 | qPCR, data analysis |
| 22-28 | Repeat for biological replicates 2 and 3 |
| 30 | Final data analysis and reporting |

## Materials List

### Cell Culture
- C2C12 cells (ATCC CRL-1772)
- DMEM high glucose (Gibco 11965)
- FBS (Gibco 26140)
- Horse serum (Gibco 16050)
- Pen/Strep (Gibco 15140)
- 35 mm Petri dishes (Corning 430588)
- 0.25% Trypsin-EDTA (Gibco 25200)

### Compounds
- HMB calcium monohydrate (Sigma 12637)
- Coenzyme Q10 (Sigma C9538)
- Resveratrol (Sigma R5010)
- Dexamethasone (Sigma D4902)
- IGF-1 recombinant human (PeproTech AF-100-11)
- DMSO (Sigma D8418)

### Readout Reagents
- MF20 anti-myosin heavy chain antibody (DSHB)
- Secondary antibody: Anti-mouse IgG Alexa Fluor 488 (Invitrogen A11001)
- DAPI (Invitrogen D1306)
- TRIzol (Invitrogen 15596) or RNeasy Mini Kit (Qiagen 74104)
- SuperScript IV (Invitrogen 18091)
- SYBR Green Master Mix (Applied Biosystems 4309155)
- Primers for Fbxo32, Trim63, Foxo3, Gapdh, Rpl13a (see primer sequences below)
- JC-1 dye (Invitrogen T3168)
- RIPA buffer + protease/phosphatase inhibitors
- Western blot antibodies:
  - p-S6K1 (Thr389) (Cell Signaling 9234)
  - total S6K1 (Cell Signaling 2708)
  - p-4E-BP1 (Thr37/46) (Cell Signaling 2855)
  - total 4E-BP1 (Cell Signaling 9644)
  - p-AMPK (Thr172) (Cell Signaling 2535)
  - total AMPK (Cell Signaling 5831)
  - SIRT1 (Cell Signaling 9475)
  - p-FoxO1 (Ser256) (Cell Signaling 9461)
  - total FoxO1 (Cell Signaling 2880)
  - beta-actin (Cell Signaling 4970)

### Primer Sequences (Mouse)
| Gene | Forward (5'->3') | Reverse (5'->3') |
|------|-----------------|------------------|
| Fbxo32 | TGAGTGGAGCGACTTGATGA | AATCCGTTCTTCTGCCTTCC |
| Trim63 | CCGAGGTGCAAGAGAGAAGA | GGTCCAGGAGGAGAAGGTG |
| Foxo3 | AGAAGTTACTGTTGGCGAGC | TTGAGCAGGGTAGAGGTGGT |
| Gapdh | AGGTCGGTGTGAACGGATTTG | TGTAGACCATGTAGTTGAGGTCA |
| Rpl13a | GGGCCTGGAGAGGTTTG | CTTGGCCTTTTCCTTCCGTT |

## Quality Control

### Inclusion/Exclusion Criteria
- Exclude dishes with <50% myotube coverage at harvest
- Exclude dishes with evidence of contamination
- Exclude qPCR samples with Ct > 35 for reference genes
- Exclude Western blot lanes with uneven loading (>20% variation in loading control)

### Validation Checkpoints
1. **Differentiation quality**: >70% nuclei in multi-nucleated myotubes (day 7)
2. **Atrophy induction**: Confirm >=20% diameter reduction in atrophy-only vs. vehicle (1g)
3. **Atrogene induction**: Confirm >=2-fold upregulation of Fbxo32 in atrophy-only vs. vehicle
4. **Positive control**: IGF-1 should partially prevent atrophy (>=10% diameter preservation)
5. **DMSO tolerance**: No significant difference between vehicle and DMSO-only groups

## Limitations and Considerations

1. **C2C12 vs. primary myotubes**: C2C12 cells are an immortalized line; primary mouse or human myotubes may show different responses. Consider confirming key findings in primary cells.

2. **Dexamethasone vs. microgravity**: Dexamethasone models glucocorticoid-induced catabolism, which is one component of microgravity atrophy but not the complete picture. The combined dexamethasone + clinostat model better approximates the spaceflight condition.

3. **Clinostat limitations**: 3D clinostats average the gravity vector to near-zero but do not replicate true microgravity. Results should be interpreted as "simulated microgravity" not "microgravity."

4. **Single dose**: This protocol tests a single dose per compound based on literature. Dose-response studies should follow if effects are observed.

5. **48-hour treatment**: This is a short-term model. Chronic exposure (7-14 days) may reveal different effects, particularly for CoQ10 (mitochondrial adaptation).

6. **Compound solubility**: CoQ10 is poorly soluble in aqueous media. DMSO stock (10 mM) diluted to 10 uM should be verified for precipitation. Consider liposomal CoQ10 formulations if solubility is problematic.

7. **Species translation**: C2C12 cells are mouse-derived. Human myotube studies (e.g., from primary human myoblasts or iPSC-derived myotubes) would strengthen translational relevance for astronaut nutrition planning.

## References

1. Aversa Z, et al. beta-Hydroxy-beta-methylbutyrate (HMB) prevents dexamethasone-induced myotube atrophy. Biochem Biophys Res Commun. 2012. doi:10.1016/j.bbrc.2012.06.029
2. Baptista IL, et al. Leucine and HMB differentially modulate proteasome system in skeletal muscle under different sarcopenic conditions. PLoS ONE. 2013. doi:10.1371/journal.pone.0076752
3. Wagner AE, et al. A combination of lipoic acid plus coenzyme Q10 induces PGC1alpha and improves stress response in C2C12 cells. Oxid Med Cell Longev. 2012. doi:10.1155/2012/835970
4. Alamdari N, et al. Resveratrol prevents dexamethasone-induced expression of the muscle atrophy-related ubiquitin ligases atrogin-1 and MuRF1 in cultured myotubes through a SIRT1-dependent mechanism. Biochem Biophys Res Commun. 2011. doi:10.1016/j.bbrc.2011.11.154
5. Wang DT, et al. Resveratrol prevents TNF-alpha-induced muscle atrophy via regulation of Akt/mTOR/FoxO1 signaling in C2C12 myotubes. Int Immunopharmacol. 2014. doi:10.1016/j.intimp.2014.02.002
6. Cohen-Or M, et al. HMB leads to phospholipase D2 activation and alters circadian rhythms in myotubes. Food Funct. 2024. doi:10.1039/d3fo04174c
7. Liu J, et al. Mitochondrial dysfunction launches dexamethasone-induced skeletal muscle atrophy via AMPK/FOXO3 signaling. Mol Pharm. 2016. doi:10.1021/acs.molpharmaceut.5b00516
8. Momken I, et al. Resveratrol prevents the wasting disorders of mechanical unloading by acting as a physical exercise mimetic in the rat. FASEB J. 2011. doi:10.1096/fj.10-177295
9. Mortreux M, et al. A moderate daily dose of resveratrol mitigates muscle deconditioning in a Martian gravity analog. Front Physiol. 2019. doi:10.3389/fphys.2019.00899
10. Standley RA, et al. Effects of HMB on skeletal muscle mitochondrial content and dynamics after 10 days of bed rest in older adults. J Appl Physiol. 2017. doi:10.1152/japplphysiol.00192.2017
