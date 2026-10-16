#!/usr/bin/env python3
"""
05_lincs_screening.py — LINCS L1000 connectivity mapping for counter-catabolic compounds.

Given the consensus catabolic signature (77 up / 64 down genes from meta-analysis),
query iLINCS API + local GMT fallback to find compounds that REVERSE the atrophy signature.

Strategy:
  1. Query iLINCS API (findConcordancesSC) with up-regulated atrophy genes + known atrophy markers
     -> Compounds with HIGH z-scores MIMIC atrophy (bad)
  2. Query iLINCS API with down-regulated genes + known anti-atrophy/MPS genes
     -> Compounds with HIGH z-scores PROMOTE muscle growth genes (good = reversers)
  3. Local GMT fallback: hypergeometric reversal scoring using LINCS L1000 drug perturbation GMT
  4. Composite reversal score: 50% promote_muscle + 30% anti_mimic + 20% GMT evidence

Inputs:
  /workspace/osdr_pipeline/results/signatures/consensus_catabolic_signature.csv
  /workspace/osdr_pipeline/results/signatures/ranked_gene_list.csv

Outputs:
  lincs_reverser_ranking_full.csv     — all ranked compounds
  lincs_reverser_top100.csv           — top 100 candidates
  ilincs_up_gene_query_raw.csv        — raw iLINCS results (up gene query)
  ilincs_down_gene_query_raw.csv      — raw iLINCS results (down gene query)
  lincs_gmt_reversal_raw.csv          — local GMT reversal results
  robustness_summary.json             — analysis summary

Usage:
  python 05_lincs_screening.py
"""
import argparse, json, os, re, sys, time
import numpy as np
import pandas as pd
import requests
from scipy.stats import hypergeom

# --- Constants ---
ILINCS_API = "https://www.ilincs.org/api"
GMT_DRUG = "/mnt/datalake/LINCS1000/RNAseq_transcriptomics_genesets/single_drug_perturbations-v1.0.gmt"
HUB_MOL = "/mnt/datalake/broad_drug_repurposing_hub/broad_repurposing_hub_molecule_with_smiles.parquet"
HUB_INFO = "/mnt/datalake/broad_drug_repurposing_hub/broad_repurposing_hub_phase_moa_target_info.parquet"

# Known atrophy/MPS genes to augment the query
KNOWN_ATROPHY = ['FOXO3', 'FBXO32', 'TRIM63', 'PPARGC1A', 'MTOR', 'AKT1', 'RPS6KB1',
                 'EIF4EBP1', 'ULK1', 'SQSTM1', 'BNIP3', 'GABARAPL1',
                 'HIF1A', 'NFKB1', 'TNF', 'IL6', 'MYOD1', 'MYOG', 'PAX7', 'MSTN',
                 'IGF1', 'PIK3CA', 'TSC2', 'RHEB', 'DDIT4']

KNOWN_ANTI_ATROPHY = ['MYOD1', 'MYOG', 'PAX7', 'IGF1', 'PIK3CA', 'AKT1', 'MTOR', 'RPS6KB1',
                      'RHEB', 'PPARGC1A', 'MYH1', 'MYH2', 'MYH4', 'ACTA1', 'MYL1', 'MYL2',
                      'TNNI1', 'TNNC1', 'CKM', 'CKB', 'GYS1', 'PHKA1', 'GPD1', 'GAPDH']


def read_signature(path):
    """Read consensus catabolic signature CSV."""
    df = pd.read_csv(path)
    up = df[df['signature_direction'] == 'up']['gene'].str.upper().tolist()
    down = df[df['signature_direction'] == 'down']['gene'].str.upper().tolist()
    return up, down


def parse_gmt(path):
    """Parse GMT file into dict of gene sets."""
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = {g.split(",")[0].strip().upper() for g in parts[2:] if g.strip()}
    return sets


def ilincs_query(genes, timeout=120):
    """Query iLINCS findConcordancesSC API with a gene list."""
    body = {
        "mode": "geneList",
        "metadata": True,
        "signatureProfile": {"genes": genes}
    }
    r = requests.post(f"{ILINCS_API}/ilincsR/findConcordancesSC",
                      json=body, timeout=timeout,
                      headers={'Accept': 'application/json', 'Content-Type': 'application/json'})
    r.raise_for_status()
    data = r.json()
    if 'sigScores' in data and data['sigScores']:
        return pd.DataFrame(data['sigScores'])
    return pd.DataFrame()


def gmt_reversal_scoring(drug_sets, up_set, down_set, universe):
    """Local GMT hypergeometric reversal scoring."""
    rows = []
    for name, genes in drug_sets.items():
        low = name.lower()
        is_up = ("-up" in low) or low.endswith(" up") or ("_up" in low)
        is_dn = ("-dn" in low) or ("down" in low) or ("_dn" in low) or ("-down" in low)
        if not (is_up or is_dn):
            continue
        if is_dn:
            query = up_set
            direction = "drug_down_reverses_query_up"
        else:
            query = down_set
            direction = "drug_up_reverses_query_down"
        k = len(genes & query)
        if k == 0:
            continue
        p = hypergeom.sf(k - 1, universe, len(query), len(genes))
        compound = re.split(r"[-_ ]", name)[0]
        rows.append({
            'signature': name, 'compound': compound,
            'overlap_genes': k, 'signature_size': len(genes),
            'query_size': len(query), 'p_value': p,
            'direction': direction,
            'overlapping_genes': ','.join(sorted(genes & query))
        })
    return pd.DataFrame(rows)


def safe_zscore(s):
    s = np.asarray(s, float)
    if np.std(s) < 1e-9:
        return np.zeros_like(s)
    return (s - np.mean(s)) / np.std(s)


def base_brd(bid):
    m = re.match(r"(BRD-[A-Z0-9]+)", str(bid))
    return m.group(1) if m else None


def load_brd_map():
    brd_map = {}
    if os.path.exists(HUB_MOL):
        hub = pd.read_parquet(HUB_MOL)
        brd_col = next((c for c in hub.columns if re.search(r"broad|brd|deprecated_broad_id|pert_id", c, re.I)), None)
        name_col = next((c for c in hub.columns if re.search(r"pert_iname|name", c, re.I)), None)
        if brd_col and name_col:
            for b, n in zip(hub[brd_col], hub[name_col]):
                bb = base_brd(b)
                if bb and isinstance(n, str) and n:
                    brd_map.setdefault(bb, n)
    return brd_map


def main():
    sig_path = "/workspace/osdr_pipeline/results/signatures/consensus_catabolic_signature.csv"
    out_dir = "/workspace/osdr_pipeline/results/lincs"
    os.makedirs(out_dir, exist_ok=True)

    # 1. Read signature
    up_genes, down_genes = read_signature(sig_path)
    print(f"Signature: {len(up_genes)} up, {len(down_genes)} down genes")

    # 2. Parse local GMT
    drug_sets = parse_gmt(GMT_DRUG)
    all_lincs_genes = set()
    for genes in drug_sets.values():
        all_lincs_genes.update(genes)
    universe = len(all_lincs_genes)
    print(f"LINCS GMT: {len(drug_sets)} signatures, {universe} unique genes")

    up_set, down_set = set(up_genes), set(down_genes)
    up_in_lincs = sorted(up_set & all_lincs_genes)
    down_in_lincs = sorted(down_set & all_lincs_genes)
    print(f"Up genes in LINCS: {len(up_in_lincs)}, Down genes in LINCS: {len(down_in_lincs)}")

    # 3. iLINCS API queries
    known_atrophy_in_lincs = [g for g in KNOWN_ATROPHY if g in all_lincs_genes]
    known_anti_in_lincs = [g for g in KNOWN_ANTI_ATROPHY if g in all_lincs_genes]

    combined_up = sorted(set(up_in_lincs + known_atrophy_in_lincs))
    combined_down = sorted(set(down_in_lincs + known_anti_in_lincs))

    print(f"\nQuerying iLINCS with {len(combined_up)} up genes...")
    df_up = ilincs_query(combined_up)
    print(f"  -> {len(df_up)} signatures, {df_up['compound'].nunique() if len(df_up) else 0} compounds")

    print(f"Querying iLINCS with {len(combined_down)} down/anti-atrophy genes...")
    df_dn = ilincs_query(combined_down)
    print(f"  -> {len(df_dn)} signatures, {df_dn['compound'].nunique() if len(df_dn) else 0} compounds")

    # 4. Local GMT reversal scoring
    print("\nRunning local GMT reversal scoring...")
    fb = gmt_reversal_scoring(drug_sets, up_set, down_set, universe)
    print(f"  -> {len(fb)} signatures with overlap, {fb['compound'].nunique()} compounds")

    # 5. Aggregate and score
    up_agg = df_up.groupby('compound').agg(
        z_mimic_atrophy_max=('zScores', 'max'),
        z_mimic_atrophy_median=('zScores', 'median'),
        n_signatures_up=('zScores', 'count'),
        cell_lines_up=('cellLine', lambda x: ','.join(sorted(set(x)))),
        protein_targets_up=('ProteinTargets', lambda x: next((v for v in x if isinstance(v, str) and v != 'NA'), None))
    ).reset_index()

    dn_agg = df_dn.groupby('compound').agg(
        z_promote_muscle_max=('zScores', 'max'),
        z_promote_muscle_median=('zScores', 'median'),
        n_signatures_dn=('zScores', 'count'),
        cell_lines_dn=('cellLine', lambda x: ','.join(sorted(set(x)))),
        protein_targets_dn=('ProteinTargets', lambda x: next((v for v in x if isinstance(v, str) and v != 'NA'), None))
    ).reset_index()

    gmt_agg = fb.groupby('compound').agg(
        gmt_overlap_count=('overlap_genes', 'sum'),
        gmt_best_p=('p_value', 'min'),
        gmt_n_signatures=('p_value', 'count'),
        gmt_directions=('direction', lambda x: ','.join(sorted(set(x)))),
        gmt_overlapping_genes=('overlapping_genes', lambda x: ','.join(sorted(set(','.join(x)))))
    ).reset_index() if len(fb) > 0 else pd.DataFrame(columns=['compound', 'gmt_overlap_count', 'gmt_best_p', 'gmt_n_signatures'])

    # Merge
    candidates = dn_agg.merge(up_agg[['compound', 'z_mimic_atrophy_max', 'z_mimic_atrophy_median']],
                               on='compound', how='left')
    candidates = candidates.merge(gmt_agg, on='compound', how='left')
    candidates['z_mimic_atrophy_max'] = candidates['z_mimic_atrophy_max'].fillna(0)
    candidates['z_mimic_atrophy_median'] = candidates['z_mimic_atrophy_median'].fillna(0)
    candidates['gmt_overlap_count'] = candidates['gmt_overlap_count'].fillna(0)
    candidates['gmt_best_p'] = candidates['gmt_best_p'].fillna(1)
    candidates['gmt_n_signatures'] = candidates['gmt_n_signatures'].fillna(0)

    # Composite score
    candidates['anti_mimic_score'] = 1 / (1 + candidates['z_mimic_atrophy_max'])
    candidates['gmt_evidence'] = -np.log10(candidates['gmt_best_p'].clip(lower=1e-10))
    candidates['z_promote_norm'] = safe_zscore(candidates['z_promote_muscle_max'])
    candidates['anti_mimic_norm'] = safe_zscore(candidates['anti_mimic_score'])
    candidates['gmt_evidence_norm'] = safe_zscore(candidates['gmt_evidence'])
    candidates['reversal_score'] = (
        0.50 * candidates['z_promote_norm'] +
        0.30 * candidates['anti_mimic_norm'] +
        0.20 * candidates['gmt_evidence_norm']
    )
    candidates = candidates.sort_values('reversal_score', ascending=False).reset_index(drop=True)
    candidates.insert(0, 'rank', candidates.index + 1)

    # 6. Annotate with BRD names and clinical info
    brd_map = load_brd_map()
    candidates['compound_name'] = candidates['compound'].apply(
        lambda x: brd_map.get(base_brd(x), x) if isinstance(x, str) else x
    )

    if os.path.exists(HUB_INFO):
        hub_info = pd.read_parquet(HUB_INFO)
        hub_info['pert_iname_lower'] = hub_info['pert_iname'].str.lower()
        candidates['compound_lower'] = candidates['compound_name'].str.lower()
        candidates = candidates.merge(
            hub_info[['pert_iname_lower', 'clinical_phase', 'moa', 'target', 'disease_area', 'indication']],
            left_on='compound_lower', right_on='pert_iname_lower', how='left'
        )
        candidates = candidates.drop(columns=['compound_lower', 'pert_iname_lower'])

    # 7. Save results
    candidates.to_csv(f"{out_dir}/lincs_reverser_ranking_full.csv", index=False)
    candidates.head(100).to_csv(f"{out_dir}/lincs_reverser_top100.csv", index=False)
    df_up.to_csv(f"{out_dir}/ilincs_up_gene_query_raw.csv", index=False)
    df_dn.to_csv(f"{out_dir}/ilincs_down_gene_query_raw.csv", index=False)
    fb.to_csv(f"{out_dir}/lincs_gmt_reversal_raw.csv", index=False)

    # Summary
    summary = {
        "engine": "ilincs_api + local_gmt",
        "signature": {"up_genes": len(up_genes), "down_genes": len(down_genes),
                      "up_genes_in_lincs": len(up_in_lincs), "down_genes_in_lincs": len(down_in_lincs)},
        "results": {"total_ranked_compounds": len(candidates),
                     "clinically_annotated": int(candidates['clinical_phase'].notna().sum()) if 'clinical_phase' in candidates else 0},
        "scoring": {"method": "composite_zscore",
                    "weights": {"promote_muscle": 0.50, "anti_mimic": 0.30, "gmt_evidence": 0.20}},
    }
    with open(f"{out_dir}/robustness_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nDone! {len(candidates)} compounds ranked. Top 5:")
    for _, row in candidates.head(5).iterrows():
        print(f"  {row['rank']}. {row['compound_name']} (score={row['reversal_score']:.3f})")


if __name__ == "__main__":
    main()
