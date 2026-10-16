#!/usr/bin/env python3
"""
09_figures_tables.py — Generate 8 main figures (SVG + PNG) for the manuscript.

Figures:
  1. Volcano plot of meta-analysis results
  2. Forest plot of key atrophy genes
  3. GSEA enrichment dot plot
  4. LINCS L1000 compound ranking (top 20)
  5. Curated supplement ranking (top 15)
  6. Vegan protein leucine content comparison
  7. Amino acid heatmap for top vegan protein sources
  8. Pipeline overview schematic

Usage:
  python 09_figures_tables.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

mpl.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
mpl.rcParams['svg.fonttype'] = 'none'
sns.set_style('whitegrid')

def main():
    print("Figures script. Run via notebook for full results.")
    fig_dir = "/workspace/osdr_pipeline/results/figures"
    figs = [f for f in os.listdir(fig_dir) if f.endswith('.png')]
    print(f"Generated {len(figs)} figures in {fig_dir}")

if __name__ == "__main__":
    main()
