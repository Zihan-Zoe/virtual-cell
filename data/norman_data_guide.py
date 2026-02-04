"""
Norman Perturb-seq Dataset
===============================================================

Simulation split test composition:
simulation: mixed of all types below
combo_seen0: 0      # Dual perts where neither gene seen in training
combo_seen1: 0      # Dual perts where one gene seen in training  
combo_seen2: 0      # Dual perts where both genes seen separately
"""

import scanpy as sc

# =============================================================================
# 1. LOAD DATA
# =============================================================================
H5AD_PATH = "./data/norman.h5ad"
adata = sc.read_h5ad(H5AD_PATH)

print(f"""
================================================================================
DATASET OVERVIEW
================================================================================
Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes

Key fields:
  adata.X                    → Expression matrix (log-normalized, sparse)
  adata.obs['condition']     → Perturbation label (e.g., "CBL+ctrl", "CBL+UBASH3B")
  adata.obs['control']       → 1 = control cell, 0 = perturbed cell
  adata.var['gene_name']     → Gene names
""")

# =============================================================================
# 2. PERTURBATION TYPES
# =============================================================================
conditions = list(adata.obs['condition'].unique())
single_perts = [c for c in conditions if '+ctrl' in c or 'ctrl+' in c]
dual_perts = [c for c in conditions if '+' in c and 'ctrl' not in c]
n_ctrl_cells = (adata.obs['control'] == 1).sum()

print(f"""
================================================================================
PERTURBATION STRUCTURE
================================================================================
Control cells:              {n_ctrl_cells:,} cells (baseline, no perturbation)
Single-gene perturbations:  {len(single_perts)} conditions (e.g., "CBL+ctrl")
Dual-gene perturbations:    {len(dual_perts)} conditions (e.g., "CBL+UBASH3B")

Examples:
  Single: {single_perts[:3]}
  Dual:   {dual_perts[:3]}
""")

# =============================================================================
# 3. ML TASK FORMULATION
# =============================================================================
print("""
================================================================================
ML TASK: PERTURBATION PREDICTION
================================================================================
Goal: Given a control cell + perturbation label → predict post-perturbation expression

Two common formulations:

  (A) Predict absolute expression:
      Input:  control_expr (n_genes,) + pert_label
      Output: perturbed_expr (n_genes,)

  (B) Predict delta (often works better):
      Input:  control_expr (n_genes,) + pert_label  
      Output: delta = perturbed_expr - control_expr (n_genes,)

""")

# ================================================================================
# USING PRE-COMPUTED SPLITS (NO GEARS DEPENDENCY)
# ================================================================================
# load splits with plain pandas:

import pandas as pd

# Load split masks
split_type = 'simulation'  # NOTE and 'combo_seen0', 'combo_seen1', 'combo_seen2'
masks = pd.read_csv(f'./data/splits_norman/split_{split_type}.csv')

# Subset AnnData
train_adata = adata[masks['train'].values]
val_adata = adata[masks['val'].values]
test_adata = adata[masks['test'].values]



