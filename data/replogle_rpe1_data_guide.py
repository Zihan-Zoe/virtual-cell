"""
Replogle RPE1 Essential Perturb-seq Dataset
===============================================================

simulation: predict the effect of perturbing gene X when it is not seen in training
"""

import scanpy as sc

# =============================================================================
# 1. LOAD DATA
# =============================================================================
H5AD_PATH = "./data/replogle_rpe1.h5ad"
adata = sc.read_h5ad(H5AD_PATH)

print(f"""
================================================================================
DATASET OVERVIEW
================================================================================
Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes

Key fields:
  adata.X                    → Expression matrix (log-normalized, sparse)
  adata.obs['condition']     → Perturbation label (e.g., "ACTB+ctrl", "ctrl")
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
Single-gene perturbations:  {len(single_perts)} conditions
Dual-gene perturbations:    {len(dual_perts)} conditions

NOTE: Replogle RPE1 dataset contains ONLY single-gene perturbations (no dual-gene combos)

Examples:
  Single: {single_perts[:5]}
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

# Load split masks (only 'simulation' split for Replogle RPE1)
split_type = 'simulation'
masks = pd.read_csv(f'./data/splits_replogle_rpe1/split_{split_type}.csv')

# Subset AnnData
train_adata = adata[masks['train'].values]
val_adata = adata[masks['val'].values]
test_adata = adata[masks['test'].values]
