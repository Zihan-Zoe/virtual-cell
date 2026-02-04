# Virtual Cell for Perturb-seq

This project studies gene-level and combinatorial perturbation effects
using single-cell Perturb-seq data.

## Datasets
- Norman Perturb-seq (dual-gene perturbations)
- Replogle K562 / RPE1 (single-gene perturbations)

Raw `.h5ad` files are not included.
Please download them separately and place under `data/`.

## Environment
```bash
conda env create -f environment.yml
conda activate virtualcell
