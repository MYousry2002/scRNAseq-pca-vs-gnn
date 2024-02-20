import anndata as ad
import pandas as pd
import numpy as np

adata = ad.read_h5ad('~/standard-pipeline/GSE198623_raw_data.h5ad')

# Calculate the number of cells to sample for 33%
num_cells_to_sample = int(adata.n_obs * 0.33)

# Generate random indices
random_indices = np.random.choice(adata.obs_names, size=num_cells_to_sample, replace=False)

# Subset the anndata object
adata_subset = adata[random_indices, :]

# save as a csv file
adata_subset.to_df().T.to_csv('~/gnn-pipeline/GSE198623/GSE198623_submatrix.csv')
