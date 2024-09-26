import scanpy as sc
import numpy as np
import pandas as pd
import h5py as h5
from sklearn.linear_model import LinearRegression

def NNLSDeconvolution(selected, scfile, spfile):
  adata_sc = sc.read_h5ad(scfile)
  adata_sp = sc.read_h5ad(spfile)
  adata_sc.obs['annotation'] = list(adata_sc.obs['annotation'])
  idx = adata_sc.obs.columns.get_loc('annotation')
  adata_sc.obs.iloc[selected, idx]  = "Selected"
  # generate A
  ucells_sc = [s for s in adata_sc.obs.columns if s.startswith('UCell_')]
  ucells_sp = [s for s in adata_sp.obs.columns if s.startswith('UCell_')]
  ucells = list(set(ucells_sc).intersection(ucells_sp))
  df_sc = adata_sc.obs[ucells].astype("float")
  df_sc['annotation'] = adata_sc.obs['annotation']
  g = df_sc.groupby(by='annotation')
  anno = list(g.groups.keys())
  A = g.mean().T

  # generate b
  df_sp = adata_sp.obs[ucells]
  B = df_sp.T

  # NNLS model
  reg_nnls = LinearRegression(positive=True)
  pred_nnls = reg_nnls.fit(A, B)
  
  return list(pred_nnls.coef_[:, anno.index("Selected")])
