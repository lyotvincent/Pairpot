import scanpy as sc
import numpy as np
import pandas as pd
import h5py as h5
from sklearn.linear_model import LinearRegression
# import AUCell as al

def NNLSDeconvolution(selected, scdata, spdata, function="anndata"):
  # adata_sc = al.AUCell_UCAssign(sc.read_h5ad(scfile),n_jobs=8)
  # adata_sp = al.AUCell_UCAssign(sc.read_h5ad(spfile),n_jobs=8)
  adata_sc=1
  adata_sp=1
  if function=="anndata":
    adata_sc=scdata
    adata_sp=spdata
  elif function=="h5adfile":
    adata_sc = sc.read_h5ad(scdata)
    adata_sp = sc.read_h5ad(spdata)
  else:
    print('No this function, use "anndata" or "h5adfile" instead.')
    return
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

