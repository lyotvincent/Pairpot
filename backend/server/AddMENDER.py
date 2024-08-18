import scanpy as sc
from normalize import mender

def AddMender4SPdata(spH5adFile):
  adata = sc.read_h5ad(spH5adFile)
  if 'leiden-1' in adata.obs.columns:
    adata = mender(adata)
    adata.write_h5ad(spH5adFile)