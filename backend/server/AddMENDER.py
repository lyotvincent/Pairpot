import scanpy as sc
from normalize import mender
import  anndata as ad
def AddMender4SPdata(spH5adFile):
  adata = sc.read_h5ad(spH5adFile)
  if 'leiden-1' in adata.obs.columns:
    if adata and isinstance(adata,ad.AnnData) and 'spatial' in adata.obsm:
      adata = mender(adata)
      adata.write_h5ad(spH5adFile)