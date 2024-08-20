import scanpy as sc
from normalize import rank
import anndata as ad
def AddRank4SPdata(spH5adFile,organs):
    adata=sc.read_h5ad(spH5adFile)
    adata = rank(adata=adata,organs=organs)
    adata.write_h5ad(spH5adFile)