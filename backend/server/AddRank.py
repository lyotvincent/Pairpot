import scanpy as sc
from normalize import rank
import numpy as np
def AddRank4SPdata(spH5adFile,organs, clu_key="leiden-1"):
    adata=sc.read_h5ad(spH5adFile)
    col = adata.obs.columns
    dorank = np.sum(list(map(lambda x: x.startswith("UCell"), col)))
    if dorank == 0: # no UCell
        adata = rank(adata=adata,organs=organs, clu_key=clu_key)
        adata.write_h5ad(spH5adFile)