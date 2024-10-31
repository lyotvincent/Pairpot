import scanpy as sc
from normalize import rank
import numpy as np
def AddRank4SPdata(spH5adFile,organs, clu_key="leiden-1"):
    adata=sc.read_h5ad(spH5adFile)
    if "annotation" not in adata.obs.columns:
        adata.obs["annotation"]=adata.obs["leiden-1"]
    col = adata.obs.columns
    dorank = np.sum(list(map(lambda x: x.startswith("UCell"), col)))
    if dorank == 0: # no UCell
        adata = rank(adata=adata,organs=organs, clu_key=clu_key)
        adata.write_h5ad(spH5adFile)

def AddRank4Output(dataset_id, scdata_id):
    adata_spFile = f"/data/rzh/RawUrls/{dataset_id}/STDS0000{dataset_id}/New_{dataset_id}.h5ad"
    adata_scFile = f"/data/rzh/RawUrls/{dataset_id}/SCDS0000{scdata_id}/New_{scdata_id}.h5ad"
    adata_outPath = f"/data/rzh/RawUrls/{dataset_id}"

    # add rank to sc_sampled.h5ad
    adata_sc_sampled = sc.read_h5ad(f"{adata_outPath}/sc_sampled.h5ad")
    adata_sc = sc.read_h5ad(adata_scFile)
    adata_sc_sampled.obs = adata_sc.obs.loc[adata_sc_sampled.obs_names, :]
    adata_sc_sampled.write_h5ad(f"{adata_outPath}/sc_sampled.h5ad")
    print("** Add rank to sc_sampled.h5ad")

    # add rank to sp_deconv.h5ad
    adata_sp_deconv = sc.read_h5ad(f"{adata_outPath}/sp_deconv.h5ad")
    adata_sp = sc.read_h5ad(adata_spFile)
    adata_sp_deconv.obs = adata_sp.obs.loc[adata_sp_deconv.obs_names, :]
    adata_sp_deconv.write_h5ad(f"{adata_outPath}/sp_deconv.h5ad")
    print("** Add rank to sp_deconv.h5ad")