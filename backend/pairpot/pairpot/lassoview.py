import numpy as np
import scanpy as sc
import pandas as pd
import random
from collections import Counter
from sklearn.semi_supervised import LabelPropagation, LabelSpreading
import h5py
import scipy.sparse
import time
import random
import numpy as np
import pkg_resources
import sys
import os
import anndata as ad
# lpa_so_path = pkg_resources.resource_filename('pairpot', 'label_propagation.cpython-310-x86_64-linux-gnu.so')
# print(lpa_so_path)
# sys.path.append(os.path.dirname(lpa_so_path))
# import label_propagation as LPA
import pairpotlpa as LPA

from .normalize import *

def EagerRefine(selected, file):
    adata = sc.read_h5ad(file)
    obs = adata.obs.copy()
    obs = obs.reset_index()
    selectedObs = obs.iloc[selected, :]
    selectedClu = np.unique(selectedObs['leiden'])
    refinedObs = obs.loc[obs['leiden'].isin(selectedClu), :]
    return list(refinedObs.index)


def LazyRefine(selected, file):
    adata = sc.read_h5ad(file)
    obs = adata.obs.copy()
    obs = obs.reset_index()
    obsCounter = Counter(obs['leiden'])

    selectedObs = obs.iloc[selected, :]
    selectedCounter = Counter(selectedObs['leiden'])

    selectedClu = []
    for key in selectedCounter:
        if selectedCounter[key] / obsCounter[key] > 0.5:
            selectedClu.append(key)

    refinedObs = obs.loc[obs['leiden'].isin(selectedClu), :]
    return list(refinedObs.index)


# def LPARefine(selected,  file, use_model=LabelPropagation):
#     adata = sc.read_h5ad(file)
#     X = adata.obsm['X_umap']
#     adata.obs['test'] = adata.obs['annotation'].iloc[
#         np.random.choice(len(adata), int(len(adata) * 0.1), replace=False)]
#     y = adata.obs['test'].copy()
#     oriCat = list(adata.obs['test'].cat.categories)
#     oriCat.append("Selected")
#     y = y.cat.rename_categories(list(range(len(oriCat)-1)))
#     y = np.array([-1 if s is np.NaN else s for s in y])
#     y[selected] = len(oriCat)
#     # do LabelPropagation
#     y = pd.Series(y, index=adata.obs_names)
#     model = use_model().fit(X, y)
#     y_pred = model.predict(X)
#     y_pred = pd.Series(y_pred)
#     y_pred = y_pred[y_pred == len(oriCat)]
#     return list(y_pred.index)

def LPARefine(adata, selected, function="anndata", use_model=LabelPropagation, do_correct=True):
    # mat=1
    # with h5py.File(file,'r') as f:
    #     group=f['obsp']['connectivities']

    #     data=group['data'][:]
    #     indices=group['indices'][:]
    #     indptr=group['indptr'][:]
    #     shape=(f['obsp']['connectivities'].attrs['shape'][0],f['obsp']['connectivities'].attrs['shape'][1])

    #     mat=scipy.sparse.csr_matrix((data,indices,indptr),shape=shape)
    
    # adata = ad.read_h5ad(file)
    # if "connectivities" not in adata.obsp:
    #     print("No probability transition matrix, attempting to generate.")
    #     #Generate probability transition matrix
        
    #     directory = os.path.dirname(file)
    #     # print(directory+"/New_"+directory[-3:]+".h5ad")
    #     # samples=filepath
    #     # sampleNames=filename
    #     # adata=concat_adata(filepath, filename, inputFunc=input_adata_h5ad)
    #     # print(adata)
    #     adata = pp(adata)
    #     if adata.shape[0]>100000 and type=='sc':
    #         sampleIdx = np.random.choice(len(adata), 50000)
    #         adata = adata[sampleIdx,:].copy()
    #     adata = clu(adata)
    #     if type=='sp':
    #         print("doing MENDER.")
    #         adata = mender(adata)
    #         clu_key = 'annotation'
    #     else:
    #         clu_key = 'leiden-1'
        
    #     try:
    #         adata = rank(adata, organs, clu_key=clu_key)
    #     except:
    #         try:
    #             print("Try kruskal")
    #             adata=rank(adata,organs,test_func=kruskal)
    #         except:
    #             print("Error in to kruskal")
    #     # print(adata.uns['UCell_Assign'])
    #     adata = marker(adata, groupby="leiden-1", method='wilcoxon')
    #     adata.obs = adata.obs.astype('str').fillna("")
    #     # adata = anno(adata, annoDict)
    #     # sc.pl.umap(adata, color='annotation')
    #     # adata = marker(adata)
    #     # print(adata)
    #     adata.write_h5ad(directory+"/New_"+directory[-3:]+".h5ad")
    #     print("Processed data has saved to"+directory+"/New_"+directory[-3:]+".h5ad")
    mat=1
    if function=="anndata":
        mat=adata.obsp['connectivities']
        if not scipy.sparse.issparse(mat):
            # print("1111111111111111")
            mat=scipy.sparse.csr_matrix(mat)
    elif function=="h5adfile":
        with h5py.File(adata,'r') as f:
            group=f['obsp']['connectivities']

            data=group['data'][:]
            indices=group['indices'][:]
            indptr=group['indptr'][:]
            shape=(f['obsp']['connectivities'].attrs['shape'][0],f['obsp']['connectivities'].attrs['shape'][1])

            mat=scipy.sparse.csr_matrix((data,indices,indptr),shape=shape)
    else:
        print('No this function, use "anndata" or "h5adfile" instead.')
        return
    coo=mat.tocoo()

    rows=coo.row
    cols=coo.col
    data=coo.data

    if function=="anndata":
        obs_col = 'annotation'
        if obs_col not in adata.obs:
            obs_col = 'leiden-1'

        if "codes" in adata.obs[obs_col]:
            mat = adata.obs[obs_col]['codes'].values
        else:
            mat = adata.obs[obs_col].values
    elif function=="h5adfile":
        with h5py.File(adata, 'r') as h5file:
            obs_group = h5file['obs']
            obs_col = 'annotation'
            if obs_col not in obs_group:
                obs_col = 'leiden-1'

            if "codes" in obs_group[obs_col]:
                mat = obs_group[obs_col]['codes'][:]
            else:
                mat = obs_group[obs_col][:]
    else:
        print('No this function, use "anndata" or "h5adfile" instead.')
        return
    val={}

    for i in np.unique(mat):
        val[i]=len(val)
    val[len(val)] = len(val)
    X = LPA.matCoo(mat.shape[0], mat.shape[0])
    for i in range(len(data)):
        X.append(rows[i], cols[i], data[i])
                
    y_label = LPA.mat(mat.shape[0], len(val))
    random_list=random.sample(range(mat.shape[0]), int(mat.shape[0] * 0.1))
    select_list=np.zeros(mat.shape[0])
    y_label.setneg()
    select_list[random_list] = 1

    # add selected item
    select_list[selected] = 1
    selected_val = len(val) - 1

    # if function=="h5adfile":
    #     mat[selected] = selected_val
    # elif function=="anndata":
    mat_list = mat.tolist()
    for t in range(len(selected)):
        mat_list[selected[t]]=selected_val
    mat = pd.Categorical(mat_list)
    for i in range(mat.shape[0]):
        if select_list[i]:
            y_label.editval2(i,val[mat[i]])
            # print(f"i:{i}, val:{y_label.getval(i,selected_val)}")
    # print("y_label")
    # for i in range(100):
    #     print(f"i:{i},val:{y_label.getval(i,selected_val)}")
    y_pred = LPA.mat(mat.shape[0], len(val))
    y_new = LPA.mat(mat.shape[0], len(val))
    LPA.labelPropagation(X, y_label, y_pred, y_new, 0.5,1000)
    y_res = np.zeros(mat.shape[0])
    if do_correct:
        for i in range(mat.shape[0]):
            y_res[i] = y_new.getval(i,0)
            # if i in selected:
                # print(f"i:{i}, val:{y_res[i]}")
    else:
        for i in range(mat.shape[0]):
            y_res[i] = y_pred.getval(i,0)
    # print(f"selected_val:{selected_val}")
    # for i in range(mat.shape[0]):
    # print("y_new")
    # for i in range(100):
    #     print(f"i:{i},val:{y_new.getval(i,0)}")
    # print("y_pred")
    # for i in range(100):
    #     print(f"i:{i},val:{y_pred.getval(i,0)}")
    y_res = pd.Series(y_res)
    y_res = y_res[y_res == selected_val]
    # print("1111111111111111")
    return list(y_res.index)


def Gen_maskSet(candidate: pd.DataFrame, errRate=0.20):
    sele_can = candidate[candidate==True]
    cell_len = len(sele_can)
    mask_can = candidate.copy()
    errArray = random.sample(range(cell_len), int(cell_len*errRate))
    for cell in errArray:
        print(sele_can.index[cell])
        mask_can.loc[sele_can.index[cell]] = not mask_can.loc[sele_can.index[cell]]
    return mask_can

def train_LPA(candidate, adata, use_model=LabelPropagation, errRate=0.05):
    X = adata.obsm['X_pca']
    y = Gen_maskSet(candidate, errRate)
    model = use_model().fit(X, y)
    y_pred = model.predict(X)

    def acc(y_true, y_pred):
        return np.sum(np.equal(y_true, y_pred))/len(y_true)
    print("acc:{}".format(acc(candidate, y_pred)))

