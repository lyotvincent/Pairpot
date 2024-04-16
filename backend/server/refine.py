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
import label_propagation as LPA

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

def LPARefine(selected,  file, use_model=LabelPropagation, do_correct=True):
    with h5py.File(file,'r') as f:
        group=f['obsp']['connectivities']

        data=group['data'][:]
        indices=group['indices'][:]
        indptr=group['indptr'][:]
        shape=(f['obsp']['connectivities'].attrs['shape'][0],f['obsp']['connectivities'].attrs['shape'][1])

        mat=scipy.sparse.csr_matrix((data,indices,indptr),shape=shape)
    
    coo=mat.tocoo()

    rows=coo.row
    cols=coo.col
    data=coo.data

    with h5py.File(file, 'r') as h5file:
        obs_group = h5file['obs']
        if "codes" in obs_group['annotation']:
            mat = obs_group['annotation']['codes'][:]
        else:
            mat = obs_group['annotation'][:]
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
    mat[selected] = selected_val
    for i in range(mat.shape[0]):
        if select_list[i]:
            y_label.editval2(i,val[mat[i]])

    y_pred = LPA.mat(mat.shape[0], len(val))
    y_new = LPA.mat(mat.shape[0], len(val))
    LPA.labelPropagation(X, y_label, y_pred, y_new, 0.5,1000)
    y_res = np.zeros(mat.shape[0])
    if do_correct:
        for i in range(mat.shape[0]):
            y_res[i] = y_new.getval(i,0)
    else:
        for i in range(mat.shape[0]):
            y_res[i] = y_pred.getval(i,0)
    y_res = pd.Series(y_res)
    y_res = y_res[y_res == selected_val]
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

