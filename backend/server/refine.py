import numpy as np
import scanpy as sc
import pandas as pd
import random
from collections import Counter
from sklearn.semi_supervised import LabelPropagation, LabelSpreading

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


def LPARefine(selected, file, use_model=LabelPropagation):
    adata = sc.read_h5ad(file)
    X = adata.obsm['X_umap']
    adata.obs['test'] = adata.obs['annotation'].iloc[
        np.random.choice(len(adata), int(len(adata) * 0.1), replace=False)]
    y = adata.obs['test'].copy()
    oriCat = list(adata.obs['test'].cat.categories)
    oriCat.append("Selected")
    y = y.cat.rename_categories(list(range(len(oriCat)-1)))
    y = np.array([-1 if s is np.NaN else s for s in y])
    y[selected] = len(oriCat)
    # do LabelPropagation
    y = pd.Series(y, index=adata.obs_names)
    model = use_model().fit(X, y)
    y_pred = model.predict(X)
    y_pred = pd.Series(y_pred, index=adata.obs_names)
    y_pred = y_pred[y_pred == len(oriCat)]
    return y_pred.index


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

