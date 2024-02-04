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


def LPARefine(selected, refered, file, use_model=LabelPropagation):
    adata = sc.read_h5ad(file)
    X = adata.obsm['X_pca']
    y = np.ones_like(adata.obs_names, dtype=int) * -1
    y[selected] = 1
    y[refered] = 0
    y = pd.Series(y, index=adata.obs_names)
    print(y)
    y_orig = np.zeros_like(adata.obs_names, dtype=int)
    y_orig[selected] = 1
    model = LabelPropagation().fit(X, y)
    y_pred = model.predict(X)
    print(y_pred)
    def acc(y_true, y_pred):
        return np.sum(np.equal(y_true, y_pred))
    print("similar:{}".format(acc(y_orig, y_pred)))
    return y_pred


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

