from sklearn.metrics import roc_auc_score
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from scipy.stats import kruskal,mannwhitneyu
from joblib import Parallel, delayed

def partition_arg_topK(matrix, K, axis=1):
    """
    perform topK based on np.argpartition
    :param matrix: to be sorted
    :param K: select and sort the top K items
    :param axis: 0 or 1. dimension to be sorted.
    :return:
    """
    a_part = np.argpartition(matrix, K, axis=axis)
    if axis == 0:
        row_index = np.arange(matrix.shape[1 - axis])
        a_sec_argsort_K = np.argsort(matrix[a_part[0:K, :], row_index], axis=axis)
        return a_part[0:K, :][a_sec_argsort_K, row_index]
    else:
        column_index = np.arange(matrix.shape[1 - axis])[:, None]
        a_sec_argsort_K = np.argsort(matrix[column_index, a_part[:, 0:K]], axis=axis)
        return a_part[:, 0:K][column_index, a_sec_argsort_K]
    
def naive_topK(matrix, K, axis=1):
    mat = np.argsort(-matrix, axis=axis)
    return mat[:, :K]
    

def AUCell_buildRankings(adata:ad.AnnData, top=0.05):
  k = int(len(adata.var_names)*top)
  adata.obsm['AUCell_rankings'] = pd.DataFrame(naive_topK(adata.X.todense(), k), index=adata.obs_names)
  adata.obsm['AUCell_rankings'].columns = np.array(adata.obsm['AUCell_rankings'].columns, dtype=str)
  return adata


def AUCell_calcAUC(adata:ad.AnnData, markerList:list, cellType:str, rankings="AUCell_rankings"):
  markerSet = list(set(markerList).intersection(set(adata.var_names)))
  if rankings in adata.obsm:
    y_score = list(range(len(adata.obsm['AUCell_rankings'].columns)))
    y_score.reverse()
    aucell = np.zeros_like(adata.obs_names)
    for i in range(len(adata.obsm['AUCell_rankings'])):
      y_test = adata.var_names[adata.obsm['AUCell_rankings'].iloc[i,:]].isin(markerSet)
      if sum(y_test) == 0:
        aucell[i] = 0
      else:
        aucell[i] = roc_auc_score(y_test, y_score)
    adata.obs[f"AUCell_{cellType}"] = aucell
  else:
    print(f"{rankings} not found in adata.obsm, run AUCell_buildRankings first.")
  return adata


def AUCell_exploreThreshold(adata:ad.AnnData, cellType:str, assign=True, index="AUCell"):
  aucell = adata.obs[f'{index}_{cellType}']
  bins = np.array(range(10))/10 * np.max(aucell)
  hist = np.histogram(aucell, bins=bins)[0]
  total = len(aucell)
  mean = np.mean(aucell)

  w0, u0, w1, u1, u = 0, 0, 0, 0, 0
  max_variance = 0.0
  threshold = 0
  for i,t in enumerate(bins):
    # 阈值为t时的类间方差计算
    w0 = np.sum(hist[:i]) / total
    w1 = 1 - w0
    if w0 == 0 or w1 == 0:
        continue
    
    u0 = np.sum(hist[:i] * bins[1:i+1]) / w0
    u1 = np.sum(hist[i:] * bins[i+1:]) / w1
    u = u0 * w0 + u1 * w1
    # 类内方差
    var_b = w0 * (u0 - mean) ** 2 + w1 * (u1 - mean) ** 2
    if var_b > max_variance:
        max_variance = var_b
        threshold = t
  
  # add to adata.uns
  if 'AUCThreshold' not in adata.uns:
    adata.uns['AUCThreshold'] = {}
  adata.uns['AUCThreshold'][cellType] = threshold
  if assign:
    assign = aucell[aucell >= threshold] 
    assign.iloc[:] = cellType
    adata.obs[f'{index}_{cellType}_Assignment'] = assign
  print(f"threshold of {cellType} is {threshold}")
  return adata


def AUCell_calcUC(adata:ad.AnnData, markerList:list, cellType:str, rankings="AUCell_rankings"):
  varList = list(adata.var_names)
  markerIdx = [varList.index(s) for s in markerList]
  rankMat = adata.obsm[rankings]
  maxRank = len(adata.obsm[rankings].columns)
  n = len(markerIdx)
  smin = n*(n+1)/2
  smax = n*maxRank
  umax = smax - smin
  ucell = np.zeros_like(adata.obs_names)
  for i in range(len(rankMat)):
    mat = rankMat.iloc[i, :]
    intagIdx = mat.isin(markerIdx)
    sum_intagIdx = np.sum(intagIdx)
    if sum_intagIdx == 0:
       ucell[i]=0
    else:
      u = np.sum(np.where(intagIdx)[0])+ (n-sum_intagIdx)*maxRank - smin
      ucell[i] = 1 - u / umax
  return ucell



def AUCell_UCAssign(adata:ad.AnnData, db:pd.DataFrame, celltype:str, alpha=10e-30, gene_col='official gene symbol'):
  annotation = {}
  for ct in celltype:
    candidates = []
    markerList = np.array(db[db['cell type'] ==ct][gene_col])
    markerList = list(set(markerList).intersection(set(adata.var_names)))
    adata = AUCell_calcUC(adata, markerList, ct)
    ucell = adata.obs[["leiden-1", f"UCell_{ct}"]]
    rank = ucell.groupby("leiden-1", observed=False).mean()
    rank = rank.sort_values(by=f"UCell_{ct}", ascending=False)
    rank = rank.reset_index()

    for i in range(len(rank)):
      anno = rank.iloc[i,0]
      sample1 = ucell.loc[ucell['leiden-1'].isin([anno]), f"UCell_{ct}"]
      sample2 = ucell.loc[~ucell['leiden-1'].isin([anno]), f"UCell_{ct}"]
      u_stat,p_val = kruskal(sample1, sample2)
      if p_val < alpha:
        candidates.append(anno)
        ucell = ucell[~ucell['leiden-1'].isin([anno])]
      else:
        break
    if len(candidates) > 0:
      annotation[ct] = candidates
  adata.uns['UCell_Assign'] = annotation
  adata.obsm['AUCell_rankings'].columns = np.array(adata.obsm['AUCell_rankings'].columns, dtype=str)
  return adata


def AUCell_UCAssign(adata, 
                    db:pd.DataFrame, 
                    celltype:str, 
                    alpha=10e-30, 
                    gene_col='official gene symbol', 
                    test_func=kruskal,
                    clu_key="leiden-1",
                    n_jobs=8):
  
  # calculate UCell
  def AUCell_calcUC_thread(i):
    ct = celltype[i]
    markerList = np.array(db[db['cell type'] ==ct][gene_col])
    markerList = list(set(markerList).intersection(set(adata.var_names)))
    ucell = AUCell_calcUC(adata, markerList, ct)
    return ucell
  ucells = Parallel(n_jobs=n_jobs)(delayed(AUCell_calcUC_thread)(i) for i in range(len(celltype)))
  ucells_col = [f"UCell_{ct}" for ct in celltype]
  ucells_df = pd.DataFrame(ucells, index=ucells_col, columns=adata.obs_names).T
  for ct in ucells_col:
    adata.obs[ct] = ucells_df[ct].astype("float")
  
  # Assign annotations by UCell
  def UCAssign_Thread(i):
    ct = celltype[i]
    candidates = []
    ucell = adata.obs[[clu_key, f"UCell_{ct}"]].copy()
    rank = ucell.groupby(clu_key, observed=False).mean()
    rank = rank.sort_values(by=f"UCell_{ct}", ascending=False)
    rank = rank.reset_index()
    assert len(rank) == len(adata.obs[clu_key].unique())

    for i in range(len(rank)):
      anno = rank.iloc[i,0]
      sample1 = ucell.loc[ucell[clu_key] == anno, f"UCell_{ct}"]
      sample2 = ucell.loc[ucell[clu_key] != anno, f"UCell_{ct}"]
      if len(sample1) > 0  and len(sample2) > 0:
        u_stat,p_val = test_func(sample1, sample2)
      else:
        p_val=1
      if p_val < alpha:
        candidates.append(anno)
        ucell = ucell[ucell[clu_key] != anno]
      else:
        break
    return candidates
  can = Parallel(n_jobs=n_jobs)(delayed(UCAssign_Thread)(i) for i in range(len(celltype)))
  annotation = dict(zip(celltype, can))
  adata.uns['UCell_Assign'] = annotation
  return adata