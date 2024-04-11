from sklearn.metrics import roc_auc_score
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from scipy.stats import kruskal

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
    intagIdx = mat[mat.isin(markerIdx)]
    if len(intagIdx) == 0:
       ucell[i]=0
    else:
      u = np.sum([list(mat).index(s) for s in intagIdx]) + (n-len(intagIdx))*maxRank - smin
      ucell[i] = 1 - u / umax

  # do knn-smooth
  # knnMat = np.argsort(-adata.obsp["connectivities"].todense(), axis=1)[:,:5]
  # for i in range(len(ucell)):
  #   knnCells = list(knnMat[i])
  #   ucell[i] = np.mean(ucell[knnCells])
  adata.obs[f"UCell_{cellType}"] = pd.Series(ucell, index=adata.obs_names, dtype=float)
  return adata



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
      w = kruskal(sample1, sample2)
      if w.pvalue < alpha:
        candidates.append(anno)
        ucell = ucell[~ucell['leiden-1'].isin([anno])]
      else:
        break
    if len(candidates) > 0:
      annotation[ct] = candidates
  adata.uns['UCell_Assign'] = annotation
  adata.obsm['AUCell_rankings'].columns = np.array(adata.obsm['AUCell_rankings'].columns, dtype=str)
  return adata