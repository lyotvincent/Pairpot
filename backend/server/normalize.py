# 将所有的输入数据标准化成统一格式，其中标准化格式如下：
# 单细胞数据
# -表达矩阵
# --下采样得到3000个细胞
# --2000个与空转交集的高变基因
# -低维表示
# --PCA前两维表示
# --UMAP2维和3维表示
# --TSNE2维表示
# -细胞类型标注
# --cellTypist标注
# --SingleR标注
# --原始数据自带标注
# --细胞类型marker基因
# --KEGG结果
# --GSEA结果
# --GO结果
# -邻接矩阵
# --k近邻矩阵
import scanpy as sc
import os
import pandas as pd
import numpy as np
import anndata as ad
import celltypist

def input_adata_10X(sample):
    adata = sc.read_mtx(sample+'/matrix.mtx.gz')
    adata = adata.T
    bar = pd.read_csv(sample+'/barcodes.tsv.gz', header=None)
    fea = pd.read_csv(sample+'/features.tsv.gz', header=None, sep='\t')
    bar.columns = ['barcodes']
    fea.columns = ['ID', 'name', 'type']
    adata.obs_names = bar.iloc[:,0]
    adata.obs_names_make_unique()
    adata.var_names = fea.iloc[:,1]
    adata.var = fea
    adata.var_names_make_unique()
    adata
    return adata

def input_adata_10Xh5(sample):
    adata = sc.read_10x_h5(sample)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata

# 将所有矩阵合并
def concat_adata(samples, sampleNames, inputFunc=input_adata_10Xh5):
    adatas = []
    for i in range(len(sampleNames)):
        adata = inputFunc(samples[i])
        adatas.append(adata)
    # 进行数据合并
    adata_concat = ad.concat(adatas, label="batch", keys=sampleNames, index_unique='-')
    adata_concat
    return adata_concat

# 预处理
def pp(adata):
    mito_genes = adata.var_names.str.startswith('MT-')
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs['mt_frac'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    
    # 过滤低表达的基因
    sc.pp.filter_cells(adata, min_genes=5)  # 过滤一个细胞中表达少于五个基因的细胞样本 
    sc.pp.filter_genes(adata, min_cells=5)  # 过滤在少于五个细胞中表达的基因
    sc.pp.filter_cells(adata, min_counts=30)   # 过滤每个细胞中计数少于29个的细胞样本 

    # 过滤线粒体核糖体基因
    rp_genes = adata.var_names.str.startswith('RP')
    mt_genes = adata.var_names.str.startswith('MT-')
    adata = adata[:, ~(rp_genes + mt_genes)]
    
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)

    adata = adata[adata.obs['mt_frac'] < 0.2]
    return adata

# macro_Mye = scv.read('velo-Macro.h5ad')
def clu(adata, key_added="leiden-1", n_neighbors=50, n_pcs=30, rep='X_pca_harmony', do_har=True, max_iter=20, resolution=1, do_scrublet=True, har_key='batch'):
    # Computing the neighborhood graph
    if do_scrublet:
        n0 = adata.shape[0]
        print("{0} Cell number: {1}".format(key_added, n0))
        sc.pp.scrublet(adata, random_state=112)
        adata = adata[adata.obs['predicted_doublet']==False,:].copy()
        print("{0} Cells retained after scrublet, {1} cells reomved.".format(adata.shape[0], n0-adata.shape[0]))
    else:
        print("Ignoring processing doublet cells...")
    sc.tl.pca(adata, svd_solver='arpack')
    if do_har:
        sc.external.pp.harmony_integrate(adata, key=har_key,max_iter_harmony=max_iter)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep)
    # Run UMAP
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added=key_added, resolution=resolution)
    sc.pl.umap(adata, color=key_added, legend_fontoutline=True, palette=sc.pl.palettes.default_20, legend_loc="on data")
    return adata

def marker(adata, groupby="leiden-1", method='wilcoxon'):
    sc.tl.rank_genes_groups(adata, groupby = groupby, method = method)
    sc.tl.dendrogram(adata, groupby=groupby, use_rep='X_pca_harmony')
    sc.pl.rank_genes_groups_dotplot(adata, groupby = groupby)
    return adata

def anno(adata):
    # 采用AUCell方法进行标注
    # 对每个细胞的基因表达进行排序并且提取前5%
    # 找到marker并对重要性进行排序
    # 计算AUC
    # 计算每个类的平均AUC并展示
    pass

# 空转数据
# -表达矩阵
# --原始in_tissue细胞数
# --2000个与空转交集的高变基因
# -低维表示
# --PCA前两维表示
# --UMAP2维和3维表示
# --TSNE2维表示
# --空间位置2维表示
# -细胞类型解卷积
# --CARD解卷积
# --cell2location解卷积
# --空间区域marker基因
# --KEGG结果
# --GSEA结果
# --GO结果
# -邻接矩阵
# --k近邻矩阵

# 联合嵌入数据
# 单细胞和多个空转样本联合嵌入UMAP2维、3维表示(降采样到5k)
# 联合嵌入邻接矩阵（降采样到5k）两次标签传播算法

