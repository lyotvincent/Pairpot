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
# --原始数据自带标注
# --细胞类型marker基因
# --KEGG结果
# --GSEA结果
# --GO结果
# -邻接矩阵
# --k近邻矩阵
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from .AUCell import *
import MENDER
from .db import panglaoDB
import cell2location
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

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
    adata.var_names = [s.upper() for s in adata.var_names]
    adata.obs_names = [s.upper() for s in adata.obs_names]
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata

def input_adata_h5ad(sample):
#   print(sample)
  adata = sc.read_h5ad(sample)
  adata.var_names = [s.upper() for s in adata.var_names]
  adata.obs_names = [s.upper() for s in adata.obs_names]
  adata.obs_names_make_unique()
  adata.var_names_make_unique()
  return adata

def input_adata_txt(sample):
    adata = sc.read_text(sample)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata

def input_adata_csv(sample):
    adata = sc.read_text(sample, delimiter=',')
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    return adata

# 将所有矩阵合并
def concat_adata(samples, sampleNames, inputFunc=input_adata_10Xh5):
    adatas = []
    for i in range(len(sampleNames)):
        adata = inputFunc(samples[i])
        sc.pp.filter_cells(adata, min_genes=5)
        sc.pp.filter_cells(adata, min_counts=30)
        adatas.append(adata)
    if len(adatas) > 0:
        intersection_var = set(adatas[0].var_names)
        for a in adatas[1:]:
            intersection_var &= set(a.var_names)
        common_vars = list(intersection_var)
        for i in range(len(adatas)):
            adatas[i] = adatas[i][:, common_vars]
    # 进行数据合并
    adata_concat = adatas[0].concatenate(adatas[1:], batch_categories=sampleNames)
    return adata_concat


# 预处理
def pp(adata:ad.AnnData):
    mito_genes = adata.var_names.str.startswith('MT')
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    mt_frac = np.sum(
        adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
    adata.obs['mt_frac'] = np.array(mt_frac)
    # 过滤低表达的基因
    sc.pp.filter_cells(adata, min_genes=5)  # 过滤一个细胞中表达少于五个基因的细胞样本 
    sc.pp.filter_genes(adata, min_cells=5)  # 过滤在少于五个细胞中表达的基因
    sc.pp.filter_cells(adata, min_counts=30)   # 过滤每个细胞中计数少于29个的细胞样本 

    # 过滤线粒体核糖体基因
    rp_genes = adata.var_names.str.startswith('RP')
    mt_genes = adata.var_names.str.startswith('MT')
    adata = adata[:, ~(rp_genes + mt_genes)]
    adata = adata[adata.obs['mt_frac'] < 0.2]
    adata.layers['Raw'] = adata.X
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,min_disp=0.5, n_top_genes=2000)
    return adata

def clu(adata, key_added="leiden-1", n_neighbors=50, n_pcs=30, rep='X_pca_harmony', do_har=True, max_iter=20, resolution=1, do_scrublet=True, har_key='batch'):
    # Computing the neighborhood graph
    if do_scrublet:
        n0 = adata.shape[0]
        print("{0} Cell number: {1}".format(key_added, n0))
        sc.external.pp.scrublet(adata)
        adata = adata[adata.obs['predicted_doublet']==False,:].copy()
        print("{0} Cells retained after scrublet, {1} cells reomved.".format(adata.shape[0], n0-adata.shape[0]))
    else:
        print("Ignoring processing doublet cells...")
    sc.pp.pca(adata, svd_solver='arpack', use_highly_variable=True)
    if do_har and len(adata.obs[har_key].cat.categories) > 1:
        sc.external.pp.harmony_integrate(adata, key=har_key,max_iter_harmony=max_iter)
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    # Run UMAP
    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    sc.tl.leiden(adata, key_added=key_added, resolution=resolution)
    sc.pl.umap(adata, color=key_added, legend_fontoutline=True, palette=sc.pl.palettes.default_20, legend_loc="on data")
    return adata

def mender(adata:ad.AnnData, key_added="mender"):
    # input parameters of MENDER
    scale = 6
    radius = 15
    # main body of MENDER
    msm = MENDER.MENDER(
        adata,
        # determine which cell state to use
        # we use the cell state got by Leiden
        batch_obs='batch',
        ct_obs='leiden-1'
    )
    msm.prepare()
    # set the MENDER parameters

    msm.set_MENDER_para(
    # default of n_scales is 6
    n_scales=scale,

    # for single cell data, nn_mode is set to 'radius'
    nn_mode='radius',

    # default of n_scales is 15 um (see the manuscript for the analysis).
    # MENDER also provide a function 'estimate_radius' for estimating the radius
    nn_para=radius,
    )

    # construct the context representation
    msm.run_representation_mp()

    msm.run_clustering_normal(-0.5)
    adata.obs[key_added] = np.array(msm.adata_MENDER.obs['MENDER'])
    adata.obs['annotation'] = np.array(msm.adata_MENDER.obs['MENDER'])
    return adata

def rank(adata, organs, 
         method="AUCell", 
         top=0.05, 
         alpha=10e-40, 
         n_jobs=16,
         clu_key='leiden-1',
         test_func=mannwhitneyu,
         ):
  # 对每个细胞的基因表达进行排序并且提取前5%
  adata = AUCell_buildRankings(adata, top=top)

  # Find the potential Celltypes
  celltype = pd.unique(panglaoDB[panglaoDB['organ'].isin(organs)]['cell type'].dropna())
  celltype.sort()

  # UCell_Assign
  adata = AUCell_UCAssign(adata, 
                          db=panglaoDB, 
                          celltype=celltype, 
                          alpha=alpha, 
                          n_jobs=n_jobs, 
                          clu_key=clu_key,
                          test_func=test_func)
  return adata

def anno(adata:ad.AnnData, annoDict:dict):
  adata.obs['annotation'] = 'Unknown'
  for key in annoDict.keys():
    adata.obs.loc[adata.obs['leiden-1'].isin(annoDict[key]), 'annotation'] = key
  return adata

def marker(adata, groupby="annotation", method='wilcoxon'):
    sc.tl.rank_genes_groups(adata, groupby = groupby, method = method)
    sc.tl.dendrogram(adata, groupby=groupby)
    sc.pl.rank_genes_groups_dotplot(adata, groupby = groupby)
    return adata

def Cell2Location_rg_sc(adata_sc, max_epoches=250, batch_size=2500, train_size=1, lr=0.002,
                        num_samples=1000, use_gpu=True):
    from cell2location.models import RegressionModel

    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_sc,
                            # cell type, covariate used for constructing signatures
                            labels_key='annotation')

    # create and train the regression model
    mod = RegressionModel(adata_sc)
    # mod.view_anndata_setup()

    # Use all data for training (validation not implemented yet, train_size=1)
    mod.train(max_epochs=max_epoches, batch_size=batch_size, train_size=train_size, lr=lr)

    # plot ELBO loss history during training, removing first 20 epochs from the plot
    # mod.plot_history(20)

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_sc = mod.export_posterior(
        adata_sc, sample_kwargs={'num_samples': num_samples, 'batch_size': batch_size, 'use_gpu': use_gpu}
    )

    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
        inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_sc.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_sc.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_sc.uns['mod']['factor_names']
    return inf_aver


def Cell2Location_rg_sp(adata_sp, inf_aver, N_cells_per_location=30, detection_alpha=20,
                        max_epoches=5000, batch_size=None, train_size=1, lr=0.002,
                        num_samples=1000, use_gpu=True):
    # do spatial mapping
    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_sp.var_names, inf_aver.index)
    adata_sp = adata_sp[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_sp, batch_key="batch")

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_sp, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=N_cells_per_location,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=detection_alpha
    )
    # mod.view_anndata_setup()
    mod.train(max_epochs=max_epoches,
              # train using full data (batch_size=None)
              batch_size=batch_size,
              # use all data points in training because
              # we need to estimate cell abundance at all locations
              train_size=train_size,
              use_gpu=use_gpu)

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    # mod.plot_history(1000)

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_sp = mod.export_posterior(
        adata_sp, sample_kwargs={'num_samples': num_samples, 'batch_size': mod.adata.n_obs, 'use_gpu': use_gpu}
    )

    # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
    # to adata.obs with nice names for plotting
    adata_sp.obs[adata_sp.uns['mod']['factor_names']] = adata_sp.obsm['q05_cell_abundance_w_sf']
    return adata_sp

def Cell2Location_run(adata_sc, adata_sp, sc_max_epoches=250, sc_batch_size=2500, sc_train_size=1, sc_lr=0.002, sc_num_samples=1000,
                      N_cells_per_location=30, detection_alpha=20,
                      sp_max_epoches=5000, sp_batch_size=2500, sp_train_size=1, sp_lr=0.002, sp_num_samples=1000, use_gpu=True):
    # rename genes to ENSEMBL
    adata_sc.var['SYMBOL'] = adata_sc.var.index
    adata_sp.var['SYMBOL'] = adata_sp.var_names
    adata_sp.var['gene_name'] = adata_sp.var_names
    cell_count_cutoff=5
    cell_percentage_cutoff2=0.03
    nonz_mean_cutoff=1.12
    adata_sc = adata_sc.copy()
    adata_sp = adata_sp.copy()
    adata_sc.X = adata_sc.layers['Raw'].astype('int')
    adata_sp.X = adata_sp.layers['Raw'].astype('int')
    print(f"Single-cell shape:{adata_sc.shape}")
    selected = cell2location.utils.filtering.filter_genes(adata_sc,
                        cell_count_cutoff=cell_count_cutoff,
                        cell_percentage_cutoff2=cell_percentage_cutoff2,
                        nonz_mean_cutoff=nonz_mean_cutoff)
    adata_sc = adata_sc[:, selected]
    inf_aver = Cell2Location_rg_sc(adata_sc, sc_max_epoches, sc_batch_size, sc_train_size,sc_lr,sc_num_samples, use_gpu)
    adata_sp = Cell2Location_rg_sp(adata_sp, inf_aver,
                                   N_cells_per_location=N_cells_per_location,
                                   detection_alpha=detection_alpha,
                                   max_epoches=sp_max_epoches,
                                   batch_size=sp_batch_size,
                                   train_size=sp_train_size,
                                   lr=sp_lr,
                                   num_samples=sp_num_samples,
                                   use_gpu=use_gpu)
    weight = adata_sp.obsm['q05_cell_abundance_w_sf']
    weight.columns = adata_sp.uns['mod']['factor_names']
    return weight

def CellphoneDB_run(adata_path, cpdb_file_path="resources/cellphonedb.zip", output_path="resources/cpdb", counts_data='hgnc_symbol'):
    adata = sc.read_h5ad(adata_path)
    meta = pd.DataFrame(adata.obs['annotation'].copy())
    meta = meta.reset_index()
    meta.columns=['Cell', 'cell_type']
    meta.to_csv("resources/meta.txt", sep='\t', index=False)
    cpdb_res = cpdb_statistical_analysis_method.call(
        cpdb_file_path=cpdb_file_path,
        meta_file_path="resources/meta.txt",
        counts_file_path=adata_path,
        counts_data=counts_data,
        output_path=output_path
    )
    return cpdb_res

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

