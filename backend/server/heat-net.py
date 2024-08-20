from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import simplejson
import ktplotspy as kpy
import networkx as nx
from matplotlib import cm, colors
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import anndata as ad
import numpy as np
from normalize import Cell2Location_run

def get_network(ax):
    # 节点+边转格式
    network_data = {}
    network_data["nodes"] = []
    network_data["links"] = []
    network_data["categories"] = []

    # 遍历节点
    maxValue = max(ax['interaction_count'].iloc[:,0])
    minValue = max(ax['interaction_count'].iloc[:,0])
    frac = maxValue - minValue
    for i in range(len(ax['interaction_count'])):
        nodes_dic = {}
        name = ax['interaction_count'].iloc[i].name
        value = ax['interaction_count'].iloc[i][0]
        nodes_dic['name'] = name
        nodes_dic['id'] = name # link必须是(id,id)对应
        nodes_dic['value'] = value
        nodes_dic['category'] = name
        # if value == 0:# 孤立节点固定
        #     nodes_dic['fixed'] = True
        #     nodes_dic['x'] = 60
        #     nodes_dic['y'] = 100
        # 
         
        # # minmax norm for symbolSize
        # if frac == 0:
        #     nodes_dic['symbolSize'] = 20
        # else:
        #     nodes_dic['symbolSize'] = 20*((value - minValue)/frac +1)

        # test symbolSize
        # nodes_dic['symbolSize'] = 100
            
        # if value <20:
        #     nodes_dic['symbolSize'] = 20
        # elif value >= 100:
        #     nodes_dic['symbolSize'] = value/3
        # else :
        #     nodes_dic['symbolSize'] = value/2
        network_data["nodes"].append(nodes_dic)
        # 不能直接在前面append，指针原因
        cat_dic = {}
        cat_dic['name'] = name
        network_data['categories'].append(cat_dic)


    # 遍历边
    for i in range(len(ax['interaction_edges'])):
        # 为0的边不用添！！！
        count = ax['interaction_edges'].iloc[i]['COUNT']
        # if count <= 8:
        #     continue

        edges_dic = {}
        source = ax['interaction_edges'].iloc[i]['SOURCE']
        target = ax['interaction_edges'].iloc[i]['TARGET']
        edges_dic['source'] = source
        edges_dic['target'] = target
        edges_dic['value'] = count
        # 自连自己的边要处理下
        # if source == target:
        #     style_dic = {}
        #     style_dic ['curveness'] = 0.5
        #     edges_dic['lineStyle'] = style_dic
        network_data['links'].append(edges_dic)

    return network_data

def add_heatmap(cpdb_res: dict, adata_in: ad.AnnData, groupby: str):   
    cpdb_res['dataArray'] = {}
    cpdb_res['cellArray'] = {}
    cpdb_res['intArray'] = {}

    cell_types = adata_in.obs[groupby].cat.categories
    _ct = list(cell_types)
    for temp_cell1 in cell_types:
        df_heat = kpy.plot_cpdb(
            adata=adata_in,
            cell_type1= temp_cell1,
            cell_type2='.',
            means=cpdb_res['means'],
            pvals=cpdb_res['pvalues'],
            celltype_key=groupby,
            keep_significant_only=True # 其实默认都是true的，保留p<0.05的
        )
        df_heat = df_heat.data
        if min(df_heat['pvals']) < 0.05:
            # filter df_heat for error included celltypes
            true_ct1 = [f'{s}-{temp_cell1}' for s in cell_types]
            true_ct2 = [f'{temp_cell1}-{s}' for s in cell_types]
            true_ct = true_ct1 + true_ct2
            df_heat = df_heat[df_heat['celltype_group'].isin(true_ct)]
            df_heat = df_heat[['interaction_group','celltype_group','scaled_means','pvals', 'neglog10p','significant']]
            if groupby == 'leiden-1': # strict for sp data
                df_heat = df_heat[(df_heat['pvals'] < 0.5) |(df_heat['significant'] == 'yes')]
            df_heat = df_heat.astype(str).fillna("")
            cpdb_res['dataArray'][temp_cell1]=df_heat

            y = df_heat['interaction_group'].unique()
            cpdb_res['intArray'][temp_cell1] = y.T

            x = df_heat['celltype_group'].unique()
            cpdb_res['cellArray'][temp_cell1]=x.T
            print(temp_cell1,"heatmap if OK")
        else: 
            _ct.remove(temp_cell1)
            print(f"{temp_cell1} has no L-R pairs.")
        cpdb_res['cellType'] = _ct

    return cpdb_res
def write_heat(adata_out,cpdb_res):
    # del adata_sp.uns['dataArray']
    adata_out.uns['dataArray'] = cpdb_res['dataArray']
    adata_out.uns['intArray'] = cpdb_res['intArray']
    adata_out.uns['cellArray'] = cpdb_res['cellArray']
    adata_out.uns['cellType'] = cpdb_res['cellType']
    adata_out.uns['dataKeys'] = [
    "interaction_group",
    "celltype_group",
    "scaled_means",
    "pvals",
    "neglog10p",
    "significant",
    ]


def convert_to_string(data):
    if isinstance(data, dict):
        return {key: convert_to_string(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [convert_to_string(item) for item in data]
    elif not isinstance(data, str):
        return str(data)
    else:
        return data

def write_net(network_data, adata_out):
    # 处理 NA 值并将数据存储到 network_data1 字典中
    network_data1 = convert_to_string(network_data)
    network_data2 = {}
    for key in network_data1:
        network_data2[key] = pd.DataFrame(network_data1[key]).fillna("")
        # print(network_data[key])
    # network_data2 = pd.DataFrame(network_data2).fillna("")
    # network_data2['categories']
    adata_out.uns['network'] = network_data2

def save_delete(adata, loc, column):
  if loc == 'var' and column in adata.var.columns:
    del adata.var[column]
  if loc == 'obs' and column in adata.obs.columns:
    del adata.obs[column]
  if loc == 'obsm' and column in adata.obsm.keys():
    del adata.obsm[column]
  if loc == 'varm' and column in adata.varm.keys():
    del adata.varm[column]
  print(f"deleted adta.{loc}['{column}']")

def run_cpdb(adata_inFile, adata_outPath, type='sc'):
    if type == 'sp':
        groupby = 'leiden-1'  # adata.obs['leiden-1'] for sp data
    else:
        groupby = 'annotation'  # adata.obs['annotation'] for sc data
    adata_in = sc.read_h5ad(adata_inFile) 
    # make obs unique, or error occurs in cpdb
    adata_in.obs_names_make_unique()
    adata_in.var_names= [s.upper() for s in adata_in.var_names]
    adata_in = adata_in[:, adata_in.var['highly_variable']].copy() # use highly variable genes
    adata_in.write_h5ad(f"{adata_outPath}/{type}_sampled.h5ad")
    meta = pd.DataFrame(adata_in.obs[groupby].copy())
    meta = meta.reset_index()
    meta.columns=['Cell', 'cell_type']
    meta.to_csv("../resources/meta.txt", sep='\t', index=False)

    # run cpdb
    cpdb_res = cpdb_statistical_analysis_method.call(
    cpdb_file_path="../resources/cellphonedb.zip",
    meta_file_path="../resources/meta.txt",
    counts_file_path=f"{adata_outPath}/{type}_sampled.h5ad",
    counts_data='hgnc_symbol',
    output_path="../resources/cpdb",
    threads=16,
    )
    print("calculating cpdb network...")
    ax = kpy.plot_cpdb_heatmap(pvals=cpdb_res['pvalues'],return_tables=True)
    network_data = get_network(ax)
    cpdb_res['network'] = network_data
    print("calculating cpdb heatmap...")
    add_heatmap(cpdb_res, adata_in=adata_in, groupby=groupby)
    # generate adata_out1(sampled, with X, stored in database, only for download)
    if type == 'sc':  # downsample for single-cell data
        sampleIdx = np.random.choice(len(adata_in), 3000)
        adata_out1 = adata_in[sampleIdx,adata_in.var['highly_variable']].copy()
    else:
        adata_out1 = adata_in[:,adata_in.var['highly_variable']].copy()
    adata_out1.obs.columns = [s.replace('/', ' or ') for s in adata_in.obs.columns]
    adata_out1.uns['CellPhoneDB_pvalues'] = cpdb_res['pvalues'].fillna("")
    adata_out1.uns['CellPhoneDB_means'] =  cpdb_res['means'].fillna("")
    sc.tl.rank_genes_groups(adata_out1, groupby=groupby,method='wilcoxon')
    sc.tl.dendrogram(adata_out1, groupby=groupby)
    if type == 'sc': 
        adata_out1.write_h5ad(f"{adata_outPath}/{type}_sampled.h5ad")
        print(f"Info:: Write adata_out1 to {adata_outPath}/{type}_sampled.h5ad")

    # generate adata_out2(meta, without X, interface data with frontend)
    adata_out2 = ad.AnnData(obs=adata_out1.obs, obsm=adata_out1.obsm, var=adata_out1.var)
    adata_out2.uns['dendrogram'] = adata_out1.uns[f'dendrogram_{groupby}']
    cell_types = adata_out1.obs[groupby].cat.categories
    nameRank = pd.DataFrame(adata_out1.uns['rank_genes_groups']['names'], columns=cell_types)
    expr_df = pd.DataFrame(columns=cell_types)
    for ct in cell_types:
        adata_cp = adata_out1[:, list(nameRank[ct])]
        expr = adata_cp[adata_cp.obs[groupby]==ct].X.mean(axis=0).A1
        expr_df[ct] = expr

    frac_df = pd.DataFrame(columns=cell_types)
    for ct in cell_types:
        adata_cp = adata_out1[:, list(nameRank[ct])]
        frac = np.array(np.mean(adata_cp[adata_cp.obs[groupby]==ct].X > 1, axis=0)).flatten()
        frac_df[ct] = frac

    adata_out2.uns['rank_genes_groups'] = adata_out1.uns['rank_genes_groups']
    del adata_out2.uns['rank_genes_groups']['pvals']
    del adata_out2.uns['rank_genes_groups']['scores']
    adata_out2.uns['rank_genes_groups']['names'] = np.array(pd.DataFrame(adata_out1.uns['rank_genes_groups']['names'])).T
    adata_out2.uns['rank_genes_groups']['expr'] = np.array(expr_df).T
    adata_out2.uns['rank_genes_groups']['pvals_adj'] = np.array(pd.DataFrame(adata_out1.uns['rank_genes_groups']['pvals_adj'])).T
    adata_out2.uns['rank_genes_groups']['frac'] = np.array(frac_df).T
    adata_out2.uns['rank_genes_groups']['logfoldchanges'] = np.array(pd.DataFrame(adata_out1.uns['rank_genes_groups']['logfoldchanges'])).T
    save_delete(adata_out2, 'obsm', 'AUCell_rankings')
    save_delete(adata_out2, 'obsm', 'X_pca')
    save_delete(adata_out2, 'obsm', 'X_pca_harmony')
    if type=='sp':
        adata_out2.obs['annotation'] = adata_out2.obs['leiden-1']
        save_delete(adata_out2, 'obs', 'in_tissue')
        save_delete(adata_out2, 'var', 'feature_types')
        save_delete(adata_out2, 'var', 'n_cells')
        save_delete(adata_out2, 'var', 'means')
        save_delete(adata_out2, 'var', 'gene_name')
        save_delete(adata_out2, 'var', 'SYMBOL')
        Ucells = adata_out2.obs.columns[[s.startswith("UCell") for s in adata_out2.obs.columns]]
        for uc in Ucells: # del Ucells for meta data
            save_delete(adata_out2, 'obs', uc)
    save_delete(adata_out2, 'var', 'genome')
    save_delete(adata_out2, 'var', 'highly_variable')
    save_delete(adata_out2, 'var', 'dispersions')
    save_delete(adata_out2, 'var', 'dispersions_norm')
    write_heat(adata_out2,cpdb_res)
    write_net(network_data,adata_out2)
    adata_out2.write_h5ad(f"{adata_outPath}/{type}_meta.h5ad")
    print(f"Info:: Write adata_out2 to {adata_outPath}/{type}_meta.h5ad")

def run_deconv(adata_scFile, adata_spFile, adata_outPath):
    adata_sc = sc.read_h5ad(adata_scFile)
    adata_sp = sc.read_h5ad(adata_spFile)
    adata_sc.var_names = [s.upper() for s in adata_sc.var_names]
    adata_sp.var_names = [s.upper() for s in adata_sp.var_names]
    adata_sc = adata_sc[:, adata_sc.var['highly_variable']].copy()
    adata_sp = adata_sp[:, adata_sp.var['highly_variable']].copy()
    weight = Cell2Location_run(adata_sc, adata_sp)
    adata_sp.obs = pd.concat([adata_sp.obs, weight], axis=1)
    adata_sp = adata_sp[:, adata_sp.var['highly_variable']].copy()
    adata_sp.write_h5ad(f"{adata_outPath}/sp_deconv.h5ad")
    print(f"Write dcv results to {adata_outPath}/sp_deconv.h5ad")

if __name__ == '__main__':
    adata_spFile = "/data/rzh/RawUrls/153/STDS0000153/New_153.h5ad"
    adata_scFile = "/data/rzh/RawUrls/153/SCDS0000153/New_153.h5ad"
    adata_outPath = "/data/rzh/RawUrls/153"
    run_deconv(adata_scFile, adata_spFile, adata_outPath)
    run_cpdb(adata_scFile, adata_outPath, type="sc")
    run_cpdb(f"{adata_outPath}/sp_deconv.h5ad", adata_outPath, type="sp")