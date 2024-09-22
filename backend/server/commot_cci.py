import commot as ct
import numpy as np
import pandas as pd

def scale_array(arr):
  min_val = 0
  max_val = np.max(arr)
  denominator = max_val - min_val
  scale_arr = (arr - min_val) / denominator
  return scale_arr

def LegRec_commot(adata, groupby='annotation'):
  df_cellchat = ct.pp.ligand_receptor_database(species='mouse',signaling_type =None, database='CellPhoneDB_v4.0')
  adata_dis500 = adata.copy()
  adata_dis500.obsm['spatial'] = adata_dis500.obsm['spatial'].astype('int')
  adata_dis500.var_names = [s.title() for s in adata_dis500.var_names]
  df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata_dis500, min_cell_pct=0.05)
  ct.tl.spatial_communication(adata_dis500, 
                              database_name='CellPhoneDB_v4.0', 
                              df_ligrec=df_cellchat_filtered, 
                              dis_thr=500, 
                              heteromeric=True, 
                              pathway_sum=True)
  
  ct_res = {}
  ct_res = {}
  ct_res['dataArray'] = {}
  ct_res['cellArray'] = {}
  ct_res['intArray'] = {}
  df_commot_all = pd.DataFrame()
  for i in df_cellchat_filtered.index:
      sender, receiver = df_cellchat_filtered[0][i], df_cellchat_filtered[1][i]
      ct.tl.cluster_communication(adata_dis500, 
                                  database_name='CellPhoneDB_v4.0', 
                                  clustering='annotation', 
                                  pathway_name=f'{sender}-{receiver}',
                                  n_permutations=100)
      sm_mat =np.array(adata_dis500.uns[f'commot_cluster-annotation-CellPhoneDB_v4.0-{sender}-{receiver}']['communication_matrix'])
      pval_mat =np.array(adata_dis500.uns[f'commot_cluster-annotation-CellPhoneDB_v4.0-{sender}-{receiver}']['communication_pvalue'])
      rowIdx = adata_dis500.uns[f'commot_cluster-annotation-CellPhoneDB_v4.0-{sender}-{receiver}']['communication_matrix'].index
      colIdx = adata_dis500.uns[f'commot_cluster-annotation-CellPhoneDB_v4.0-{sender}-{receiver}']['communication_matrix'].columns
      p_rows, p_cols = np.where(pval_mat < 0.5)
      p_indices = list(zip(p_rows, p_cols))
      s_rows, s_cols = np.where(sm_mat > 0.001)
      s_indices = list(zip(s_rows, s_cols))
      indices = list(set(s_indices).intersection(p_indices))
      df_commot = pd.DataFrame()
      df_commot["source"] = [rowIdx[r] for r,c in indices]
      df_commot["target"] = [colIdx[c] for r,c in indices]
      df_commot["interaction_group"] = [f'{sender}-{receiver}'.upper()]*len(indices)
      df_commot["celltype_group"] = [f"{rowIdx[r]}-{colIdx[c]}" for r, c in indices]
      df_commot["scores"] = np.array([sm_mat[r][c] for r, c in indices])
      df_commot["pvals"] = np.array([pval_mat[r][c] for r, c in indices])
      df_commot_all = pd.concat([df_commot_all, df_commot])
  df_commot_all = df_commot_all.reset_index(drop=True)
  df_commot_all['scaled_means'] = scale_array(df_commot_all["scores"])
  df_commot_all['neglog10p'] = -np.log10(df_commot_all['pvals']+1e-3)
  df_commot_all['significant'] = ['yes' if df_commot_all['pvals'][i] < 0.05 else 'nan' for i in df_commot_all.index ]
  df_commot_all['source'] = df_commot_all['source'].astype('str')
  df_commot_all['target'] = df_commot_all['target'].astype('str')
  cell_types = adata.obs[groupby].cat.categories
  _ct = list(cell_types)
  for temp_cell1 in cell_types:
      df_heat = df_commot_all.iloc[np.array(df_commot_all['source'] == temp_cell1) + 
                              np.array(df_commot_all['target'] == temp_cell1),: ]
      df_heat = df_heat[['interaction_group',
                              'celltype_group',
                              'scaled_means',
                              'pvals', 
                              'neglog10p',
                              'significant']]
      if len(df_heat) > 0:
          if min(df_heat['pvals']) < 0.05:
              ct_res['dataArray'][temp_cell1]=df_heat
              y = df_heat['interaction_group'].unique()
              ct_res['intArray'][temp_cell1] = y.T

              x = df_heat['celltype_group'].unique()
              ct_res['cellArray'][temp_cell1]=x.T
              print(temp_cell1,"heatmap if OK")
          else: 
              _ct.remove(temp_cell1)
              print(f"{temp_cell1} has no L-R pairs.")
      else:
          _ct.remove(temp_cell1)
          print(f"{temp_cell1} has no L-R pairs.")
  ct_res['cellType'] = np.array(_ct)
  return ct_res