import scanpy as sc
import numpy as np
from sklearn.linear_model import LinearRegression

adata_sc = sc.read_h5ad("sc1_sampled.h5ad")
adata_sp = sc.read_h5ad("sp1-deconv.h5ad")

# generate A, 用UCell指标，感觉效果很好
ucells = [s for s in adata_sc.obs.columns if s.startswith('UCell_')]
ucells.append('annotation')
df_sc = adata_sc.obs[ucells]
A = df_sc.groupby(by='annotation').mean().T

# generate b
df_sp = adata_sp.obs[[s for s in adata_sc.obs.columns if s.startswith('UCell_')]]
B = df_sp.T

# NNLS model
reg_nnls = LinearRegression(positive=True)  # 此处使用
pred_nnls = reg_nnls.fit(A, B)

# show results，Interneurons是ground truth，Interneurons1是NNLS得到的，看起来再做个knn平滑会更好
adata_sp.obs['Interneurons1'] = pred_nnls.coef_[:, 2]
sc.pl.spatial(adata_sp[adata_sp.obs['batch']=="Control02"], color=['Interneurons', 'Interneurons1'], spot_size=70)