import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

def layerView(adata, batch="batch", cluster="leiden-1",color_key='Pastel1'):
    spatial_data = adata.obsm['spatial']  # (n_cells, 2)
    z_values = adata.obs[batch].cat.codes # (n_cells,)

    if cluster in adata.obs:
        color_data = adata.obs[cluster]
    elif cluster in adata.obsm:
        color_data = adata.obsm[cluster]
    else:
        raise ValueError(f"{cluster} not found in adata.obs")

    palette = sns.color_palette(color_key, len(color_data.astype('category').cat.categories))
    colors = [palette[i] for i in color_data.astype('category').cat.codes]  # 使用颜色依据的编码

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(spatial_data[:, 0], spatial_data[:, 1], z_values, c=colors)
    ax.set_xlabel('X coordinate')
    ax.set_ylabel('Y coordinate')
    ax.set_zlabel('Batch category')

    return fig