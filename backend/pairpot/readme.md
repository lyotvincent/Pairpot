# Usage of Pairpot Library
`pairpot` library provides some powerful tools to analyze single-cell data and spatial transcriptomics data.

`pairpotlpa` library provides support of `lassoView` function of `pairpot`.

Run in console to install `pairpot` and `pairpotlpa`.

```bash
pip install pairpot
```

And use this code to import.
```python
import pairpot as pt
```

You can also find our library code in https://github.com/lyotvincent/Pairpot/tree/master/backend/pairpot.
<!-- import pairpotlpa as LPA -->
## Download Dataset from Pairpot

Corresponding adata can be downloaded in MetaInfo or use the code as following.

`pt.download` acceptes 3 parameters: `dataset_id`, `type` and `file`.`dataset_id` acceptes a string denoting the id of dataset in pairpot, `type` denotes the data type (`'sc'`: single-cell dataand `'sp'`: spatial transcriptomics data) and `file` denotes which file(`'complete'`: the complete file and `'meta'`: the processed file) needs to be download.

`pt.download` returns `adata` with anndata format. Following code shows the usage of `pt.download`
```python
# Use complete dataset_id to download
adata = pt.download("STDS0000235", type='sc', file='complete')

# Also, you can use the last 3 digits instead.
adata = pt.download("235", type='sc', file='complete')
```

## Perform Lasso-View Using Python Offline
`pt.lassoView` can perform Lasso-View for anndata offline.

`pt.lassoView` acceptes parameters `adata`, `selected`, `format`and `do_correct`.`adata` denotes the anndata format data or the path of h5ad file to be analyzed. `selected` acceptes a list to denote the lassoed cells indexs, which can also be gotten by lasso tools on the pairpot website. `format` denotes the format of `adata`, `'h5adfile'` denotes that parameter `adata` is a path of h5ad file and `'anndata'` denotes that parameter `adata` is an anndata format data.  `do_correct` denotes that after the label-propagation, does Lasso-View conduct `knn` rectification or not. In default, `format` is set to `'anndata'` and `do_correct` is set to `True`. `pt.lassoView` return a list of cell indexs, including ones are considered as the same cell type of the lassoed parts.

Following code shows the usage of `pt.lassoView`. 
```python
# Perform Pair-View using Python
# Run in console: pip install pairpot
import pairpot as pt

# Corresponding adata can be downloaded in MetaInfo or use the code as following.
adata = ad.read_h5ad("sample.h5ad")

lassoed_index = # Which can be generated in the pairpot website.
[70,153,166,290,343,423,513,537,
863,957,963,1111,1207,1252,1308,
1337,1354,1547,1750,1908,1928,1962,
1984,2043,2178,2211,2305,2343,2346,
2350,2427,2866,2932,2969]

refined_index = pt.lassoView(lassoed_index, adata)
```

## Perform Pair-View Using Python Offline
`pt.pairView` can perform PairView for anndata offline.

Similar to `pt.lassoView`, `pt.pairView` acceptes parameters `selected`, `scdata`, `spdata`, `format`. 

`selected` acceptes a list to denote the selected cells indexs, which can also be gotten on the pairpot website. `scdata` and `spdata` acceptes the anndata format single-cell and spatial transcriptomic data or the path of sc and sp h5ad file. `format` controls the format of `scdata`, and `spdata`: `'h5adfile'` denotes that parameter `scdata` and `spdata` are both path of h5ad files and `'anndata'` denotes that parameter `scdata` and `spdata` are both anndata format datas.

In default, `format` is set to `'anndata'`. `pt.pairView` return a list of proportion o f selected cells in each spot.

Following code shows the usage of `pt.pairView`. 
```python
# Perform Pair-View using Python
# Run in console: pip install pairpot
import pairpot as pt

# Corresponding adata can be downloaded in MetaInfo or use the code as following.
scdata = pt.download("STDS0000235", type='sc', file='complete')
spdata = pt.download("STDS0000235", type='sp', file='complete')
    
# Code for generating cell proportions for "your annotation" in SRT data offline
lassoed_index = 
[993,1991,2136,2685,2731,70,112,
153,166,290,343,423,513,537,863,
957,963,1111,1132,1207,1252,1308,
1337,1354,1547,1750,1908,1928,
1962,1984,2043,2131,2178,2211,
2292,2305,2343,2346,2350,2382,
2427,2486,2715,2866,2932,2969,2992]
props = pt.pairView(lassoed_index, scdata, spdata)

# the cell proportions inferred online can be found through 'Save' -> 'JSON'.
spdata.obs["your_annotation"] = props
```
## Using Online Refined Results Directly
As follow, pairpot website also provides refined_index in Lasso-View tool. which can be used directly to conducted further research.

```python
refined_index = 
[42,46,70,112,153,166,290,343,405,
423,513,518,537,710,863,942,957,
963,1085,1111,1132,1207,1252,1308,
1337,1354,1547,1703,1712,1750,1908,
1928,1962,1984,2043,2065,2131,2178,
2211,2292,2305,2343,2346,2350,2382,
2427,2715,2747,2866,2932,2944,2969,2992]

adata.obs['annotation'] = list(adata.obs['annotation'])
adata.obs['annotation'].iloc[refined_index] = "annotation 1"
adata.obs['annotation'] = adata.obs['annotation'].astype("category")
```
## Offline Data Processing
Pairpot library also provides some functions to curate data before conducting analysis.
### Preprocessing, Data Integration and Neighbor Graph Construction
`pt.lassoProc` function conducts data preprocess as follows:
- Removing minority of cells 
- Performing PCA
- Establishes similarity metrix $W$ and probability transition matrix $P$
- Perform dimensional reduction for visulization (using UMAP and t-SNE based on the neighbor graphs)
- Perform multi-scale celluar context representations of SRT data by MENDER

`pt.lassoProc` accepts only one parameter `path`, denotes the directory path. Then `lassoProc` process all h5 or h5ad file in the directory and save the processed file as `New_xxx.h5ad`.

Following code shows the usage of `pt.pairView`. 
```python
pt.lassoProc('./download')
```

### UCell Scores Evaluation
`pt.pairProc` function conducts offline clustering as follows:
- Segmenting the spatial domains in multiple slices using MENDER based on the multi-scale celluar context representations
- For scRNA-seq data, clusting cells using Leiden algorithm with the default resolution
- Using UCell to assign signature scores and evaluate cell types for datasets lacking original annotations

`pt.pairProc` accepts parameters: `adata`, `organs`, `top`, `alpha`, `n_jobs` and `clu_key`. `adata` denotes an anndata format data to be processed, `organs` denotes the organs and species from existing cell marker databases(like _PanglaoDB_ and _CellMarker_ ), `top` denotes the topK ratio, `alpha` denotes the resolution level, `n_jobs` denotes the number of threads for parallel processing and `clu_key` denotes the cluster approach.

In default, `top`, `alpha` and `n_jobs` are separatedly set to 0.05, 10e-40 and 16, and `clu_key` are set to  'leiden-1'.

Following code shows the usage of `pt.pairProc`. 
```python
adata=ad.read_h5ad('sample.h5ad')
organs=['Brain','Blood']
proc_adata=pt.pairProc(adata,organs)
```