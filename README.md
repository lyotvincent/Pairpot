# Pairpot: A Database with Real-time Lasso-based Analysis Tailored for Paired Single-cell and Spatial Transcriptomics

Paired single-cell and spatial resolved transcriptomics (SRT) data supplement each other, providing in-depth insights into biological processes and disease mechanisms.
Previous SRT databases have limitations in curating sufficient single-cell and SRT pairs (SC-SP pairs) and providing real-time heuristic analysis, which hinder the effort to uncover potential biological insights.
Here, we developed Pairpot (http://pairpot.bioxai.cn), a database tailored for paired single-cell and SRT data with real-time heuristic analysis. 
Pairpot curates 99 high-quality pairs including 1,425,656 spots from 299 datasets, and creates the association networks.
It constructs the curated pairs by integrating multiple slices and establishing potential associations between single-cell and SRT data. 
On this basis, Pairpot adopts semi-supervised learning that enables  real-time heuristic analysis for SC-SP pairs where Lasso-View refines the user-selected SRT domains within milliseconds, Pair-View infers cell proportions of spots based on user-selected cell types in real-time, and Layer-View displays SRT slices using a 3D hierarchical layout. 
Experiments demonstrated Pairpot's efficiency in identifying heterogeneous domains and cell proportions. 

<img src="./cover figure.png" width="60%">

## Offline Deployment
### FrontEnd

This project was bootstrapped with [Create React App](https://github.com/facebook/create-react-app).


Built using:

- Front-end library: React Ant-design Echarts
- CSS animations library: Animate.css

In the /frontend, you can run:

#### `npm start`

Runs the front-end in the development mode.\
Open [http://localhost:6634](http://localhost:6634) to view it in your browser.

The page will reload when you make changes.\
You may also see any lint errors in the console.

#### `npm test`

Launches the test runner in the interactive watch mode.\
See the section about [running tests](https://facebook.github.io/create-react-app/docs/running-tests) for more information.

#### `npm run build`

Builds the app for production to the `build` folder.\
It correctly bundles React in production mode and optimizes the build for the best performance.

The build is minified and the filenames include the hashes.\
Your app is ready to be deployed!

See the section about [deployment](https://facebook.github.io/create-react-app/docs/deployment) for more information.

## Back-End

This back-end was powered by Flask, Sqlite3, and neo4j via python. For environments installation:
```
    pip install -r requirements.txt # install pagkages
    python setup.py clean --all
    python setup.py build_ext --inplace # install LassoView server
```

The, in the ./backend you can run:

#### `python main.py`

Runs the back-end in the development mode.\
Open [http://localhost:5522](http://localhost:5522) to access the api for the browser.

## Dataset collection
Spatial transcriptomics and single-cell datasets were acquired and downloaded from databases such as the National Center for Biotechnology Information (NCBI), European Bioinformatics Institute (EBI), China National Center for Bio-information (CNCB), and 10X Genomics([https://www.10xgenomics.com](https://www.10xgenomics.com)).

Our strategy of dataset collection is below:

1. **SRT studies mentioned single-cell datasets**:
    - *SRT studies themselves provide single-cell data* : collect them directly
    - *SRT studies lack their own single-cell data but use previous single-cell datasets for integration analysis*:  collect these corresponding single-cell datasets instead
2. **SRT studies didn't mention single-cell datasets**: search for paired single-cell datasets according to the consistency of features such as species, tissues, and diseases among the SRT and single-cell studies.

If there exist multiple single-cell studies with the same features, Pairpot chooses the single-cell data that gives annotation files by default. When there exist multiple studies that provide annotation files, Pairpot chooses the study containing the widest variety of cell types by default.

## Data
The resource website for Pairpot is([http://src.bioxai.cn](http://src.bioxai.cn)). In the resource website, STpair.db contains all the paired metadata. Each folder for processed data is named with the last three id of dataset_ids (e.g., the processed data of STDS0000235 is in the folder 235.). There are also some additional files, including PanglaoDB, CellMarker2.0, and CellPhoneDB.

## Tools
We have developed several real-time heuristic analysis and data exploration tools for SC-SP pairs including Lasso-View, Layer-View, Pair-view, marker tables, interaction networks, and L-R pairs heatmap.

### Lasso-View
Users select cells (spots) of interest as a cluster, then Lasso-View would discover unselected cells similar to those in the user-defined cluster and remove the cells mistakenly selected into the cluster.

### Layer-View
Users can explore multiple slices of a study inthe left chart, or focus on a single slice in the right chart. In the 3D layout, the x-axis and y-axis denote the spatial coordinates of the slices, while the z-axis represents different batches of these slices. Users can rotate the axis to switch perspectives, click ’Inverse’ to hide all annotations, and highlight specific domains by clicking their legends.

### Pair-view
Users can use lasso tools to select cells of their interests in the single-cell chart, similar to Lasso-View. They can also use options in ’scConfigs’ to select customized cell types under different annotations and embeddings. After selecting customized cell types, users can click the ’Deconv’ buttons to call the Pair-View API. 


### More Alternative Tools

We added more alternative tools for deconvolution in ./server/Rsrc. Each of these following tools requires a single-cell and an SRT file with .h5ad format as input, and the results will append to the target file with .h5ad format. 

- Cell2location (by default in heat-net.py)
- RCTD (Rscript _RCTD.R -c your_single-cell_file -p your_SRT_file -o your_SRT_file)
- CARD (Rscript _CARD.R -c your_single-cell_file -p your_SRT_file -o your_SRT_file)
- SpaTalk (Rscript _SpaTalk-deconv.R -c your_single-cell_file -p your_SRT_file -o your_SRT_file)
- Seurat (Rscript _Seurat-deconv.R -c your_single-cell_file -p your_SRT_file -o your_SRT_file)
- CellTrek (Rscript _CellTrek.R -c your_single-cell_file -p your_SRT_file -o your_SRT_file)

We also added more alternative tools for inferring cell-cell interactions, so that users can perform their analytical tasks using different bioinformatic tools according to their preferences. 

- CellPhoneDB (by default in heat-net.py)
- CellChat (Rscript _RCTD.R -i your_h5ad_file -o your_h5ad_file -s your_species)
- iTALK ([https://github.com/Coolgenome/iTALK](https://github.com/Coolgenome/iTALK)) (Rscript _iTALK.R -i your_h5ad_file -o your_h5ad_file -s your_species)
- COMMOT (by default in ./server/heat-net.py)

## Environment
### Python Tools
*Python 3.10.12*

- scanpy 1.9.5
- pandas 2.0.0
- numpy 1.22.4
- anndata 0.9.1
- sklearn 1.2.2
- matplotlib 3.8.2
- sqlite3 3.37.2
- h5py 3.10.0
- torch 2.1.2
- torchvision 0.16.2
- joblib 1.3.2
- scipy 1.10.1

### R Tools
*R 4.1.2*
- Seurat 5.0.3
- getopt 1.20.4
- scater 1.22.0
- dplyr 1.1.4
- rhdf5 2.38.1
- viridis 0.6.5
- magrittr 2.0.3
- SeuratObject 5.0.1
- ConsensusClusterPlus 1.58.0
- dbscan 1.2.0
- akima 0.6.3.4
- randomForestSRC 3.2.0
- scales 1.3.0
- packcircles 0.3.6