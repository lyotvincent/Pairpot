import React from 'react';
import Markdown from 'react-markdown';
import ReactMarkdown from 'react-markdown';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { solarizedlight } from 'react-syntax-highlighter/dist/esm/styles/prism';

var Rcode = 
`
### Seurat
**Source:** https://github.com/satijalab/seurat\n
**Citation:** Butler, A., Hoffman, P., Smibert, P., Papalexi, E., and Satija, R. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat. Biotechnol. 36, 411-420.

### Key features and functions of SeuratDeconv include:
- Bulk RNA-seq Deconvolution
- Spatial Transcriptomics Integration

### Rcode
**********
\`\`\`r 
# fetching params
library(getopt)
spec <- matrix(c("help", 'h', 0, 'logical', 'help',
                 "scpath", "c", 1, "character", "path to single-cell h5ad file",
                 "sppath", "p", 1, "character", "path to spatial h5ad file",
                 "outpath", "o", 2, "character", "path to output h5ad file, default as sp path",
                 "raw", "r", 2, "logical", "whether use raw data in h5ad file"),byrow=T,ncol=5)

args <- getopt(spec)
if (!is.null(args$help) || is.null(args$scpath) || is.null(args$sppath) ) {
  cat(paste(getopt(spec, usage = T), "\\n"))
  q(status=1)
}
outpath <- args$outpath
if(is.null(outpath)){
  outpath <- args$sppath
}
raw <- args$raw
if(is.null(raw)){
  raw <- TRUE
}

source("_loading.R")
library(spatstat.explore)
library(Seurat)
library(scater)
library(SeuratObject)
library(Matrix)
library(Rcpp)
library(rhdf5)

Preprocess_Seurat <- function(seu, assay = "Spatial"){
  seu <- tryCatch(expr = {
    seu <- SCTransform(seu, assay = assay, return.only.var.genes = FALSE, verbose = FALSE)
    seu <- RunPCA(seu, assay = "SCT", verbose = FALSE, features = VariableFeatures(seu))
    return(seu)
  }, error = function(e){
    print(e)
    print("Using LogNorm, FindVariableFeatures, Scale Pipline.")
    seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    seu <- FindVariableFeatures(seu, selection.method = "vst")
    seu <- ScaleData(seu)
    seu <- RunPCA(seu, assay = assay, verbose = FALSE, features = VariableFeatures(seu))
    ppAssay <- assay
    return (seu)
  })
  return (seu)
}

Cluster_Seurat <- function(seu){
  seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
  seu <- FindClusters(seu, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:30)
  return(seu)
}

Estimation_Seurat <- function(seu){
   de_markers <-FindAllMarkers(seu, assay = "SCT")
   return(de_markers)
}

Deconvolution_Seurat <- function(seu.sp, seu.sc){
  anchors <- tryCatch({
    FindTransferAnchors(reference = seu.sc, query = seu.sp, normalization.method = "SCT")
  }, error = function(e){
    print(e)
    anchors <- FindTransferAnchors(reference = seu.sc, query = seu.sp, normalization.method = "LogNormalize")
    return(anchors)
  })
  predictions.assay <- TransferData(anchorset = anchors, refdata = seu.sc$annotation, prediction.assay = TRUE,
                                    weight.reduction = seu.sp[["pca"]], dims = 1:30)
  return(predictions.assay)
}

FindMarkers_Seurat <- function(seu, ident){
  de_markers <- FindMarkers(seu, ident.1 = ident)
  return(de_markers)
}

Run_Seurat <- function(scpath, sppath, outpath, raw){
  sce <- Load_h5adsc_to_SCE(scpath, raw = raw)
  spe <- Load_h5adsp_to_SCE(sppath, raw = raw)
  rownames(sce) <- toupper(rownames(sce))
  rownames(spe) <- toupper(rownames(spe))
  rownames(sce@assays@data@listData$counts) <- toupper(rownames(sce@assays@data@listData$counts))
  rownames(spe@assays@data@listData$counts) <- toupper(rownames(spe@assays@data@listData$counts))
  spe <- logNormCounts(spe)
  colnames(spe) <- Make_unique(colnames(spe))
  spseu <- SeuratObject::as.Seurat(x = spe, counts = "counts", data = "logcounts")
  sce <- logNormCounts(sce)
  colnames(sce) <- Make_unique(colnames(sce))
  scseu <- SeuratObject::as.Seurat(x = sce, counts = "counts", data = "logcounts")

  spseu@images$image <- new(
    Class = "SlideSeq",
    assay = "originalexp",
    key = "image_",
    coordinates = data.frame(
      spatial_X = spe$spatial_X,
      spatial_Y = spe$spatial_Y,
      row.names = colnames(spe)
    )
  )

  scseu <- Preprocess_Seurat(scseu, assay = "originalexp")
  spseu <- Preprocess_Seurat(spseu, assay = "originalexp")
  mat <- Deconvolution_Seurat(spseu, scseu)
  results <- as.matrix(mat@data[-dim(mat@data)[1],])
  results <- results[sort(rownames(results)),]
  results <- results[,colnames(spe)]
  h5write(results, outpath, '/obsm/Seurat')
  message(paste0("Deconvolution with 'Seurat' finished, idents saved at /obsm/Seurat in", outpath))
}

cat("#### Running Seurat-deconv ####\\n")
cat(paste0("Single-cell h5ad file:", args$scpath, "\\n"))
cat(paste0("Spatial h5ad file:", args$sppath, "\\n"))
cat(paste0("Output h5ad file:", outpath, "\\n"))
cat(paste0("Using raw:", as.character(raw), "\\n"))
Run_Seurat(args$scpath, args$sppath, outpath, raw)
cat("#### Seurat-deconv Finished ####\\n")
\`\`\`
`

const Seurat = () => {
    return (
        <div style={{height:'25rem',overflow: 'auto'}}>
            {/* <ReactMarkdown>{Rcode}</ReactMarkdown> */}
            <ReactMarkdown
                children={Rcode}
                components={{
                code({ node, inline, className, children, ...props }) {
                    return inline ? (
                    <code className={className} {...props}>
                        {children}
                    </code>
                    ) : (
                    <SyntaxHighlighter
                        customStyle={{ backgroundColor: 'white' }}
                        style={solarizedlight}
                        language="r"
                        PreTag="div"
                        {...props}
                    >
                        {String(children)}
                    </SyntaxHighlighter>
                    );
                },
                }}
            />
        </div>
  );
};

export default Seurat;