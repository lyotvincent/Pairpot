import React from 'react';
import Markdown from 'react-markdown';
import ReactMarkdown from 'react-markdown';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { solarizedlight } from 'react-syntax-highlighter/dist/esm/styles/prism';

var Rcode = 
`
### Card
Card (Cellular Annotation for RNA-seq Data) is a tool used for cell type annotation in single-cell RNA sequencing (scRNA-seq) data, utilizing reference datasets for accurate classification.

### Key features and functions of Card include:
- Cell Type Annotation
- Reference-based Annotation
### Rcode
**********
\`\`\`r 
# fetching params
library(getopt)
spec <- matrix(c(
  "help", "h", 0, "logical", "help",
  "scpath", "c", 1, "character", "path to single-cell h5ad file",
  "sppath", "p", 1, "character", "path to spatial h5ad file",
  "outpath", "o", 2, "character", "path to output h5ad file,
  default as sp path",
  "raw", "r", 2, "logical", "whether use raw data in h5ad file"
), byrow = TRUE, ncol = 5)

args <- getopt(spec)
if (!is.null(args$help) || is.null(args$scpath) || is.null(args$sppath)) {
  cat(paste(getopt(spec, usage = TRUE), "\\n"))
  q(status = 1)
}
outpath <- args$outpath
if (is.null(outpath)) {
  outpath <- args$sppath
}
raw <- args$raw
if (is.null(raw)) {
  raw <- TRUE
}

source("_loading.R")
Check_Load_BiocPackages("TOAST")
Check_Load_GithubPackages("MuSiC", "xuranw/MuSiC@7c58348")
Check_Load_GithubPackages("CARD", "YingMa0107/CARD")

# install TOAST, MuSiC


Deconvolution_CARD <- function(spe, sce) {
  rownames(sce@assays@data$counts) <- Make_unique(rownames(sce))
  colnames(sce@assays@data$counts) <- Make_unique(colnames(sce))
  sc_count <- sce@assays@data@listData$counts
  sc_meta <- data.frame(
    annotation = sce$annotation,
    sampleInfo = "sample1",
    batch = sce$batch,
    row.names = colnames(sc_count)
  )
  spatial_location <- data.frame(
    x = spe@colData$spatial_X,
    y = spe@colData$spatial_Y,
    row.names = colnames(spe)
  )
  spatial_count <- spe@assays@data@listData$counts
  card <- createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct.varname = "annotation",
    ct.select = unique(sc_meta$annotation),
    sample.varname = "batch",
    minCountSpot = 0,
    minCountGene = 0
  )
  card <- CARD_deconvolution(card)
  return(card)
}

Run_CARD <- function(scpath, sppath, outpath, raw) {
  sce <- Load_h5adsc_to_SCE(scpath, raw = raw)
  spe <- Load_h5adsp_to_SCE(sppath, raw = raw)
  rownames(sce) <- toupper(rownames(sce))
  rownames(spe) <- toupper(rownames(spe))
  rownames(sce@assays@data@listData$counts) <- toupper(rownames(sce@assays@data@listData$counts))
  rownames(spe@assays@data@listData$counts) <- toupper(rownames(spe@assays@data@listData$counts))
  card <- Deconvolution_CARD(spe, sce)
  res <- card@Proportion_CARD
  diff <- setdiff(colnames(spe), rownames(res))
  nrows <- length(diff)
  ncols <- length(colnames(res))
  # length of CARD results may less than colnames(spe)
  makeup <- matrix(data = 0, nrow = nrows, ncol = ncols)
  rownames(makeup) <- diff
  colnames(makeup) <- colnames(card@Proportion_CARD)
  mres <- rbind(makeup, res)
  mres <- mres[colnames(spe), ]
  mat <- t(mres)
  h5write(mat, outpath, "/obsm/CARD")
  message(paste0("Deconvolution with 'CARD' finished,
  idents saved at /obsm/CARD in", outpath))
}

cat("#### Running CARD ####\\n")
cat(paste0("Single-cell h5ad file:", args$scpath, "\\n"))
cat(paste0("Spatial h5ad file:", args$sppath, "\\n"))
cat(paste0("Output h5ad file:", outpath, "\\n"))
cat(paste0("Using raw:", as.character(raw), "\\n"))
Run_CARD(args$scpath, args$sppath, outpath, raw)
cat("#### CARD Finished ####\\n")
\`\`\`
`

const Card = () => {
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

export default Card;