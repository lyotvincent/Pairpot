import React from 'react';
import Markdown from 'react-markdown';
import ReactMarkdown from 'react-markdown';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter';
import { solarizedlight } from 'react-syntax-highlighter/dist/esm/styles/prism';

var Rcode = 
`
### RCTD
**Source:** https://github.com/dmcable/spacexr\n
**Citation:** Cable, D.M., Murray, E., Zou, L.S., Goeva, A., Macosko, E.Z., Chen, F., and Irizarry, R.A. (2022). Robust decomposition of cell type mixtures in spatial transcriptomics. Nat. Biotechnol. 40, 517-526.

### Key features and functions of spacexr include:
- Deconvolution of Spatial Data
- Inferring Cell Type Composition

### Rcode
**********
\`\`\`r 
# fetching params
library(getopt)
spec <- matrix(c(
  "help", "h", 0, "logical", "help",
  "scpath", "c", 1, "character", "path to single-cell h5ad file",
  "sppath", "p", 1, "character", "path to spatial h5ad file",
  "outpath", "o", 2, "character", "path to output h5ad file, default as sp path",
  "raw", "r", 2, "logical", "whether use raw data in h5ad file",
  "float_in_raw", "f", 2, "logical", "whether contains float in raw data."
), byrow = T, ncol = 5)

args <- getopt(spec)
if (!is.null(args$help) || is.null(args$scpath) || is.null(args$sppath)) {
  cat(paste(getopt(spec, usage = T), "\\n"))
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

float_in_raw <- args$float_in_raw
if (is.null(float_in_raw)) {
  float_in_raw <- FALSE
}

source("_loading.R")
library(rhdf5)
library(SingleCellExperiment)
library(Matrix)
Check_Load_GithubPackages("spacexr", "dmcable/spacexr")


Deconvolution_spacexr <- function(sr, ref) {
  rctd <- create.RCTD(sr, ref,
    max_cores = 7,
    UMI_min = 0,
    counts_MIN = 0,
    CELL_MIN_INSTANCE = 0
  )
  rctd <- run.RCTD(rctd)
  return(rctd)
}

GenerateRef_spacexr <- function(sce, float_in_raw) {
  cell_type <- as.factor(sce@colData@listData$annotation)
  names(cell_type) <- colnames(sce)
  counts <- sce@assays@data$counts
  if (float_in_raw) {
    counts <- floor(exp(as.matrix(counts)) - 1)
  }
  counts <- as(counts, "sparseMatrix")
  rownames(counts) <- Make_unique(rownames(sce))
  colnames(counts) <- Make_unique(colnames(sce))

  ref <- tryCatch({
    ref.0 <- spacexr::Reference(counts,
      cell_types = cell_type
    )
    return(ref.0)
  }, error = function(e) {
    cat("**Fetching error", e$message, '\\n')
    cat("**Trying do exp and floor.\\n")
    counts <- floor(exp(as.matrix(counts)) - 1)
    counts <- as(counts, "sparseMatrix")
    rownames(counts) <- Make_unique(rownames(sce))
    colnames(counts) <- Make_unique(colnames(sce))
    return(spacexr::Reference(counts,
      cell_types = cell_type,
      require_int = FALSE
    ))
  }, finally = {
    cat("**Ref for spacexr generated...")
  })
  return(ref)
}

Run_spacexr <- function(scpath, sppath, outpath, raw, float_in_raw) {
  sce <- Load_h5adsc_to_SCE(scpath, raw = raw)
  rownames(sce) <- toupper(rownames(sce))
  anno <- sce$annotation
  sce$annotation <- gsub("/", "^", anno)
  anno <- table(anno)
  spe <- Load_h5adsp_to_SCE(sppath, raw = raw)
  rownames(spe) <- toupper(rownames(spe))
  rownames(sce@assays@data@listData$counts) <- toupper(rownames(sce@assays@data@listData$counts))
  rownames(spe@assays@data@listData$counts) <- toupper(rownames(spe@assays@data@listData$counts))
  spatial_locs <- data.frame(x = spe$spatial_X, y = spe$spatial_Y, row.names = colnames(spe))
  sr <- SpatialRNA(
    coords = spatial_locs,
    counts = spe@assays@data@listData$counts,
    use_fake_coords = FALSE,
    require_int = FALSE
  )
  ref <- GenerateRef_spacexr(sce, float_in_raw)
  rctd <- Deconvolution_spacexr(sr, ref)
  rctd@results$weights@Dimnames[[2]] <- attr(anno, "names")
  results <- rctd@results
  norm_weights <- normalize_weights(results$weights)
  mat <- t(as.matrix(norm_weights))
  h5write(mat, outpath, "/obsm/RCTD")
  message(paste0("Deconvolution with 'spacexr' finished, idents saved at /obsm/RCTD in", outpath))
}

cat("#### Running RCTD ####\\n")
cat(paste0("Single-cell h5ad file:", args$scpath, "\\n"))
cat(paste0("Spatial h5ad file:", args$sppath, "\\n"))
cat(paste0("Output h5ad file:", outpath, "\\n"))
cat(paste0("Using raw:", as.character(raw), "\\n"))
Run_spacexr(args$scpath, args$sppath, outpath, raw, float_in_raw)
cat("#### RCTD Finished ####\\n")
\`\`\`
`

const Spacexr = () => {
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

export default Spacexr;