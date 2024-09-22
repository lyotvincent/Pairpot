# fetching params
library(getopt)
spec <- matrix(c(
  "help", "h", 0, "logical", "help",
  "scpath", "c", 1, "character", "path to single-cell h5ad file",
  "sppath", "p", 1, "character", "path to spatial h5ad file",
  "outpath", "o", 2, "character", "path to output h5ad file, default as sp path",
  "species", "s", 2, "character", "Species for SpaTalk, Human or Mouse.",
  "raw", "r", 2, "logical", "whether use raw data in h5ad file"
), byrow = TRUE, ncol = 5)

args <- getopt(spec)
if (!is.null(args$help) || is.null(args$scpath) || is.null(args$sppath)) {
  cat(paste(getopt(spec, usage = TRUE), "\n"))
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
species <- args$species
if (is.null(species)) {
  species <- "Mouse"
}

library(spatstat.explore)
library(SpaTalk)
source("_loading.R")

Deconv_SpaTalk <- function(spe, sce, species = "Mouse", is_cell = F) {
  # create SpaTalk Object from SpatialExperiment and SingleCellExperiment object
  st_data <- spe@assays@data@listData$counts
  colnames(st_data) <- Make_unique(colnames(spe))
  rownames(st_data) <- Make_unique(rownames(spe))
  sc_data <- sce@assays@data@listData$counts
  rownames(sc_data) <- Make_unique(rownames(sce))
  colnames(sc_data) <- Make_unique(colnames(sce))
  sc_celltype <- as.character(sce$annotation)
  if (is_cell) {
    st_meta <- data.frame(
      cell = colnames(spe),
      x = spe$spatial_X,
      y = spe$spatial_Y
    )
  } else {
    st_meta <- data.frame(
      spot = colnames(spe),
      x = spe$spatial_X,
      y = spe$spatial_Y
    )
  }

  spaObj <- createSpaTalk(
    st_data = st_data,
    st_meta = st_meta,
    species = species,
    if_st_is_sc = is_cell,
    spot_max_cell = 5
  )
  spaObj <- dec_celltype(spaObj, sc_data, sc_celltype, use_n_cores = 16, iter_num = 500)
  return(spaObj)
}

Run_SpaTalk <- function(scpath, sppath, outpath, species, raw) {
  sce <- Load_h5adsc_to_SCE(scpath, raw = raw)
  spe <- Load_h5adsp_to_SCE(sppath, raw = raw)
  rownames(sce) <- toupper(rownames(sce))
  rownames(spe) <- toupper(rownames(spe))
  rownames(sce@assays@data@listData$counts) <- toupper(rownames(sce@assays@data@listData$counts))
  rownames(spe@assays@data@listData$counts) <- toupper(rownames(spe@assays@data@listData$counts))
  spaObj <- Deconv_SpaTalk(spe, sce, species)
  mat <- t(as.matrix(spaObj@meta$rawmeta[, 8:ncol(spaObj@meta$rawmeta)]))
  h5write(mat, outpath, "/obsm/SpaTalk")
  message(paste0("Deconvolution with `SpaTalk` finished, 
  idents saved at /obsm/SpaTalk in", outpath))
}

cat("#### Running SpaTalk-deconv ####\n")
cat(paste0("Single-cell h5ad file:", args$scpath, "\n"))
cat(paste0("Spatial h5ad file:", args$sppath, "\n"))
cat(paste0("Output h5ad file:", outpath, "\n"))
cat(paste0("Using raw:", as.character(raw), "\n"))
cat(paste0("Species:", species, "\n"))
Run_SpaTalk(args$scpath, args$sppath, outpath, species, raw)
cat("#### SpaTalk-deconv Finished ####\n")
