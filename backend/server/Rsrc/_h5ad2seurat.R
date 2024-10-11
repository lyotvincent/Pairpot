library(rhdf5)
library(SingleCellExperiment)
library(Matrix)
library(scater)

Make_unique <- function(arr) {
  dep_arr <- list()
  dep_idx <- which(duplicated(arr))
  for (idx in dep_idx) {
    item <- arr[idx]
    if (!(item %in% attr(dep_arr, "names"))) {
      dep_arr[item] <- 1
    }
    count <- dep_arr[[item]]
    arr[idx] <- paste0(arr[idx], ".", count)
    dep_arr[[item]] <- dep_arr[[item]] + 1
  }
  return(arr)
}

Load_h5adsp_to_SCE <- function(sppath, raw = TRUE) {
  spmat <- h5read(sppath, "/")
  if ("Raw" %in% attr(spmat$layers, "names") && raw) {
    cat("**Using raw data...\n")
    if (class(spmat$layers$Raw)[[1]] == "list") {
      cat("**Raw data is a sparse matrix...\n")
      dat <- sparseMatrix(
        i = spmat$layers$Raw$indices[] + 1,
        p = spmat$layers$Raw$indptr[],
        x = as.numeric(spmat$layers$Raw$data[]),
        repr = "C"
      )
    } else {
      cat("**Raw data is a dense matrix...\n")
      dat <- as(spmat$layers$Raw, "dgCMatrix")
    }
    h5ad.var <- spmat$var
  } else {
    cat("**Using processed data...\n")
    X <- spmat$X
    dat <- sparseMatrix(
      i = X$indices[] + 1,
      p = X$indptr[],
      x = as.numeric(X$data[]),
      repr = "C"
    )
    h5ad.var <- spmat$var
  }
  cat("**Core data loaded...\n")

  # generate factors using categories
  var <- list()
  if ("__categories" %in% attr(h5ad.var, "names")) { # old anndata
    for (name in attr(h5ad.var[["__categories"]], "names")) {
      if (length(h5ad.var[[name]]) >= length(h5ad.var[["__categories"]][[name]])) {
        var[[name]] <- factor(h5ad.var[[name]], labels = h5ad.var[["__categories"]][[name]])
      }
    }
  } else {
    for (name in attr(h5ad.var, "names")) {
      if (name != "_index") {
        if (class(var[[name]]) == "list") {
          var[[name]] <- factor(h5ad.var[[name]]$codes, labels = h5ad.var[[name]]$categories)
        } else {
          var[[name]] <- h5ad.var[[name]]
        }
      } else {
        var[[name]] <- h5ad.var[[name]]
      }
    }
  }
  cat("**Var loaded...\n")

  h5ad.obs <- spmat$obs
  obs <- list()
  if ("__categories" %in% attr(h5ad.obs, "names")) { # old anndata
    for (name in attr(h5ad.obs[["__categories"]], "names")) {
      if (length(h5ad.obs[[name]]) >= length(h5ad.obs[["__categories"]][[name]])) {
        obs[[name]] <- factor(h5ad.obs[[name]], labels = h5ad.obs[["__categories"]][[name]])
      }
    }
  } else { # new anndata
    for (name in attr(h5ad.obs, "names")) {
      if (name != "_index") {
        if (class(h5ad.obs[[name]]) == "list") {
          obs[[name]] <- factor(h5ad.obs[[name]]$codes, labels = h5ad.obs[[name]]$categories)
        } else {
          obs[[name]] <- h5ad.obs[[name]]
        }
      } else {
        obs[[name]] <- h5ad.obs[[name]]
      }
    }
  }
  cat("**Obs loaded...\n")

  h5ad.obsm <- spmat$obsm
  coord <- data.frame(t(h5ad.obsm$spatial[]))
  obs[["spatial_X"]] <- coord[[1]]
  obs[["spatial_Y"]] <- coord[[2]]
  obs$`_index` <- Make_unique(obs$`_index`)
  print(length(obs$`_index`))
  print(length(obs$annotation))
  print(length(obs$batch))
  obs_use <- data.frame(
    annotation = obs$annotation,
    batch = obs$batch,
    spatial_X = as.numeric(coord[[1]]),
    spatial_Y = as.numeric(coord[[2]]),
    row.names = obs$`_index`
  )
  cat("**Obsm loaded...\n")

  dims <- c(length(h5ad.var[["_index"]]), length(obs$`_index`))
  dat@Dim <- dims
  dat@Dimnames <- list(as.character(h5ad.var[["_index"]]), as.character(obs$`_index`))
  cat("**Dims loaded...\n")

  spe <- SingleCellExperiment(assays = list(counts = dat),
                              rowData = DataFrame(data.frame(var[['_index']])),
                              colData = DataFrame(obs_use))
  cat("**SPE constructed...\n")
  return(spe)
}

Load_h5adsp_to_Seurat <- function(sppath, raw = TRUE) {
  spe <- Load_h5adsp_to_SCE(sppath, raw)
  rownames(spe) <- toupper(rownames(spe))
  spe <- logNormCounts(spe)
  colnames(spe) <- Make_unique(colnames(spe))
  spseu <- SeuratObject::as.Seurat(x = spe, counts = "counts", data = "logcounts")
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
  return(spseu)
}