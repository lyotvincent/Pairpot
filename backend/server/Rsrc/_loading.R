library(rhdf5)
library(SingleCellExperiment)
library(Matrix)

Check_Load_BiocPackages <- function(pkgName) {
  if (!require(pkgName, character.only = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    library("BiocManager")
    BiocManager::install(pkgName, update = FALSE)
    library(pkgName, character.only = TRUE)
  }
}

Check_Load_GithubPackages <- function(pkgName, URL) {
  if (!require(pkgName, character.only = TRUE)) {
    if (!require("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    library("devtools")
    install_github(URL, upgrade = "default", quietly = TRUE, build_vignettes = F)
    library(pkgName, character.only = TRUE)
  }
}

Check_Load_InstallPackages <- function(pkgName) {
  if (!require(pkgName, character.only = TRUE)) {
    install.packages(pkgName, quiet = TRUE)
    library(pkgName, character.only = TRUE)
  }
}

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

Load_h5adsc_to_SCE <- function(scpath, raw = TRUE) {
  scmat <- h5read(scpath, "/")
  if (("Raw" %in% attr(scmat$layers, "names")) && raw) {
    cat("**Using raw data...\n")
    if (class(scmat$layers$Raw)[[1]] == "list") {
      cat("**Raw data is a sparse matrix...\n")
      dat <- sparseMatrix(
        i = scmat$layers$Raw$indices[] + 1,
        p = scmat$layers$Raw$indptr[],
        x = as.numeric(scmat$layers$Raw$data[]),
        repr = "C"
      )
    } else {
      cat("**Raw data is a dense matrix...\n")
      dat <- as(scmat$layers$Raw, "dgCMatrix")
    }
    h5ad.var <- scmat$var
  } else {
    cat("**Using processed data...\n")
    X <- scmat$X
    dat <- sparseMatrix(
      i = X$indices[] + 1,
      p = X$indptr[],
      x = as.numeric(X$data[]),
      repr = "C"
    )
    h5ad.var <- scmat$var
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
          var[[name]] <- factor(h5ad.var[[name]]$codes,
            labels = h5ad.var[[name]]$categories
          )
        } else {
          var[[name]] <- h5ad.var[[name]]
        }
      } else {
        var[[name]] <- h5ad.var[[name]]
      }
    }
  }
  cat("**Var loaded...\n")
  h5ad.obs <- scmat$obs
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
          obs[[name]] <- factor(h5ad.obs[[name]]$codes,
            labels = h5ad.obs[[name]]$categories
          )
        } else {
          obs[[name]] <- h5ad.obs[[name]]
        }
      } else {
        obs[[name]] <- h5ad.obs[[name]]
      }
    }
  }
  # some datasets uses barcodes instead of `_index`
  if (!("_index" %in% attr(h5ad.obs, "names"))) {
    obs$`_index` <- obs$`barcodes`
  }
  obs$`_index` <- Make_unique(obs$`_index`)
  obs_use <- data.frame(
    annotation = obs$annotation,
    batch = obs$batch,
    row.names = as.character(obs$`_index`)
  )
  cat("**Obs loaded...\n")
  dims <- c(length(h5ad.var[["_index"]]), length(obs$`_index`))
  dat@Dim <- dims
  dat@Dimnames <- list(
    as.character(h5ad.var[["_index"]]),
    as.character(obs$`_index`)
  )
  cat("**Dims loaded...\n")
  sce <- SingleCellExperiment(
    assays = list(counts = dat),
    rowData = DataFrame(data.frame(var[['_index']])),
    colData = DataFrame(data.frame(obs_use))
  )
  cat("**SCE constructed...\n")
  return(sce)
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
