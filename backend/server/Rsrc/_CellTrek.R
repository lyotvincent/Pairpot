# fetching params
library(getopt)
spec <- matrix(c(
  "help", "h", 0, "logical", "help",
  "scpath", "c", 1, "character", "path to single-cell h5ad file",
  "sppath", "p", 1, "character", "path to spatial h5ad file",
  "outpath", "o", 2, "character", "path to output h5ad file, default as sp path",
  "raw", "r", 2, "logical", "whether use raw data in h5ad file"
), byrow = T, ncol = 5)

args <- getopt(spec)
if (!is.null(args$help) || is.null(args$scpath) || is.null(args$sppath)) {
  cat(paste(getopt(spec, usage = T), "\n"))
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

options(stringsAsFactors = F)
library(scater)
source("_loading.R")
library(dplyr)
library(rhdf5)
library(spatstat.explore)
library(Seurat)
library(viridis)
library(magrittr)
library(SeuratObject)
library(ConsensusClusterPlus)

mytraint <- function(st_data, sc_data, st_assay = "Spatial", sc_assay = "scint", norm = "LogNormalize", nfeatures = 2000,
                     cell_names = "cell_names", coord_xy = c("spatial_X", "spatial_Y"), gene_kept = NULL, ...) {
  st_data$id <- names(st_data$orig.ident)
  sc_data$id <- names(sc_data$orig.ident)
  sc_data$cell_names <- make.names(sc_data@meta.data[, cell_names])
  st_data$type <- "st"
  sc_data$type <- "sc"
  st_data$coord_x <- st_data@images[[1]]@coordinates[, coord_xy[1]]
  st_data$coord_y <- st_data@images[[1]]@coordinates[, coord_xy[2]]
  DefaultAssay(st_data) <- st_assay
  DefaultAssay(sc_data) <- sc_assay

  cat("Finding transfer anchors... \n")
  st_idx <- st_data$id
  sc_idx <- sc_data$id

  ## Integration features ##
  sc_st_list <- list(st_data = st_data, sc_data = sc_data)
  sc_st_features <- Seurat::SelectIntegrationFeatures(sc_st_list, nfeatures = nfeatures)
  if (!is.null(gene_kept)) {
    sc_st_features <- union(sc_st_features, gene_kept)
  }

  sc_st_features <- sc_st_features[(sc_st_features %in% rownames(st_data[[st_assay]]@data)) &
    (sc_st_features %in% rownames(sc_data[[sc_assay]]@data))]
  cat("Using", length(sc_st_features), "features for integration... \n")
  ###

  sc_st_anchors <- Seurat::FindTransferAnchors(
    reference = sc_data, query = st_data,
    reference.assay = sc_assay, query.assay = st_assay,
    normalization.method = norm, features = sc_st_features, reduction = "cca", ...
  )

  cat("Data transfering... \n")
  st_data_trans <- Seurat::TransferData(
    anchorset = sc_st_anchors,
    refdata = GetAssayData(sc_data, assay = sc_assay, slot = "data")[sc_st_features, ], weight.reduction = "cca"
  )
  st_data@assays$transfer <- st_data_trans

  cat("Creating new Seurat object... \n")
  sc_st_meta <- dplyr::bind_rows(st_data@meta.data, sc_data@meta.data)
  counts_temp <- cbind(data.frame(st_data[["transfer"]]@data), data.frame(sc_data[[sc_assay]]@data[sc_st_features, ] %>% data.frame()))
  rownames(sc_st_meta) <- make.names(sc_st_meta$id)
  colnames(counts_temp) <- make.names(sc_st_meta$id)
  sc_st_int <- CreateSeuratObject(counts = counts_temp, assay = "traint", meta.data = sc_st_meta)
  # sc_st_int[['traint']]@data <- sc_st_int[['traint']]@counts
  # sc_st_int[['traint']]@counts <- matrix(NA, nrow = 0, ncol = 0)

  cat("Scaling -> PCA -> UMAP... \n")
  sc_st_int <- NormalizeData(sc_st_int)
  sc_st_int <- ScaleData(sc_st_int, features = sc_st_features) %>%
    RunPCA(features = sc_st_features)
  sc_st_int <- RunUMAP(sc_st_int, dims = 1:30)
  sc_st_int@images <- st_data@images
  sc_st_int@images[[1]]@coordinates <- data.frame(
    imagerow = sc_st_int@meta.data$coord_x,
    imagecol = sc_st_int@meta.data$coord_y,
    row.names = rownames(sc_st_int@meta.data)
  )
  return(sc_st_int)
}

celltrek_repel <- function(celltrek_inp, repel_r = 5, repel_iter = 10) {
  celltrek_dr_raw <- Embeddings(celltrek_inp, "celltrek_raw")
  celltrek_dr <- Embeddings(celltrek_inp, "celltrek")
  repel_input <- data.frame(celltrek_dr_raw, repel_r = repel_r)

  ## Add noise ##
  theta <- runif(nrow(celltrek_dr), 0, 2 * pi)
  alpha <- sqrt(runif(nrow(celltrek_dr), 0, 1))
  repel_input[, 1] <- repel_input[, 1] + sin(theta) * alpha * repel_r
  repel_input[, 2] <- repel_input[, 2] + cos(theta) * alpha * repel_r

  ## Repelling ##
  cat("Repelling points...\n")
  celltrek_repel <- circleRepelLayout(repel_input, sizetype = "radius", maxiter = repel_iter)
  celltrek_dr[, 1] <- celltrek_repel$layout$x
  celltrek_dr[, 2] <- celltrek_repel$layout$y
  celltrek_out <- celltrek_inp
  celltrek_out@reductions$celltrek@cell.embeddings <- celltrek_dr
  return(celltrek_out)
}

#' Calculate the RF-distance between sc and st
#'
#' @param st_sc_int Seurat traint object
#' @param int_assay Name of integration assay
#' @param reduction Dimension reduction method used, usually pca
#' @param intp If TRUE, do interpolation
#' @param intp_pnt Interpolation point number
#' @param intp_lin If TRUE, use linear interpolation
#' @param nPCs Number of PCs used for CellTrek
#' @param ntree Number of trees in random forest
#' @param keep_model If TRUE, return the trained random forest model
#'
#' @return A list of 1. celltrek_distance matrix; 2. trained random forest model (optional)
#'
#' @import dbscan
#' @importFrom akima interpp
#' @import magrittr
#' @import dplyr
#' @import randomForestSRC
#'
#' @examples dist_test <- celltrek_dist(st_sc_int = st_sc_int, int_assay = "traint", reduction = "pca", intp = T, intp_pnt = 10000, intp_lin = F, nPCs = 30, ntree = 1000, keep_model = T)
celltrek_dist <- function(st_sc_int, int_assay = "traint", reduction = "pca", intp = T, intp_pnt = 10000, intp_lin = F, nPCs = 30, ntree = 1000, keep_model = T) {
  DefaultAssay(st_sc_int) <- int_assay
  kNN_dist <- dbscan::kNN(na.omit(st_sc_int@meta.data[, c("coord_x", "coord_y")]), k = 6)$dist
  spot_dis <- median(kNN_dist) %>% round()
  cat("Distance between spots is:", spot_dis, "\n")

  st_sc_int$id <- names(st_sc_int$orig.ident)
  st_idx <- st_sc_int$id[st_sc_int$type == "st"]
  sc_idx <- st_sc_int$id[st_sc_int$type == "sc"]
  meta_df <- data.frame(st_sc_int@meta.data)

  st_sc_int_pca <- st_sc_int@reductions[[reduction]]@cell.embeddings[, 1:nPCs] %>%
    data.frame() %>%
    mutate(
      id = st_sc_int$id, type = st_sc_int$type, class = st_sc_int$cell_names,
      coord_x = st_sc_int$coord_x, coord_y = st_sc_int$coord_y
    )
  st_pca <- st_sc_int_pca %>%
    dplyr::filter(type == "st") %>%
    dplyr::select(-c(id:class))

  ## Interpolation ##
  ## Uniform sampling ##
  if (intp) {
    cat("Interpolating...\n")
    spot_ratio <- intp_pnt / nrow(st_pca)
    st_intp_df <- apply(st_pca[, c("coord_x", "coord_y")], 1, function(row_x) {
      runif_test <- runif(1)
      if (runif_test < spot_ratio %% 1) {
        theta <- runif(ceiling(spot_ratio), 0, 2 * pi)
        alpha <- sqrt(runif(ceiling(spot_ratio), 0, 1))
        coord_x <- row_x[1] + (spot_dis / 2) * sin(theta) * alpha
        coord_y <- row_x[2] + (spot_dis / 2) * cos(theta) * alpha
      } else {
        theta <- runif(floor(spot_ratio), 0, 2 * pi)
        alpha <- sqrt(runif(floor(spot_ratio), 0, 1))
        coord_x <- row_x[1] + (spot_dis / 2) * sin(theta) * alpha
        coord_y <- row_x[2] + (spot_dis / 2) * cos(theta) * alpha
      }
      data.frame(coord_x, coord_y)
    }) %>% Reduce(rbind, .)

    st_intp_df <- apply(st_pca[, 1:nPCs], 2, function(col_x) {
      akima::interpp(
        x = st_pca$coord_x, y = st_pca$coord_y, z = col_x,
        linear = intp_lin, xo = st_intp_df$coord_x, yo = st_intp_df$coord_y
      ) %>%
        magrittr::extract2("z")
    }) %>%
      data.frame(., id = "X", type = "st_intp", st_intp_df) %>%
      na.omit()
    st_intp_df$id <- make.names(st_intp_df$id, unique = T)
    st_sc_int_pca <- bind_rows(st_sc_int_pca, st_intp_df)
  }

  cat("Random Forest training... \n")
  ## Training on ST ##
  data_train <- st_sc_int_pca %>%
    dplyr::filter(type == "st") %>%
    dplyr::select(-c(id:class))
  rf_train <- randomForestSRC::rfsrc(Multivar(coord_x, coord_y) ~ ., data_train, block.size = 5, ntree = ntree)

  cat("Random Forest prediction...  \n")
  ## Testing on all ##
  data_test <- st_sc_int_pca
  rf_pred <- randomForestSRC::predict.rfsrc(rf_train, newdata = data_test[, c(1:nPCs)], distance = "all")

  cat("Making distance matrix... \n")
  rf_pred_dist <- rf_pred$distance[data_test$type == "sc", data_test$type != "sc"] %>%
    set_rownames(data_test$id[data_test$type == "sc"]) %>%
    set_colnames(data_test$id[data_test$type != "sc"])

  output <- list()
  output$spot_d <- spot_dis
  output$celltrek_dist <- rf_pred_dist
  output$coord_df <- st_sc_int_pca[, c("id", "type", "coord_x", "coord_y")] %>%
    dplyr::filter(type != "sc") %>%
    magrittr::set_rownames(.$id) %>%
    dplyr::select(-id)
  if (keep_model) {
    output$model <- rf_train
  }
  return(output)
}

#'
#'
#' @param dist_mat Distance matrix of sc-st (sc in rows and st in columns)
#' @param coord_df Coordinates data frame of st (must contain coord_x, coord_y columns, barcode rownames)
#' @param dist_cut Distance cutoff
#' @param top_spot Maximum number of spots that one cell can be charted
#' @param spot_n Maximum number of cells that one spot can contain
#' @param repel_r Repelling radius
#' @param repel_iter Repelling iterations
#'
#' @return SC coordinates
#'
#' @import data.table
#' @import scales
#' @import dplyr
#' @importFrom packcircles circleRepelLayout
#'
#' @examples
celltrek_chart <- function(dist_mat, coord_df, dist_cut = 500, top_spot = 10, spot_n = 10, repel_r = 5, repel_iter = 10) {
  cat("Making graph... \n")
  dist_mat[dist_mat > dist_cut] <- NA
  dist_mat_dt <- data.table::data.table(Var1 = rownames(dist_mat), dist_mat)
  dist_edge_list <- data.table::melt(dist_mat_dt, id = 1, na.rm = T)
  colnames(dist_edge_list) <- c("Var1", "Var2", "value")
  dist_edge_list$val_rsc <- scales::rescale(dist_edge_list$value, to = c(0, repel_r))
  dist_edge_list$Var1 %<>% as.character
  dist_edge_list$Var2 %<>% as.character
  dist_edge_list$Var1_type <- "sc"
  dist_edge_list$Var2_type <- "non-sc"

  cat("Pruning graph...\n")
  dist_edge_list_sub <- dplyr::inner_join(
    dist_edge_list %>% group_by(Var1) %>% top_n(n = top_spot, wt = -value),
    dist_edge_list %>% group_by(Var2) %>% top_n(n = spot_n, wt = -value)
  ) %>% data.frame()

  cat("Spatial Charting SC data...\n")
  sc_coord <- sc_coord_raw <- data.frame(id_raw = dist_edge_list_sub$Var1, id_new = make.names(dist_edge_list_sub$Var1, unique = T)) # nolint
  sc_coord$coord_x <- sc_coord_raw$coord_x <- coord_df$coord_x[match(dist_edge_list_sub$Var2, rownames(coord_df))]
  sc_coord$coord_y <- sc_coord_raw$coord_y <- coord_df$coord_y[match(dist_edge_list_sub$Var2, rownames(coord_df))]
  ## Add noise ##
  theta <- runif(nrow(dist_edge_list_sub), 0, 2 * pi)
  alpha <- sqrt(runif(nrow(dist_edge_list_sub), 0, 1))
  sc_coord$coord_x <- sc_coord$coord_x + dist_edge_list_sub$val_rsc * sin(theta) * alpha
  sc_coord$coord_y <- sc_coord$coord_y + dist_edge_list_sub$val_rsc * cos(theta) * alpha
  ## Point repelling ##
  cat("Repelling points...\n")
  sc_repel_input <- data.frame(sc_coord[, c("coord_x", "coord_y")], repel_r = repel_r)
  sc_repel <- packcircles::circleRepelLayout(sc_repel_input, sizetype = "radius", maxiter = repel_iter)
  sc_coord$coord_x <- sc_repel$layout$x
  sc_coord$coord_y <- sc_repel$layout$y
  return(list(sc_coord_raw, sc_coord))
}

mycelltrek <- function(st_sc_int, int_assay = "traint", sc_data = NULL, sc_assay = "RNA", reduction = "pca",
                       intp = T, intp_pnt = 10000, intp_lin = F, nPCs = 30,
                       ntree = 1000, dist_thresh = .4, top_spot = 10, spot_n = 10, repel_r = 5, repel_iter = 10, keep_model = F, ...) {
  dist_res <- celltrek_dist(st_sc_int = st_sc_int, int_assay = int_assay, reduction = reduction, intp = intp, intp_pnt = intp_pnt, intp_lin = intp_lin, nPCs = nPCs, ntree = ntree, keep_model = T)
  spot_dis_intp <- median(unlist(dbscan::kNN(dist_res$coord_df[, c("coord_x", "coord_y")], k = 4)$dist)) # nolint # nolint: infix_spaces_linter.
  if (is.null(repel_r)) {
    repel_r <- spot_dis_intp / 4
  }
  sc_coord_list <- celltrek_chart(dist_mat = dist_res$celltrek_dist, coord_df = dist_res$coord_df, dist_cut = ntree * dist_thresh, top_spot = top_spot, spot_n = spot_n, repel_r = repel_r, repel_iter = repel_iter) # nolint: line_length_linter, infix_spaces_linter.
  sc_coord_raw <- sc_coord_list[[1]]
  sc_coord <- sc_coord_list[[2]]
  return(sc_coord)
}


Run_CellTrek <- function(scpath, sppath, outpath, raw) {
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

  batches <- as.character(unique(spseu$batch))
  for (bat in batches) {
    cat(paste0("Spatial Mapping at ", bat))
    spseu0 <- subset(x = spseu, subset = batch == bat)
    print(colnames(spseu0@images[[1]]@coordinates))
    st_sc_int <- mytraint(
      st_data = spseu0,
      sc_data = scseu,
      st_assay = "originalexp",
      sc_assay = "originalexp",
      coord_xy = c("x", "y"),
      cell_names = "annotation"
    )
    sc_coord <- mycelltrek(st_sc_int)
    h5write(sc_coord, outpath, paste0("/uns/mapping_", bat))
  }
  message(paste0("Spatial mapping with `CellTrek` finished, coords saved at /uns/mapping_batches in", outpath)) # nolint
}

cat("#### Running CellTrek #### \n")
cat(paste0("Single-cell h5ad file:", args$scpath, "\n"))
cat(paste0("Spatial h5ad file:", args$sppath, "\n"))
cat(paste0("Output h5ad file:", outpath, "\n"))
cat(paste0("Using raw:", as.character(raw), "\n"))
Run_CellTrek(args$scpath, args$sppath, outpath, raw)
cat("#### CellTrek Finished #### \n")
