#' Title Sometimes the slice is not in the correct position on the chip. This
#' function can make the data point rotate at a specific angle
#'
#' @param data the data has x and y column represent coordinate position.
#' @param angle the rotate angle, must be a pi number, like pi/2
#'
#' @return data
#' @export adjust_position
#'
#' @examples #
adjust_position <- function(data = NULL, angle = NULL){
  rotate_mat <- matrix( c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2 )
  adjust_mat = as.matrix(data[,c("x","y")]) %*% rotate_mat
  data[,1] = adjust_mat[,1]
  data[,2] = adjust_mat[,2]
  return(data)
}

#' Title Split rownames into the coordinate position of X and Y.
#'
#' @param data the data need its rownames including the postion x and y
#' @param pattern merge coord x and y use pattern as a split character,default is "a"
#' @return data
#' @export add_loc
#'
#' @examples #
add_loc <- function(data = NULL, pattern = "_"){
  #data = as.data.frame(data)
  rownames(data) = gsub(".*_","",rownames(data))
  data_loc <- cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(data),split=pattern),"[",1)),
                               y=as.numeric(sapply(strsplit(rownames(data),split=pattern),"[",2)),
                               data)
  return(data_loc)
}


#' Title read data from BGI stereo-seq, and creat a seurat object
#'
#' @param file the data csv file path
#' @param pattern merge coord x and y use pattern as a split character,default is "a"
#' @param seurat TRUE/FALSE creat or not creat a seurat object
#' @param name the object name
#' @return object
#' @export read_spat
#' @importFrom data.table fread
#' @importFrom Matrix sparseMatrix
#' @importFrom Seurat CreateSeuratObject
#'
#' @examples #
read_spat <- function(file,pattern = "_",seurat = TRUE,name = "Spatial"){
  data_BGI_raw = data.table::fread(file)
  splitline()
  messageline("Data imported but waiting for convert")
  data_BGI_raw$coord  = with(data_BGI_raw, paste0(x,pattern,y))
  # the gene expression
  exp = data_BGI_raw$MIDCounts
  #create gene index for sparse.matrix
  gene_id = unique(data_BGI_raw$geneID)
  index_gid = data.frame(row.names = gene_id, value = seq(length(unique(gene_id))))
  gid = index_gid[data_BGI_raw$geneID, "value"]
  #create coordinate index for sparse.matrix
  bin_id = unique(data_BGI_raw$coord)
  index_coord = data.frame(row.names = bin_id, value = seq(length(unique(bin_id))))
  bid = index_coord[data_BGI_raw$coord,"value"]
  # create sparse.matrix
  stereo_matrix_sparse <- Matrix::sparseMatrix(i = gid, j = bid, x = exp)
  rownames(stereo_matrix_sparse) <- rownames(index_gid)
  colnames(stereo_matrix_sparse) <- rownames(index_coord)
  # NEXT
  splitline()
  messageline("Data successfully converted")
  if(seurat == TRUE){
    stereo <- Seurat::CreateSeuratObject(counts = stereo_matrix_sparse, project = name, min.cells = 3, min.features = 200)
  }
  splitline()
  messageline("Successfully converted to Seurat object")
  return(stereo)
}
#' Title Retrieves data (feature expression, PCA scores, metrics, etc.) for a
#' set of cells in a Seurat object and add the coordinate position of X and Y of spatial chip.
#'
#' @param object the Seurat object
#' @param features the var want to select
#' @param cells the cell want to subset
#' @param slot must be one of "counts","data","scale.data"
#' @param ... inherit
#' @importFrom Seurat FetchData
#'
#' @return data
#' @export get_exp_loc
#'
#' @examples#
get_exp_loc <- function(object, features, cells = NULL, slot = "scale.data",...){
  data = Seurat::FetchData(object = object, vars = features, cells = cells, slot = slot,...)
  data_loc_exp = add_loc(data)
  return(data_loc_exp)
}

#' Title the packages environment need install or library
#'
#' @param pkgs which R packages you need library
#'
#' @return environment
#' @export load_spat_env
#' @importFrom utils txtProgressBar setTxtProgressBar installed.packages install.packages
#'
#' @examples #
load_spat_env <- function(pkgs = c("Seurat","Matrix","tidyverse","patchwork","data.table")){
  paknames = rownames(utils::installed.packages())
  need_pak = pkgs
  IF = vector()
  j = 0
  pb <- utils::txtProgressBar(min = 0, max = length(need_pak), style = 3)
  for(i in 1:length(need_pak)){
    packagename = need_pak[i]
    if(packagename %in% paknames){suppressPackageStartupMessages(require(packagename,character.only = TRUE,quietly = TRUE));IF[j] = TRUE } else {utils::install.packages(packagename,quiet = TRUE);suppressPackageStartupMessages(require(packagename,character.only = TRUE,quietly = TRUE));IF[j] = FALSE}
    message("\nPackages: <",packagename,">"," is OK")
    j = j + 1
    utils::setTxtProgressBar(pb, j)
  }
  if(all(IF)){
    cat("\nAll Packages installed and library !")
  }
  close(pb)
}

#' Title cool but useless
#'
#' @return charater
#' @export splitline
#'
#' @examples splitline()
splitline <- function() {
  width <- getOption("width")
  ws <- rep("=", floor(width))
  cat("\n",ws, sep = "")
}

#' Title nothing
#'
#' @param message you want say to user
#'
#' @return charater
#' @export messageline
#'
#' @examples messageline("yesimola !")
messageline <- function(message) {
  width <- getOption("width")
  mid <- paste("^_^   ",message,"   ^_^\n",sep = "")
  ws <- rep(" ", floor((width - nchar(mid))/2))
  cat(
    ws,
    mid,
    ws,
    sep = "")
}

# #' Title the github website:https://github.com/satijalab/seurat/issues/2833
# #'
# #' @param object the Seurat object and want to analysis in monocle3
# #'
# #' @importFrom monocle3 new_cell_data_set
# #' @importFrom utils data
# #' @return object
# #' @export seurat_to_monocle3
# #'
# #' @examples #
# seurat_to_monocle3 <- function(object){
#   stopifnot("The obect must be a Seurat"=class(object) ==  "Seurat")
#   gene_annotation <- as.data.frame(rownames(object@reductions[["pca"]]@feature.loadings),
#                                    row.names = rownames(object@reductions[["pca"]]@feature.loadings))
#   colnames(gene_annotation) <- "gene_short_name"
#   cell_metadata <- as.data.frame(object@assays[["RNA"]]@counts@Dimnames[[2]],
#                                  row.names = object@assays[["RNA"]]@counts@Dimnames[[2]])
#   colnames(cell_metadata) <- "barcode"
#   New_matrix <- object@assays[["RNA"]]@counts
#   New_matrix <- New_matrix[rownames(object@reductions[["pca"]]@feature.loadings), ]
#   expression_matrix <- New_matrix
#   cds_from_seurat <- monocle3::new_cell_data_set(expression_matrix,
#                                        cell_metadata = cell_metadata,
#                                        gene_metadata = gene_annotation)
#   recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
#   names(recreate.partition) <- cds_from_seurat@colData@rownames
#   recreate.partition <- as.factor(recreate.partition)
#   cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
#   list_cluster <- object@active.ident
#   names(list_cluster) <- object@assays[["RNA"]]@data@Dimnames[[2]]
#   cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
#   cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
#   cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-object@reductions[["umap"]]@cell.embeddings
#   cds_from_seurat@preprocess_aux$gene_loadings <- object@reductions[["pca"]]@feature.loadings
#   return(cds_from_seurat)
# }

#' Title the ggplot2 ggsave function changed for the dir path
#'
#' @param dir the dirpath you want to save in
#' @param filename the filename you want to name
#' @param width the width parameter in ggplot2::ggsave
#' @param height the height parameter in ggplot2::ggsave
#' @param dpi the dpi parameter in ggplot2::ggsave
#'
#' @return ggplot2 object
#' @export gs
#' @importFrom ggplot2 ggsave
#'
#' @examples #

gs <- function(dir,filename,width,height,dpi){
  ggplot2::ggsave(paste(dir, "/",filename,sep = ""),
                  width = width,height = height,dpi = dpi)
}

#' Title the St work created by I and qiujiaying,but the workflow don't suit.so make
#' a "loci2name" function for the object transfrom.
#'
#' @param object the seurat spatialObj
#' @importFrom Seurat RenameCells
#'
#' @return object the
#'
#' @export loci2name
#'
#' @examples #
loci2name <- function(object){
  splitline()
  messageline("Add The chip coord to the cell name")
  splitline()
  rowid = object@images[["slice1"]]@coordinates[["row"]]
  colid = object@images[["slice1"]]@coordinates[["col"]]
  coord  = paste0(rowid,"x",colid)
  object <- Seurat::RenameCells(object, new.names = coord)
  return(object)
}


#' Add the image to the seurat object for some software use coord and image information
#'
#' @param seuratObject a seurat object with no image s4 obejct
#'
#' @return seuratObject
#' @export add_image
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider
#' @importFrom methods new is
#'
#' @examples #
add_image <- function(seuratObject){
  #stopifnot("The obect must be a Seurat"=class(object) ==  "Seurat")
  spotfactor = list(spot =1,fiducial=1, hires=1,lowres=1)
  class(spotfactor) = "scalefactors"
  position =as.data.frame(seuratObject@reductions[["spatial"]]@cell.embeddings)
  #----------------
  ncounts = as.data.frame(seuratObject@meta.data[["total_counts"]])
  ncounts = cbind(position,ncounts)
  colnames(ncounts) = c("y","x","ncount")
  ncounts$ncount = round((ncounts$ncount - min(ncounts$ncount))/(max(ncounts$ncount)-min(ncounts$ncount)),digits = 4)
  ncounts = round(ncounts)
  img = tidyr::pivot_wider(data = ncounts,x,names_from = y,values_from = ncount)
  img = tibble::column_to_rownames(img,"x")
  img = as.matrix(img)
  img[is.na(img)] = 0
  img = img[order(as.numeric(rownames(img))),]
  img = img[,order(as.numeric(colnames(img)))]
  img = 1 - img
  img <- array(rep(img,3),dim = c(nrow(img),ncol(img),3))
  #grid::grid.raster(img_array, interpolate=FALSE)
  #---------------------
  position = position[,c(1,2,1,2)]
  colnames(position) = c("col","row","imagecol","imagerow")
  spotradius =round(1.71/(length(unique(position$row))+length(unique(position$col))),digits = 3)
  img <- new(
    Class = 'VisiumV1',
    image = img,
    coordinates = position,
    scale.factors = spotfactor,
    spot.radius = spotradius
  )
  img2 = list(slice1 = img)
  seuratObject@images = img2
  #-------
  seuratObject@images[["slice1"]]@key = "slice1_"
  seuratObject@images[["slice1"]]@assay = "RNA"
  seuratObject@images[["slice1"]]@coordinates[["tissue"]] = rep(1,nrow(position))
  #seuratObject@assays[["RNA"]]@counts = as.sparse(seuratObject@assays[["RNA"]]@counts)
  #seuratObject@assays[["RNA"]]@data = as.sparse(seuratObject@assays[["RNA"]]@data)
  seuratObject@active.ident = seuratObject$louvain_clusters
  return(seuratObject)
}
