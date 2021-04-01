#' Title: Sometimes the slice is not in the correct position on the chip. This
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
add_loc <- function(data = NULL, pattern = "a"){
  data = as.data.frame(data)
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
#'
#' @examples #
read_spat <- function(file,pattern = "a",seurat = TRUE,name = "Spatial"){
  data_BGI_raw = data.table::fread(file)
  splitline()
  messageline("Data imported successfully")
  data_BGI_raw$loci  = with(data_BGI_raw, paste0(x,pattern,y))
  data_BGI_raw = data_BGI_raw[,c(1,4,5)]
  data_matrix_raw =  data.table::dcast(data_BGI_raw, geneID~loci,value.var = "MIDCounts")
  data_matrix_raw =  data.table::setDF(data_matrix_raw, rownames = data_matrix_raw$geneID)
  data_matrix_raw = data_matrix_raw[,-1]
  data_matrix_raw[is.na(data_matrix_raw)] <- 0
  stereo = data_matrix_raw
  splitline()
  messageline("Data successfully converted")
  if(seurat == TRUE){
    stereo_matrix_sparse <- Seurat::as.sparse(data_matrix_raw)
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
#' @param featrues the var want to select
#' @param cells the cell want to subset
#' @param slot must be one of "counts","data","scale.data"
#' @param ... inherit
#' @importFrom Seurat FetchData
#'
#' @return data
#' @export get_exp_loc
#'
#' @examples#
get_exp_loc <- function(object, featrues, cells = NULL, slot = "scale.data",...){
  data = Seurat::FetchData(object = object, vars = featrues, cells = cells, slot = slot,...)
  data_loc_exp = add_loc(data)
  return(data_loc_exp)
}

#' Title the packages environment need install or library
#'
#' @return environment
#' @export load_spat_env
#' @importFrom utils txtProgressBar setTxtProgressBar installed.packages install.packages
#'
#' @examples #
load_spat_env <- function(){
  paknames = rownames(utils::installed.packages())
  need_pak = c("Seurat","tidyverse","data.table","patchwork")
  IF = vector()
  j = 0
  pb <- utils::txtProgressBar(min = 0, max = length(need_pak), style = 3)
  for(i in need_pak){
    if(i %in% paknames){require(i);IF[j] = TRUE } else {utils::install.packages(i);require(i);IF[j] = FALSE}
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
#' @examples #
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
#' @examples #
messageline <- function(message) {
  width <- getOption("width")
  mid <- paste0("^_^   ",message,"   ^_^\n",sep = "")
  ws <- rep(" ", floor((width - nchar(mid))/2))
  cat(ws, mid, ws, sep = "")
}
