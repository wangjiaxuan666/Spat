#' Title: Draw a cell clustering map at the chip position
#'
#' @param object the Seurat object
#' @param cols the color vector to mapping
#' @param rotate the rotate angle, must be a pi number, like pi/2
#' @importFrom  ggplot2 ggplot aes geom_point theme_void labs coord_cartesian scale_color_manual
#' @return ggplot
#' @export cell_in_chip
#'
#' @examples #
cell_in_chip <- function(object = NULL, cols = NULL, rotate = NULL){
  data_loc = object@active.ident
  # 提取分群结果
  data_loc = add_loc(data = data_loc)

  if(!is.null(rotate)){
    data_loc = adjust_position(data_loc, angle = rotate)
  }

  # 设置范围-------------------
  slice_coord_x_min = min(data_loc$x)
  slice_coord_x_max = max(data_loc$x)
  slice_coord_y_min = min(data_loc$y)
  slice_coord_y_max = max(data_loc$y)

  ggplot2::ggplot(data_loc)+
    ggplot2::geom_point(ggplot2::aes(x = x,y = y,color = data),size = 3)+
    ggplot2::theme_void()+
    ggplot2::labs(color = "Cell Type")+
    ggplot2::coord_cartesian(xlim = c(slice_coord_x_min,slice_coord_x_max),
                    ylim = c(slice_coord_y_min,slice_coord_y_max)) -> plot_result
  if(!is.null(cols)){
    plot_result = plot_result + ggplot2::scale_color_manual(values = cols)
  }
  return(plot_result)
}


#' Title: heatmap on chip
#'
#' @param object the Seurat object
#' @param featrues the gene ID want to display
#' @param cells the cell want to subset
#' @param slot must be one of "counts","data","scale.data"
#' @param rotate the rotate angle, must be a pi number, like pi/2
#' @param cols the color vector to mapping
#' @param ... inherit
#' @importFrom  ggplot2 ggplot aes geom_point theme_void labs coord_cartesian scale_color_manual
#'
#' @return ggplot
#' @export exp_in_chip
#'
#' @examples #
exp_in_chip <- function(object, featrues, cells = NULL, slot = "scale.data", rotate = NULL,cols = c("white","red","black"),...){
  
  data_loc_exp = get_exp_loc(object = object, featrues = featrues, cells = cells, slot = slot,...)
  
  if(!is.null(rotate)){
    data_loc_exp = adjust_position(data_loc_exp, angle = rotate)
  }
  
  # 设置范围-------------------
  slice_coord_x_min = min(data_loc_exp$x)
  slice_coord_x_max = max(data_loc_exp$x)
  slice_coord_y_min = min(data_loc_exp$y)
  slice_coord_y_max = max(data_loc_exp$y)
  colnames(data_loc_exp)[3] = "features"
  ggplot2::ggplot(data_loc_exp)+
    ggplot2::geom_point(aes(x = x,y = y,color = features),size = 3)+
    ggplot2::theme_void()+
    ggplot2::labs(color = "Gene Expression Level")+
    ggplot2::coord_cartesian(xlim = c(slice_coord_x_min,slice_coord_x_max),
                             ylim = c(slice_coord_y_min,slice_coord_y_max)) -> plot_result
  if(!is.null(cols)){
    plot_result = plot_result + ggplot2::scale_color_gradientn(colors = cols)
  }
  return(plot_result)
}