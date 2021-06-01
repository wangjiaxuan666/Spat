#' Title: Draw a cell clustering map at the chip position
#'
#' @param object the Seurat object
#' @param cols the color vector to mapping
#' @param rotate the rotate angle, must be a pi number, like pi/2
#' @param ... parameter
#' @param size the point size
#' @importFrom  ggplot2 ggplot aes geom_point theme_void labs coord_cartesian scale_color_manual
#' @return ggplot
#' @export cell_in_chip
#'
#' @examples #
cell_in_chip <- function(object = NULL, cols = NULL, rotate = NULL, size = 2, ...){
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
    ggplot2::geom_point(ggplot2::aes(x = x,y = y,color = data),size = size, ...)+
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
#' @param features the gene ID want to display
#' @param cells the cell want to subset
#' @param slot must be one of "counts","data","scale.data"
#' @param rotate the rotate angle, must be a pi number, like pi/2
#' @param cols the color vector to mapping
#' @param size the point size
#' @param ... inherit
#' @importFrom  ggplot2 ggplot aes geom_point theme_void labs coord_cartesian scale_color_manual theme
#'
#' @return ggplot
#' @export exp_in_chip
#'
#' @examples #
exp_in_chip <- function(object, features, cells = NULL, slot = "scale.data", rotate = NULL, size = 2, cols = c("#F3EDEC","red","black"),...){
  
  data_loc_exp = get_exp_loc(object = object, features = features, cells = cells, slot = slot,...)
  
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
    ggplot2::geom_point(aes(x = x,y = y,color = features),size = size, ...)+
    ggplot2::theme_void()+
    ggplot2::theme(plot.title = element_text(color = "black",hjust = 0.5,size = 30,face = "bold"))+
    ggplot2::labs(color = "Exp",title = features)+
    ggplot2::coord_cartesian(xlim = c(slice_coord_x_min,slice_coord_x_max),
                             ylim = c(slice_coord_y_min,slice_coord_y_max)) -> plot_result
  if(!is.null(cols)){
    plot_result = plot_result + ggplot2::scale_color_gradientn(colors = cols)
  }
  return(plot_result)
}


#' Title plot the point in spatial chip
#'
#' @param data the dataframe include x and y column represent the location in chip.
#' @param features which features you want to display.
#' @param size the ggplot2 point size 
#' @param cols the colors
#' @param ... 
#'
#' @return ggplot2 object
#' @export plot_chip
#'
#' @examples #
plot_chip <- function(data = data,features = features,size = 2, cols = c("white", "#e31a1c"),...){
  # 设置范围-------------------
  slice_coord_x_min = min(data$x)
  slice_coord_x_max = max(data$x)
  slice_coord_y_min = min(data$y)
  slice_coord_y_max = max(data$y)
  colnames(data)[grep(features,colnames(data))] = "name"
  ggplot2::ggplot(data)+
    ggplot2::geom_point(aes(x = x,y = y,color = name),size = size, ...)+
    ggplot2::theme_void()+
    ggplot2::theme(plot.title = element_text(color = "black",hjust = 0.5,size = 30,face = "bold"))+
    ggplot2::labs(color = "Exp",title = features)+
    ggplot2::coord_cartesian(xlim = c(slice_coord_x_min,slice_coord_x_max),
                             ylim = c(slice_coord_y_min,slice_coord_y_max)) -> plot_result

  if(!is.null(cols)){
    plot_result = plot_result + ggplot2::scale_color_gradientn(colors = cols)
  }
  return(plot_result)
}
