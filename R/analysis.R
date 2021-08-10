#' Title To integrate the scRNAseq and ST datasets, MIA method can analysis proceeds by first delineating sets of cell type-specific and tissue region-specific genes and then determining whether their overlap is higher (enrichment) or lower (depletion) than expected by chance.
#'
#' @param sc_object the seurat scRNA object
#' @param sp_object the seurat scRNA object or spatial object
#' @param sc_diff the different expression result of seurat analysis
#' @param sp_diff the different expression result of seurat analysis
#' @importFrom stats phyper
#'
#' @return data
#' @export MIA_analysis
#'
#' @examples # d = MIA_analysis(sc_object = stereo_smn, sp_object = stereo, sc_diff = sc_diff, sp_diff = sp_diff)
MIA_analysis <- function(sc_object, sp_object, sc_diff, sp_diff){
  dt = data.frame()
  sum_gene = unique(union(rownames(sc_object),rownames(sp_object)))
  for(i in levels(sc_diff$cluster)){
    sc_gene = sc_diff[sc_diff$cluster == i,]$gene
    for(j in levels(sp_diff$cluster)){
      sp_gene = sp_diff[sp_diff$cluster == j,]$gene
      in_gene = intersect(sc_gene,sp_gene)
      p = -log10(phyper(length(in_gene),length(sc_gene),length(sum_gene)-length(sc_gene),length(sp_gene),lower.tail = F))
      dt[i,j] = p
    }
  }
  return(dt)
}
