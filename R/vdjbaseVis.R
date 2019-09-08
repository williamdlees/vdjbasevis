# vdjbaseVis package documentation and import directives

#' The vdjbaseVis package
#'
#' The \code{vdjbaseVis} package provides multiple visualizations for genotypes.
#'
#'
#' @section  Genotype visualization:
#' \code{vdjbaseVis} provides tools to infer haplotypes based on given anchor genes,
#' deletion detection based on relative gene usage, pooling v genes, and a single anchor gene.
#'
#' \itemize{
#'   \item  \link{genoHeatmap}:              Genotype comparison for multiple smaples in heatmap format.
#'   \item  \link{multipleGenoytpe}:         Genotype visualization for single and multiple smaples.
#' }
#'
#' @name     vdjbaseVis
#' @docType  package
#' @references
#' \enumerate{
#'  }
#'
#' @import   ggplot2
#' @import   graphics
#' @import   methods
#' @import   utils
#' @import   dendextend
#' @importFrom  cowplot      get_legend plot_grid ggdraw draw_label background_grid
#' @importFrom  gridExtra    arrangeGrob
#' @importFrom  plotly       ggplotly subplot
#' @importFrom  dplyr        do n desc funs %>% distinct
#'                           as_data_frame data_frame data_frame_
#'                           bind_cols bind_rows combine rowwise slice
#'                           filter filter_ select select_ arrange arrange_
#'                           group_by group_by_ ungroup
#'                           mutate mutate_ summarize summarize_
#'                           mutate_at summarize_at count_ count na_if
#'                           rename rename_ transmute transmute_ pull ungroup row_number
#' @importFrom  data.table   := rbindlist data.table .N setDT CJ setorderv setkey
#' @importFrom  reshape2     melt
#' @importFrom  gtools       ddirichlet
#' @importFrom  stats        hclust as.dendrogram as.dist binom.test p.adjust setNames weighted.mean
#' @importFrom  ggdendro     dendro_data segment
#' @importFrom  htmlwidgets  saveWidget
#' @importFrom  gtable       gtable_filter
#' @importFrom  grDevices    dev.off pdf colorRampPalette embedFonts recordPlot
#' @importFrom  alakazam     getGene
#' @importFrom  rlang        .data
#' @importFrom  tigger       sortAlleles
#' @importFrom  RColorBrewer brewer.pal
#' @importFrom  tidyr        separate_rows drop_na
#' @importFrom  stringi      stri_detect_regex stri_detect_fixed
#' @importFrom  grid         gpar textGrob
#' @importFrom  mltools      bin_data
#' @importFrom  splitstackshape cSplit
#' @importFrom  fastmatch    %fin%
NULL

