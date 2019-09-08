########################################################################################################
#' Graphical output of allele usage in given population
#'
#' The \code{alleleUsageBar_html} function generates an interactive graphical output of sum of unique allele that appeared in the given population.
#'
#'
#' @param    gene_segment         a data frame of allele usage count in a given population. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#'
#' @return
#'
#' An interactive stacked barplot visualization of the allele hetrouzygousity in a given population.
#'
#' @details
#'
#' A \code{data.frame} with the following columns.
#' \itemize{
#'   \item \code{'GENE'}:     The gene call
#'   \item \code{'COUNT'}:    Number of unique alleles that appeared in the given population.
#' }
#'
#' @example
#'  gene_segment <- data.frame(GENE = c("V1-2",'V3-3','D2-8','D3-16','J4','J6'), COUNT = c(2,9,2,2,1,2))
#'  alleleUsageBar_html(gene_segment)
#' @export
alleleUsageBar_html <- function(gene_segment, chain = c("IGH", "IGK", "IGL")){

  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)

  plot_usage <- function(sub, title){
    sub %>%
      plotly::plot_ly(
        #width = 500,
        #height = 1500,
        y = ~GENE,
        x = ~COUNT,
        name = " ",
        type = 'bar', orientation = 'h', legendgroup = "A",
        marker = list(color = 'peru') , showlegend = F) %>%
      plotly::layout(annotations = list( text = title, xref = "paper", yref = "paper",
                    yanchor = "bottom", xanchor = "center", align = "center", x = 0.5, y = 1, showarrow = FALSE),
                     xaxis = list(title = ""),hovermode = 'y')
  }

  # get the different segments
  nth <- ifelse(grepl(paste0('^',chain),gene_segment$GENE), 4, 1)
  setDT(gene_segment)[,SEGMENT:=substring(GENE, nth, nth)]

  # create the plots based on segments
  plotList <- lapply(unique(gene_segment$SEGMENT), function(x){

    plot_usage(gene_segment[gene_segment$SEGMENT == x,],
               paste0(chain,x))

  })

  # plot all segments
  return(plotly::subplot(plotList,nrows = 1, margin = 0.1))
}
