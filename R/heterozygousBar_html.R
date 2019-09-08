########################################################################################################
#' Graphical output of allele heterozygousity
#'
#' The \code{heterozygousBar_html} function generates an interactive graphical output of the alleles heterozygousity in a given population.
#'
#'
#' @param    gene_segment         a data frame of allele heterozygousity count in a given population. See details.
#' @param    chain                 the IG chain: IGH,IGK,IGL. Default is IGH.
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
#'   \item \code{'HM'}:       Count of the Individuals homozygous to the gene
#'   \item \code{'HT'}:       Count of the Individuals heterozygous to the gene
#' }
#'
#' @example
#'  gene_segment <- data.frame(GENE = c("V1-2",'V3-3','D2-8','D3-16','J4','J6'), HM = c(20,60,55,7,30,0) , HT = c(80,40,45,93,0,45))
#'  heterozygousBar_html(gene_segment)
#' @export
heterozygousBar_html <- function(gene_segment, chain = c("IGH", "IGK", "IGL")){

  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)

  plot_hetro <- function(sub, title, ...){
  sub %>%
    plotly::plot_ly(
      #width = 500,
      #height = 1500,
      y = ~GENE,
      x = ~FRAC_HT,
      text = ~paste0("<br />", HT," from ", TOTAL, 'subjects'),
      type = 'bar', orientation = 'h',
      legendgroup = "A",
      name = "Heterozygous",
      marker = list(color = 'indianred', line = list(color = 'indianred', width = 1.5)),...) %>%
    plotly::add_trace(
        x = ~FRAC_HM,
        text = ~paste0("<br />", HM," from ", TOTAL, 'subjects'),
        name = "Homozygous",
        legendgroup = "B",
        marker = list(color = 'skyblue', line = list(color = 'skyblue', width = 1.5)),...) %>%
    plotly::layout(annotations = list( text = title, xref = "paper", yref = "paper",
                       yanchor = "bottom", xanchor = "center", align = "center", x = 0.5, y = 1, showarrow = FALSE),
                   xaxis = list(title = ""),barmode = "stack",hovermode = 'y')
  }

  # create the fractions
  setDT(gene_segment)[, c("FRAC_HT","FRAC_HM", "TOTAL") := list((HT)/(HM+HT),(HM)/(HM+HT),(HM+HT)), by=GENE]
  # get the different segments
  nth <- ifelse(grepl(paste0('^',chain),gene_segment$GENE), 4, 1)
  gene_segment[,SEGMENT:=substring(GENE, nth, nth)]

  # create the plots based on segments
  plotList <- lapply(1:length(unique(gene_segment$SEGMENT)), function(x){

    plot_hetro(gene_segment[gene_segment$SEGMENT == unique(gene_segment$SEGMENT)[x],],
               paste0(chain,unique(gene_segment$SEGMENT)[x]),
               showlegend = ifelse(x>1,F, T))

  })

  # plot all segments
  return(plotly::subplot(plotList,nrows = 1, margin = 0.1))
}
