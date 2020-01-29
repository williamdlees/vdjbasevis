########################################################################################################
#' Graphical output of gene usage in given population
#'
#' The \code{geneUsage} function generates a graphical output of the gene usage in a given population.
#'
#'
#' @param    gene_segment         a data frame of gene usage in a given population. See details.
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
#'   \item \code{'FREQ'}:     Gene usage of each individual, comma delimited.
#' }
#'
#' @examples
#'  gene_segment <- data.frame(GENE = c("V1-2",'V3-3','D2-8','D3-16','J4','J6'), FREQ = rep("0.2,0.1,0.2,0.4,0.5,0.6,0.7",6))
#'  # plotting with base R boxplot
#'  p <- geneUsage(gene_segment, plot_style = "base")
#'  cowplot::ggdraw(p)
#'  # plotting with ggplot
#'  p <- geneUsage(gene_segment, plot_style = "ggplot")
#'  # plotting with plotly
#'  p <- geneUsage(gene_segment, plot_style = "plotly")
#' @export
geneUsage <- function(gene_segment, chain = c("IGH", "IGK", "IGL"), plot_style = c("base","ggplot","plotly"), genePlotOrder = c("D",'J','V')){

  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)

  if (missing(chain)) {
    plot_style = "base"
  }
  plot_style <- match.arg(plot_style)

  plot_usage_ggplot <- function(sub, title, fam_col){
    if(is.factor(sub$GENE)) sub$GENE <- droplevels(sub$GENE)
    sub$FAM <- factor(sub$FAM)
    ggplot(sub, aes_string(x = "GENE", y = "FREQ")) + geom_boxplot(outlier.size=NA,outlier.shape=NA) +
      geom_point(aes_string(fill="FAM"), colour='black',position = 'jitter',shape=21,size=1.75,stroke =0.15) +
      scale_fill_manual(name='', values = fam_col, drop=F) +
      ylab('Frequency') + xlab('Gene') + labs(title = title) +
      theme(axis.text.y = element_text(size=16), axis.title = element_text(size=16),
            axis.text.x = element_text(size=14,angle = 90, hjust = 1,vjust=0.5), legend.position = 'none',
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
  }

  plot_usage_base <- function(sub, title){
    if(is.factor(sub$GENE)) sub$GENE <- droplevels(sub$GENE)
    boxplot(FREQ ~ GENE, data = sub, lwd = 2, ylab = 'Frequency', xlab = "Gene", outline=FALSE, main = title,las = 2)
    stripchart(FREQ ~ GENE, data = sub, vertical = TRUE,
               method = "jitter", add = TRUE, pch = 20, col = unique(sub, by = c("GENE","FAM_COL"))$FAM_COL)
  }

  plot_usage_plotly <- function(sub, title, fam_col){
    if(is.factor(sub$GENE)) sub$GENE <- droplevels(sub$GENE)
    sub %>%
      plotly::plot_ly(colors = c("0"="black",fam_col)) %>%
      plotly::add_trace(x = ~as.numeric(GENE),y = ~FREQ, color = "0", fillcolor="transparent", type = "box", showlegend = FALSE) %>%
      plotly::add_markers(x = ~jitter(as.numeric(GENE)), inherit = F , y = ~FREQ, color = ~FAM,
                  marker = list(size = 7),
                  showlegend = FALSE, text = ~FREQ, hoverinfo = 'text') %>%
      plotly::layout(xaxis = list(title = "Gene", tickvals = 1:length(unique(sub$GENE)), ticktext = unique(sub$GENE)),
                     yaxis = list(title = "Frequency"))
  }

  # familiy colors
  fam_col <- setNames(c('deepskyblue','darkorange','forestgreen','red','mediumslateblue','saddlebrown','pink'),"1":"7")
  # split the frequency to rows
  gene_segment <- splitstackshape::cSplit(gene_segment, "FREQ", sep = ",", direction = "long", fixed = T, type.convert = F)
  gene_segment[,"FREQ" := as.numeric(gene_segment$FREQ)]
  # chain in gene name
  chain_gene <- grepl(paste0('^',chain),gene_segment$GENE)
  # get the gene familiy
  nth_fam <- ifelse(chain_gene, 5, 2)
  gene_segment[,"FAM":=substring(gene_segment$GENE, nth_fam, nth_fam)]
  gene_segment[,"FAM_COL":=fam_col[gene_segment$FAM]]
  # get the different segments
  nth_seg <- ifelse(chain_gene, 4, 1)
  gene_segment[,"SEGMENT":=substring(gene_segment$GENE, nth_seg, nth_seg)]
  gene_segment <- gene_segment[order(match(GENE, GENE.loc[chain]))]
  segment_order <- genePlotOrder[genePlotOrder %in% unique(gene_segment$SEGMENT)]
  n_plots <- length(unique(gene_segment$SEGMENT))
  layout_p <- if(n_plots == 3) c(1,3,2,3) else if(n_plots == 2) c(1,2) else 1
  # plot layout
  layout.matrix <- matrix(layout_p, nrow = ifelse(n_plots>1,n_plots-1,1))

  if(plot_style=="base"){
    graphics::layout(mat = layout.matrix
                     #,heights = c(3, 0.5, 0.3) # Heights of the three rows
    )
    invisible(lapply(segment_order, function(x){
      plot_usage_base(gene_segment[gene_segment$SEGMENT == x,],
                      paste0(chain,x))
    }))
    p1.base <- recordPlot()
    invisible(dev.off())

  }else{
    if(plot_style=="ggplot"){
      plotList <- lapply(segment_order, function(x){
        plot_usage_ggplot(gene_segment[gene_segment$SEGMENT == x,],
                        paste0(chain,x), fam_col)
      })
      top <- plotList[1:ifelse(n_plots>1,n_plots-1,1)]
      bottom <- ifelse(n_plots>1,plotList[n_plots],NA)
      p1.base <- cowplot::plot_grid(plotlist = top, nrow = 1, align = "hv")
      if(!is.na(bottom)) p1.base <- cowplot::plot_grid(p1.base, plotlist = bottom, nrow =  2)
    }else{
      plotList <- lapply(unique(gene_segment$SEGMENT), function(x){
        plot_usage_plotly(gene_segment[gene_segment$SEGMENT == x,],
                        paste0(chain,x), fam_col)
      })
      nrows <- ifelse(n_plots>1,n_plots-1,1)
      p1.base <- plotly::subplot(plotList, nrows = nrows, margin = 0.1)
      if(nrows)  p1.base <- p1.base %>% plotly::layout(xaxis3 = list(domain=list(x=c(0,0.5),y=c(0,0.5))))
    }
  }
  return(p1.base)
}
