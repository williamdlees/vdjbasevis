########################################################################################################
#' Graphical output of alleles appeareance
#'
#' The \code{alleleAPP} function generates a graphical output of the alleles appeareance for each gene.
#'
#'
#' @param    gene_list            list of alleles appearce for each gene.
#' @param    max_alleles_bar      the number of maximum alleles per each bar plot. Defualt is 8 alleles.
#' @param    file                 file path to render the plot.
#'
#' @return
#'
#' A bar plot for each gene with the allele appearence.
#'
#'
#'
#' @export
alleleAPP <- function(gene_list, max_alleles_bar = 8, file = file.path(normalizePath(tempdir()),"allele_appeareance.pdf")) {


  plot_app <- function(alleles, appears, plot_name, nplots, max_a, max_p, p_dis, k, gene, long_names, singles = FALSE) {
    names_presents <- NULL
    for (a in alleles) {
      if (grepl('_',a,fixed=T)) {
        d <- paste0("*", (length(long_names)+1))
        names_presents <- c(names_presents, d)
        long_names <- c(long_names, paste0(gene, "*", a))
      } else {
        if (grepl("[A,T,C,G][0-9]+[A,C,G,T]", a)) {
          a <- sub("_", "\n", a)
        }
        names_presents <- c(names_presents, a)
      }
    }

    xx <- barplot(appears, names.arg = names_presents, cex.names = 0.9, ylim=c(0,max_a *1.08), main = plot_name, beside=FALSE, las=2)
    text(x = xx, y = appears, label = appears, pos = 3, cex = 0.8, col = "black")
    box()
    if(!singles) {
      axis(4, at = seq(0, max_a*1.08, max_a/(max_p/p_dis)),labels = seq(0, max_p*1.08, p_dis),  cex.axis = 1,las=1)
    }
    if ((nplots) %% 3 == 1) {
      mtext("Appearance in population",side = 2, line = 2.5, cex = 0.9)
    } else if ((nplots) %% 3 == 0) {
      if(!singles) {
        mtext("Appearance(%)",side = 4, line = 2.5, cex = 0.9)
      }
    }
    if(!singles) {
      if ((nplots - 1) %% 9 > 5 | k + 3 >= length(data_merge)) {
        mtext("Alleles",side = 1, line = 3.5, cex = 1)
      }
    } else {
      if ((nplots - 1) %% 6 >= 3 | k + 3 >= length(data_merge)) {
        mtext("Alleles",side = 1, line = 7.5, cex = 1)
      }
    }

    if(nplots %% 9 == 0) {
      plot(0:40,0:40,type="n",xlab = "", axes = FALSE)
      i <- 1
      for (n in long_names) {
        d <-paste0("*", as.character(i))
        long_name <- paste0(d, " - ", n)
        text(as.integer((i - 1)/8) * 8 + 2.5, 37.5 - (as.integer((i - 1) %% 8)*5), long_name, cex = 0.65)
        i <- i + 1
      }
      long_names <- NULL
    }
    return(long_names)
  }


  genes <- names(gene_list)
  gene_id <- 1
  plot_id <- 1
  long_names <- NULL # indicates the length of allele annotation
  singles_genes <- grepl('^Single',genes,perl=T)

  pdf(file,height=9,width=12)
  layout(matrix(c(1:10,10,10), 4, 3, byrow = TRUE),heights = c(1,1,1,0.55))
  par(mar=c(3.5, 3, 1., 2) + 0.1,oma=c(0,1.5,0,2))
  for (alleles in gene_list[!singles_genes]) {
    sum_a <- sum(alleles$Appeared) # total appeareance
    max_a <- max(alleles$Appeared) # maximun appeareance
    max_p <- as.integer((max_a * 100)/alleles$Total[1]) # maximun frequency
    p_dis <- 20 # the bar y axis "jump" value
    if (max_p > 60) {
      p_dis <- 20
    } else if (max_p > 40) {
      p_dis <- 10
    } else if (max_p > 15) {
      p_dis <- 5
    } else if (max_p > 7) {
      p_dis <- 2
    } else {
      p_dis <- 1
    }

    num_of_plots <- ceiling((length(alleles$Allele)-1)/max_alleles_bar)

    if (num_of_plots > 1) {
      cut_a <- 1
      while(cut_a < length(alleles$Allele)) {
        plot_num <- ceiling(cut_a/max_alleles_bar)
        plot_name <- paste0(genes[gene_id], " [", plot_num, " of ", num_of_plots,"]")
        long_names <- plot_app(alleles$Allele[cut_a : min(length(alleles$Appeared), (cut_a + max_alleles_bar - 1))],
                                   alleles$Appeared[cut_a : min(length(alleles$Appeared), (cut_a + max_alleles_bar - 1))],
                                   plot_name, plot_id, max_a, max_p, p_dis, gene_id - (num_of_plots - plot_num), genes[gene_id], long_names)
        plot_id <- plot_id + 1
        cut_a <- cut_a + max_alleles_bar
      }
    } else {
      long_names <- plot_app(alleles$Allele, alleles$Appeared, genes[gene_id], plot_id, max_a, max_p, p_dis, gene_id, genes[gene_id], long_names)
      plot_id <- plot_id + 1
    }
    gene_id <- gene_id + 1
  }

  temp_index <- plot_id
  if(!is.null(long_names)) {
    while (temp_index %% 9 != 1) {
      plot(0:4,0:4,type="n",xlab = "", axes = FALSE)
      temp_index <- temp_index+1
    }
    plot(0:40,0:40,type="n",xlab = "", axes = FALSE)
    i <- 1
    for (n in long_names) {
      d <-paste0("*", i)
      long_name <- paste0(d, " - ", n)
      text(as.integer((i - 1)/8) * 8 + 2.5, 37.5 - (as.integer((i - 1) %% 8)*5), long_name, cex = 0.65)
      i <- i + 1
    }
    long_names <- NULL
  }

  if (any(singles_genes)) {
    plot_id <- 1
    par(mfrow = c(2,3), mar=c(10, 3, 1, 0.5) + 0.1,oma=c(0,1.5,0,2))
    for (alleles in gene_list[singles_genes]) {
      gene_id <- length(data_merge) - ceiling(length(alleles$Allele)/max_alleles_bar)
      max_a <- max(alleles$Appeared)
      cut_a <- 1
      plot_name <- "Single allele genes"
      while(cut_a < length(alleles$Allele)) {
        plot_app(alleles$Allele[cut_a : min(length(alleles$Appeared), (cut_a + max_alleles_bar - 1))],
                     alleles$Appeared[cut_a : min(length(alleles$Appeared), (cut_a + max_alleles_bar - 1))],
                     plot_name, plot_id, max_a, max_p, p_dis, gene_id - num_of_plots, "", long_names, singles = TRUE)
        cut_a <- cut_a + max_alleles_bar
        plot_id <- plot_id + 1
        gene_id <- gene_id + 1
      }
    }
  }
  dev.off()
}
