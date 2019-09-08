#' @include vdjbaseVis.R
#' @include internal_functions.R
NULL

########################################################################################################
#' Graphical output of alleles division by genotypes
#'
#' The \code{genoHeatmap_html} function generates an interactive graphical output of the alleles per gene in multiple samples.
#'
#'
#' @param    geno_table           genoytpe summary table. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    gene_sort            if by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#' @param    lk_cutoff            the lK cutoff value to be considerd low for texture layer. Defualt is lK<1.
#' @param    mark_low_lk          if TRUE, a texture is add for low lK values. Defualt is TRUE.
#' @param    html                 if TRUE, an interactive html visualization is produced. Defualt is FALSE.
#'
#' @return
#'
#' A heat-map visualization of the genotype inference for multiple samples.
#'
#' @details
#'
#' A \code{data.frame} created by \code{inferGenotypeBaysian}.
#'
#'
#' @export
genoHeatmap_html <- function(geno_table, chain = c("IGH", "IGK", "IGL"), gene_sort = "position", removeIGH = TRUE, lk_cutoff = 1, mark_low_lk = TRUE, html = FALSE, n_line = 4, line_length=60, pseudo_genes = FALSE, ORF_genes = FALSE) {


  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)
  lk_cutoff = as.numeric(lk_cutoff)



  samples <- unique(geno_table$SUBJECT)

  geno_db <- geno_table %>% select(SUBJECT,GENE,GENOTYPED_ALLELES,K_DIFF)
  names(geno_db)[3:4] <- c("ALLELES", "K")
  geno_db <- tidyr::separate_rows(geno_db, "ALLELES", sep = ",")
  genes <- unique(geno_db$GENE)
  geno_db$ALLELES <- gsub("Deletion","Del",geno_db$ALLELES)
  for(s in samples){
    sub <- geno_db[geno_db$SUBJECT==s,]
    # get missing genes
    id <- which(!genes %in% sub$GENE)
    if(length(id) != 0) geno_db <- bind_rows(geno_db,data.frame(SUBJECT = s, GENE = genes[id], ALLELES = 'Unk',
                                                                K = NA_integer_, stringsAsFactors = F))
  }

  color_pes_orf <- c()
  if(pseudo_genes){
    color_pes_orf <- c(grep("V",PSEUDO[[chain]],value = T),color_pes_orf)
  }
  if(ORF_genes){
    color_pes_orf <- c(unique(grep("OR|NL", geno_db$GENE,value = T)),color_pes_orf)
  }



  geno_db <- sortDFByGene(DATA = geno_db, chain = chain, method = gene_sort, removeIGH = removeIGH, geno = T,
                          peseudo_remove = pseudo_genes, ORF_remove = ORF_genes)



  geno_db$K[grep("Del",geno_db$ALLELES)] <- NA_integer_

  gene_num <- round(length(unique(geno_db$GENE))/3)

  geno_db_texture <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("SUBJECT","GENE","ALLELES",'K','points','yend','x','xend'))
  loc <- 1:length(levels(droplevels(geno_db$GENE)))
  names(loc) <- rev(levels(droplevels(geno_db$GENE)))
  for (i in 1:nrow(geno_db)) {
    if (geno_db$K[i] < lk_cutoff && !geno_db$ALLELES[i] %in% c("Unk", "Del", "NR")) {
      tmp_point <- geno_db[i, ] %>% slice(rep(1, each = n_line)) %>%
        dplyr::mutate(points = seq(0, 0.9, length.out = n_line),
                            yend = seq(0, 0.9, length.out = n_line) + 0.1, x = loc[as.character(GENE)] - 0.49,
               xend = loc[as.character(GENE)] + 0.49)
      geno_db_texture <- bind_rows(geno_db_texture, tmp_point)
    }
  }



  heatmap.df <- as.data.frame(geno_db %>% dplyr::group_by(SUBJECT, GENE) %>% dplyr::mutate(n = dplyr::n()))
  heatmap.df$freq <- ifelse(heatmap.df$n == 2, 0.5, ifelse(heatmap.df$n != 1, 0.25, 1))
  heatmap.df$GENE <- factor(heatmap.df$GENE, levels = gsub("IG[H|K|L]", "", GENE.loc[[chain]]))
  gene_loc <- 1:length(unique(heatmap.df$GENE)[order(match(unique(heatmap.df$GENE), levels(heatmap.df$GENE)))])
  names(gene_loc) <- unique(heatmap.df$GENE)[order(match(unique(heatmap.df$GENE), levels(heatmap.df$GENE)))]
  heatmap.df$GENE_LOC <- gene_loc[as.character(heatmap.df$GENE)]
  heatmap.df$title <- "Genotype"



  if (length(grep("[0-9][0-9]_[0-9][0-9]", heatmap.df$ALLELES)) != 0) {
    non_reliable_alleles_text <- nonReliableAllelesText_V2(heatmap.df[grep("[0-9][0-9]_[0-9][0-9]", heatmap.df$ALLELES), ],map = T)
  } else {
    non_reliable_alleles_text <- c()
  }

  heatmap.df$ALLELE_TEXT <- sapply(1:nrow(heatmap.df),function(i){
    if(grepl("[0-9][0-9]_[0-9][0-9]", heatmap.df$ALLELES[i])){
      unique(non_reliable_alleles_text$text_bottom[heatmap.df$GENE[i]==non_reliable_alleles_text$GENE & grepl(heatmap.df$ALLELES[i],non_reliable_alleles_text$text_bottom)])
    }else{
      heatmap.df$ALLELES[i]
    }
  })

  heatmap.df$ALLELES[grep("[0-9][0-9]_[0-9][0-9]", heatmap.df$ALLELES)] <- "NRA"
  allele_palette <- allelePalette(heatmap.df$ALLELES)
  non_reliable_alleles_text$ALLELES <- factor(non_reliable_alleles_text$ALLELES, levels = allele_palette$AlleleCol)
  non_reliable_alleles_text$GENE_LOC <- gene_loc[as.character(non_reliable_alleles_text$GENE)]

  if(any(grepl('^[0-9]+[_][A-Z][0-9]+[A-Z]',heatmap.df$ALLELES))){
  #check novel allele count
  novel <- data.frame(Novel=grep('^[0-9]+[_][A-Z][0-9]+[A-Z]',heatmap.df$ALLELES,value=T),
                      Base = sapply(grep('[A-Z][0-9]+[A-Z]',heatmap.df$ALLELES,value=T),
                                    function(x) strsplit(x,'[_]')[[1]][1])) %>%
    dplyr::distinct() %>% dplyr::group_by(Base) %>% dplyr::mutate(n = dplyr::n())
  # if any base is larger then 6 change to dagger and number
  dagger = "\u005E"
  if(any(novel$n>6)){
    id <- grep('^[0-9]+[_][A-Z][0-9]+[A-Z]',names(allele_palette$transper))
    allele_palette$transper[id] <- 1
    new_allele <- paste0(dagger,1:length(id),'-',allele_palette$AlleleCol[id])
    names(new_allele) <-allele_palette$AlleleCol[id]
    novel_allele_text <- novelAlleleAnnotation(heatmap.df[grep(paste0("^",names(new_allele),collapse = "|"), heatmap.df$ALLELES), ],
                                               new_label = new_allele)
    for(i in 1:length(new_allele)){
      allele <- names(new_allele)[i]
      heatmap.df$ALLELES[grep(allele,heatmap.df$ALLELES)] <- new_allele[i]
      heatmap.df$ALLELE_TEXT[grep(allele,heatmap.df$ALLELE_TEXT)] <- new_allele[i]
    }
    allele_palette$AlleleCol[id] <- new_allele
    names(allele_palette$transper)[id] <- new_allele
    novel <- T
  }else{
    novel <- F
    novel_allele_text <- c()
  }}else{
      novel <- F
      novel_allele_text <- c()
      }


  sub_lev <- F
  if (length(levels(geno_table$SUBJECT)) != 0) {
    sub_lev <- T
    heatmap.df$SUBJECT <- factor(heatmap.df$SUBJECT, levels = levels(geno_table$SUBJECT))
    non_reliable_alleles_text$SUBJECT <- factor(non_reliable_alleles_text$SUBJECT, levels = levels(geno_table$SUBJECT))
  }

  heatmap.df$ALLELES <- factor(heatmap.df$ALLELES, levels = allele_palette$AlleleCol)
  heatmap.df$text2 <- paste0("Individual: ",heatmap.df$SUBJECT,
                            '<br />', "Gene: ", heatmap.df$GENE, '<br />',"Allele: ", heatmap.df$ALLELE_TEXT,
                            '<br />',"lK: ", round(as.numeric(heatmap.df$K),4))

  col <- ifelse(c(1:(length(samples) * 4))%%3 == 0, "black", "transparent")
  col_x <- ifelse(names(gene_loc) %in% gsub(chain,"",color_pes_orf), "brown", "black")

  p <- ggplot() + geom_col(data = heatmap.df, mapping = aes_string(x = "GENE_LOC", y = "freq", fill = "ALLELES", text = "text2"),
                           position = "fill", width = 0.95, na.rm = T) +
    scale_fill_manual(values = alpha(names(allele_palette$AlleleCol), allele_palette$transper), name = "Alleles", drop = FALSE) +
    facet_grid(SUBJECT ~ title, as.table = FALSE, switch = "y")

  if(html){
    p <- p + scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0,0),
                                breaks = 1:length(gene_loc),
                                labels = names(gene_loc)) +
      theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1, size = 10, colour = col_x),
            strip.text.y = element_text(angle = 0, size = 10) ,strip.placement = "outside",
            strip.background = element_blank(),
            axis.ticks.y = element_line(colour = col), axis.text.y = element_blank(), strip.background.y = element_blank(),
            legend.position = "bottom",legend.justification = "center", panel.spacing.y = unit(0.9, "pt")) +
      labs(y = "", x = "Gene") + guides(fill = guide_legend(nrow = ceiling(length(allele_palette$AlleleCol)/gene_num),order = 1,
                                                            override.aes = list(color = "#DCDCDC")))
  }else{
    p <- p + scale_y_continuous(expand = c(0, 0), sec.axis =  dup_axis(name = "")) + scale_x_continuous(expand = c(0,0),
                                breaks = 1:length(gene_loc),
                                labels = names(gene_loc),
                                sec.axis = dup_axis(name = "")) +
      theme(axis.text.x.top = element_text(angle = 90,vjust = 0.5, hjust = 0, size = 10, colour = col_x),
            axis.text.x.bottom = element_text(angle = 90,vjust = 0.5, hjust = 1, size = 10, colour = col_x),
            strip.text.x = element_blank(),
            strip.text.y = element_text(angle = 180, size = 10) ,strip.placement = "outside",
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size=4, fill=NA),
            axis.ticks.y.left = element_line(colour = col), axis.ticks.y.right = element_blank(), axis.text.y = element_blank(), strip.background.y = element_blank(),
            legend.position = "bottom",legend.justification = "center", panel.spacing.y = unit(0.9, "pt")) +
      labs(y = "", x = "Gene") + guides(fill = guide_legend(nrow = ceiling(length(allele_palette$AlleleCol)/gene_num),order = 1,
                                                            override.aes = list(color = "#DCDCDC")))
  }

  legend_rows = ceiling(length(allele_palette$AlleleCol)/gene_num)

  add_geom = F
  if (mark_low_lk & nrow(geno_db_texture) != 0) {

    geno_db_texture <- geno_db_texture[!duplicated(geno_db_texture[, c("GENE", "ALLELES", "K", "points", "SUBJECT")]), ]
    geno_db_texture$col <- paste0("lk<",lk_cutoff)
    if(sub_lev) geno_db_texture$SUBJECT <- factor(geno_db_texture$SUBJECT, levels = levels(geno_table$SUBJECT))
    geno_db_texture$ALLELES[grep("[0-9][0-9]_[0-9][0-9]", geno_db_texture$ALLELES)] <- "NRA"
    geno_db_texture$ALLELES <- factor(geno_db_texture$ALLELES, levels = allele_palette$AlleleCol)

    geno_db_texture <- merge(geno_db_texture,heatmap.df,by = names(heatmap.df)[1:4])

    # Get Allele legend
    gt1 = ggplotGrob(p)

    if(html){
      p <- p + geom_segment(data = geno_db_texture,
                            mapping = aes_string(x = "x", xend = "xend", y = "points", yend = "yend", linetype = "col",text="text2"),
                            colour = "white")
      add_geom = T
    }else{
      p <- p + geom_segment(data = geno_db_texture,
                            mapping = aes_string(x = "x", xend = "xend", y = "points", yend = "yend", color = "col",text="text2"),
                            colour = "white")

      gt2 = getDigLegend(unique(geno_db_texture$col))

      # Get lK legend
      leg1 = gtable::gtable_filter(gt1, "guide-box")
      leg2 = gtable::gtable_filter(gt2, "guide-box")
      # Combine the legends
      leg <- cbind(leg1[["grobs"]][[1]], leg2[["grobs"]][[1]], size = "first")
      # Insert legend into g1 (or g2)
      gt1$grobs[gt1$layout$name == "guide-box"][[1]] <- leg
      gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[3, c("t", "b")] <- gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[1, c("t", "b")]
      gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[3, c("l", "r")] <- gt1$grobs[gt1$layout$name == "guide-box"][[1]]$layout[1, c("l", "r")] +
        3
      legend <- get_legend(gt1)}
  } else legend <- get_legend(p)

  short_reads = F
  if(is.data.frame(non_reliable_alleles_text) | novel){
  unique_text <- bind_rows(non_reliable_alleles_text, novel_allele_text) %>% group_by(GENE, SUBJECT) %>%
    dplyr::mutate(pos_f = ifelse(n == 1, 0.5,
                        ifelse(n == 2, seq(0.25,1,by = 0.5)[1:max(dplyr::row_number())],
                               ifelse(n == 3, seq(0.165,1,by = 0.33)[1:max(dplyr::row_number())],
                                      seq(0.125,1,by = 0.25)[1:max(dplyr::row_number())]))))

  unique_text$text2 <- paste0("Individual: ",unique_text$SUBJECT,
                                            '<br />', "Gene: ", unique_text$GENE,
                                            '<br />',"Allele: ", unique_text$text_bottom,
                                            '<br />',"lK: ", round(as.numeric(unique_text$K),4))


  bottom_annot <- c()
  if (is.data.frame(non_reliable_alleles_text)) {
    p <- p + geom_text(data = unique_text[unique_text$ALLELES=="NRA",],
                       aes_string(label = "text", x = "GENE_LOC", y = "pos_f",text = "text2"), angle = 0,
                       size = unique_text$size[unique_text$ALLELES=="NRA"])

    bottom_annot <- unique(non_reliable_alleles_text$text_bottom)
    names(bottom_annot) <- unique(non_reliable_alleles_text$text)
    short_reads = T
  }
  ## multiple novel alleles annot
  if (novel) {
    p <- p + geom_text(data = unique_text[unique_text$ALLELES!="NRA",],
                       aes_string(label = "text", x = "GENE_LOC", y = "pos_f", text = "text2"), angle = 0,
                       size = unique_text$size[unique_text$ALLELES!="NRA"])
  }}

  if(html){
    options(warn = -1)
    pl <- ggplotly(p + xlab(""), tooltip = "text2") %>%
      plotly::layout(xaxis = list(title = paste0(c(rep("&nbsp;", 3),
                                                   "Gene", rep("&nbsp;", 3), rep("\n&nbsp;", 1)),
                                                 collapse = "")))

    # add bg clor to legend (seashell2)
    pl$x$layout$legend$bgcolor = "#eee5de"

    # switch sample names to left yaxis
    smaples_id <- grep("Genotype|Alleles|Gene", pl$x$layout$annotations, invert = T)
    # add margin on the left
    pl$x$layout$margin$l = 40
    pl$x$layout$margin$b = 150
    pl$x$layout$margin$l = 100
    samples_names <- c()
    for(i in smaples_id){
      # move to the left side
      pl$x$layout$annotations[[i]]$x = -0.005
      # anchor to the right
      pl$x$layout$annotations[[i]]$xanchor = "right"
      samples_names <- c(samples_names, pl$x$layout$annotations[[i]]$text)
    }

    pl$x$layout$xaxis$ticktext <- sapply(pl$x$layout$xaxis$ticktext, function(x){
      ifelse(x %in% gsub(chain,"",color_pes_orf), paste0("[",x,"]"),x)
    },USE.NAMES = F)

    names(samples_names) <- 1:length(samples_names)
    # get the y-axis
    yaxis = grep('yaxis',names(pl$x$layout))
    # remove ticks
    for(i in yaxis){
      pl$x$layout[[i]]$tickcolor = "transparent"
    }


    text_for_hovertext <- function(text, annot, lk_annot, nra=F, allele_anot=c()) {
      labels <- text
      for (i in 1:length(text)) {
        label <- text[i]
        gene <- strsplit(strsplit(label, "Gene: ")[[1]][2], '[<]')[[1]][1]
        gene_text <- annot[as.numeric(gene)]
        if(nra){
          allele_nra <- strsplit(strsplit(label,'[<]')[[1]][[1]]," ")[[1]][2]
          allele <- paste0("Allele: ", allele_anot[allele_nra])
          lk = round(unique(lk_annot$K[lk_annot$GENE==gene_text]),4)
        }else{
          allele <- strsplit(grep(paste0(gene,'[<]'),text,value=T), "<")[[1]][1]
          lk = round(unique(lk_annot$K[lk_annot$GENE==gene_text]),4)
        }

        labels[i] <- paste0("Gene: ", gene_text, "<br />", allele, "<br />", "lK: ",lk)
      }
      return(labels)
    }

    # sort annotations
    # title
    pl$x$layout$annotations[[1]]$text = ""
    # legend name
    pl$x$layout$annotations[[grep("Alleles", pl$x$layout$annotations)]]$text = "Alleles"


    # Sort legend after K line addition
    # get the duplicated markers
    unique_legend_markers <- sapply(pl$x$data, function(x) gsub('[(]','',strsplit(x$name,',')[[1]][1]))
    group_legend <- sapply(pl$x$data, function(x) gsub('[)]','',strsplit(x$name,',')[[1]][2]))
    dup <- which(duplicated(unique_legend_markers))
    gene_loc_id <- names(gene_loc)
    names(gene_loc_id) <- gene_loc
    novel_allele_text_l <- novel_allele_text$text_bottom
    names(novel_allele_text_l) <- novel_allele_text$text
    for(i in 1:length(unique_legend_markers)){
      if(i %in% dup){
        # Turn off duplicated
        pl$x$data[[i]]$showlegend = FALSE
      }else{
        # sort the names of the unique markers
        pl$x$data[[i]]$name = unique_legend_markers[[i]]
      }
      if(!grepl('[<]|NRA',unique_legend_markers[[i]])){
        pl$x$data[[i]]$marker$line$color = "black"
        pl$x$data[[i]]$marker$line$width = 0.1
        if(any(grepl('[*]',pl$x$data[[i]]$text))){
          pl$x$data[[i]]$type <-  suppressWarnings("text")
          pl$x$data[[i]]$hoveron <- "fill"
        }else{
          if(!is.null(pl$x$data[[i]]$hovertext)){
            if(grepl("\\^",pl$x$data[[i]]$text, fixed = T)){
            pl$x$data[[i]]$type <-  suppressWarnings("text")
            pl$x$data[[i]]$hoveron <- "fill"
            pl$x$data[[i]]$hoverinfo <- "skip"
            }
          }
        }

      }else{
        # sort hovertext for K line
        pl$x$data[[i]]$hoverinfo <- "skip"
        if(grepl('[<]',unique_legend_markers[[i]])) pl$x$data[[i]]$name = unique_legend_markers[[i]]
        else{
          pl$x$data[[i]]$name = unique_legend_markers[[i]]
        }
      }

      # switch all to single group
      pl$x$data[[i]]$legendgroup = unique_legend_markers[i]
      pl$height = max(length(samples) * 35,400)
      pl$width = "150%"

    }

    return(pl)

  }else{
    p2 <- plot_grid(p + theme(legend.position = "none",
                              plot.background = element_rect(fill = "transparent", colour = NA),
                              panel.background = element_rect(fill = "transparent", colour = NA),
                              axis.line = element_line(color="black")),
                    legend, rel_heights = c(0.9,0.1), ncol = 1)

    if(short_reads){
      bottom_annot <- unique(non_reliable_alleles_text$text_bottom)
      # Create text for annotating the short reads labels
      annot <- splitlines(bottom_annot, line_length)
      annot <- annot[!is.na(dplyr::na_if(annot,""))]

      gp <- cowplot::ggdraw(cowplot::add_sub(p2,annot))
      #lab <- grid::textGrob(paste0(annot,collapse = '\n')
      #                      , x = unit(.1, "npc"), just = c("left"), gp = grid::gpar(fontsize = 10, col = "black"))
      short_reads_rows = length(annot)
      #gp <- ggplotGrob(p2)
      # Add a row below the 2nd from the bottom
      #gp <- gtable::gtable_add_rows(gp, unit(1, "grobheight", lab), -2)

      # Add 'lab' grob to that row, under the plot panel
      #gp <- gtable::gtable_add_grob(gp, lab, t = -2, l = gp$layout[grep('panel',gp$layout$name)[1],]$l)

      return(list(p=plot_grid(gp),legend_rows = legend_rows, short_reads_rows = short_reads_rows))
    }else{
      return(list(p=plot_grid(p2),legend_rows = legend_rows, short_reads_rows = 0))
    }
  }
}
