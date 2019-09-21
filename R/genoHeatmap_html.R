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
#'
#' @return
#'
#' An interactive heat-map visualization of the genotype inference for multiple samples.
#'
#' @details
#'
#' A \code{data.frame} created by \code{inferGenotypeBaysian}.
#'
#'
#' @export
genoHeatmap_html <- function(geno_table, chain = c("IGH", "IGK", "IGL"), gene_sort = "position", removeIGH = TRUE, lk_cutoff = 1, mark_low_lk = TRUE, n_line = 4, line_length=60, pseudo_genes = FALSE, ORF_genes = FALSE, file = file.path(normalizePath(tempdir()),"genotype_heatmap.html")) {
  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)
  lk_cutoff = as.numeric(lk_cutoff)
    # select columns
  geno_db <- geno_table[,c("SUBJECT", "GENE", "GENOTYPED_ALLELES", "K_DIFF", "Freq_by_Clone")]
  # rename the columns
  names(geno_db)[3:4] <- c("ALLELES", "K")
  # correct deletion annotations
  geno_db$ALLELES <- gsub("Deletion","Del",geno_db$ALLELES)
  # set data.table and correct missing Unk annotations and K
  geno_db <- setDT(geno_db)[CJ(SUBJECT = SUBJECT, GENE = GENE, unique=TRUE), on=.(SUBJECT, GENE)]
  geno_db[is.na(ALLELES) , c("ALLELES","K") := list("Unk", NA_integer_)]
  # set K value for deleted genes
  geno_db$K[grep("Del",geno_db$ALLELES)] <- NA_integer_
  # expand row, one allele per row
  geno_db <- splitstackshape::cSplit(geno_db, "ALLELES", sep = ",", direction = "long", fixed = T, type.convert = F)

  # add pseudo genes and orf to color base
  color_pes_orf <- c()
  if(pseudo_genes){
    color_pes_orf <- c(grep("V",PSEUDO[[chain]],value = T),color_pes_orf)
  }
  if(ORF_genes){
    color_pes_orf <- c(unique(grep("OR|NL", geno_db$GENE,value = T)),color_pes_orf)
  }

  # sort the data, remove pseudo and orf if needed
  geno_db <- sortDFByGene(DATA = geno_db, chain = chain, method = gene_sort, removeIGH = removeIGH, geno = T,
                          peseudo_remove = pseudo_genes, ORF_remove = ORF_genes)

  geno_db$GENE <- factor(geno_db$GENE, levels = gsub("IG[H|K|L]", "", GENE.loc[[chain]]))

  # rename genes to numbers
  gene_loc <- 1:length(unique(geno_db$GENE)[order(match(unique(geno_db$GENE), levels(geno_db$GENE)))])
  names(gene_loc) <- unique(geno_db$GENE)[order(match(unique(geno_db$GENE), levels(geno_db$GENE)))]
  geno_db$GENE_LOC <- gene_loc[as.character(geno_db$GENE)]

  ######sort the heatmap for plotting
  geno_db_m <- geno_db[, n:=  .N, by = list(SUBJECT,GENE)][] # count number of alleles for group
  geno_db_m$ALLELES_G <- geno_db_m$ALLELES # for grouping
  geno_db_m$text <- ''
  geno_db_m$text_bottom <- geno_db_m$ALLELES
  # change ambiguous alleles call
  id_nra <- grepl("[0-9][0-9]_[0-9][0-9]", geno_db_m$ALLELES)
  nra <- F
  if (any(id_nra)) {
    # number ambiguous alleles
    num_text <- paste0('[*',1:length(unique(geno_db_m$ALLELES[id_nra])),']')
    names(num_text) <- unique(geno_db_m$ALLELES[id_nra])
    # text for plot
    geno_db_m$text[id_nra] <- num_text[geno_db_m$ALLELES[id_nra]]
    # text for legend
    geno_db_m$text_bottom[id_nra] <- paste(num_text[geno_db_m$ALLELES[id_nra]],geno_db_m$ALLELES[id_nra])
    # change allele to NRA - non reliable allele
    geno_db_m$ALLELES[id_nra] <- "NRA"
    # indicates that nra exists
    nra <- T
  }
  # create allele palette
  allele_palette <- allelePalette(geno_db_m$ALLELES)

  # sort novel allele calls for plot
  val_novel <- grep('^[0-9]+[_][A-Z][0-9]+[A-Z]',geno_db_m$ALLELES, value = T)
  novel <- F
  novel_allele_text <- c()
  novel_symbol <- "\u005E"
  if(length(val_novel)!=0){
    # sort the palettle colors for novel alleles
    id <- grep('^[0-9]+[_][A-Z][0-9]+[A-Z]',names(allele_palette$transper))
    allele_palette$transper[id] <- 1
    # cerate code index for novel allele
    code_allele <- paste0(novel_symbol,1:length(id))
    names(code_allele) <-allele_palette$AlleleCol[id]
    new_allele <- paste0(novel_symbol,1:length(id),'-',allele_palette$AlleleCol[id])
    names(new_allele) <-allele_palette$AlleleCol[id]
    # change the text for plot
    ids <- geno_db_m$ALLELES %fin% names(new_allele)
    rep <- new_allele[geno_db_m$ALLELES[ids]]
    rep2 <- code_allele[geno_db_m$ALLELES[ids]]
    # add new allele code to data
    geno_db_m[ids, c("ALLELES","text_bottom","text") := list(rep,rep,rep2)]
    # change annotation in legend colors
    allele_palette$AlleleCol[id] <- new_allele
    names(allele_palette$transper)[id] <- new_allele
    # indicates that novel exists
    novel <- T
  }

  geno_db_m$ALLELES <- factor(geno_db_m$ALLELES, levels = allele_palette$AlleleCol)


  # samples names and number
  samples <- unique(geno_table$SUBJECT)
  samples_n <- length(samples)
  # genes names and number
  genes <- unique(geno_db_m$GENE)
  genes_n <- length(genes)

  # order the data by gene loc
  setorderv(geno_db_m, c("SUBJECT","GENE_LOC"))
  setkey(geno_db_m, SUBJECT)

  # sort data for matrix
  geno_db_m[,line:=12/n]
  allele_code <- sapply(strsplit(gsub("\\^[0-9]+[-]","",allele_palette$AlleleCol),"_",fixed = T), "[[",1, USE.NAMES = T)
  last <- as.numeric(allele_code[length(allele_code)-3])+1
  ids_a <- last:(last+2)
  allele_code[(length(allele_code)-2):length(allele_code)] <- ids_a
  names(allele_code) <- gsub("\\^[0-9]+[-]","",allele_palette$AlleleCol)
  # sort the alleles in gene box
  geno_db_m[,A_CODE:=as.numeric(allele_code[ALLELES])]
  geno_db_m[grep("[0-9]_[0-9]",geno_db_m$ALLELES,perl = T),A_CODE:=allele_code["NRA"]]
  setorderv(geno_db_m, c("SUBJECT","GENE_LOC","A_CODE"))

  # duplicate the data by 12 box to gene
  geno_db_m[,id := 1:.N, by = .(SUBJECT, GENE)]
  geno_db_f = geno_db_m[,.(n_line = 1:line), by = .(SUBJECT, GENE, GENE_LOC, ALLELES_G, A_CODE,text_bottom), nomatch = 0]

  # transform allele codes to matrix, 12 box for each gene. each row is an individual
  m <- matrix(geno_db_f[[5]],ncol = 12*genes_n,byrow = T,dimnames = list(unique(geno_db_f[[1]]),geno_db_f[[2]][1:(12*genes_n)]))

  allele_code_t <- allele_palette$AlleleCol
  names(allele_code_t) <- allele_code

  geno_db_f[,text:=paste("Individual:",SUBJECT,"<br />Gene:",GENE,"<br />Allele:",text_bottom)]

  conditions.text <- matrix(geno_db_f[[8]], ncol = 12*genes_n, byrow = TRUE)
  #conditions.cols <- matrix(geno_db_f[[9]], ncol = 12*genes_n, byrow = TRUE)

  vline <- function(x = 0, color = "white") {
    list(
      type = "line",
      y0 = 0,
      y1 = 1,
      yref = "paper",
      x0 = x,
      x1 = x,
      line = list(color = color)
    )
  }

  kline <- function(NR,NC,X,Y, color = "white") {
    STEP_X<-1/(NC-1)
    STEP_Y<-1/(NR-1)
    list(list(
      type = "line",
      y0 = (Y-0.5)*STEP_Y,
      y1 = (Y)*STEP_Y,
      yref = "paper",
      x0 = (X+0.5)*STEP_X,
      x1 = (X+3.5)*STEP_X,
      xref = "paper",
      line = list(color = color)
    ),list(
      type = "line",
      y0 = (Y-0.5)*STEP_Y,
      y1 = (Y)*STEP_Y,
      yref = "paper",
      x0 = (X+4.5)*STEP_X,
      x1 = (X+7.5)*STEP_X,
      xref = "paper",
      line = list(color = color)
    ),list(
      type = "line",
      y0 = (Y-0.5)*STEP_Y,
      y1 = (Y)*STEP_Y,
      yref = "paper",
      x0 = (X+8.5)*STEP_X,
      x1 = (X+11.5)*STEP_X,
      xref = "paper",
      line = list(color = color)
    ))
  }

  col_names <- unique(sapply(strsplit(gsub("\\^[0-9]+[-]","",allele_palette$AlleleCol),"_",fixed = T), "[[",1, USE.NAMES = T))
  mypal <- colorRampPalette(unique(names(allele_palette$AlleleCol)))
  ncols = length(unique(names(allele_palette$AlleleCol)))#+1
  cols <- mypal(ncols)
  zseq <- seq(0,1,length.out=ncols+1)
  colorScale <- data.frame(
    z = c(0,rep(zseq[-c(1,length(zseq))],each=2),1),
    col=rep(cols,each=2)
  )
  colorScale$col <- as.character(colorScale$col)

  zmx <- round(max(m))
  zmn <- round(min(m))
  colorbar=list(tickmode='array', tick0=-zmn, dtick=1,tickvals= seq(0.5, nrow(m), length.out = (ncols+1)), ticktext=c("",col_names),
                len = 1, outlinecolor="white",bordercolor="white",borderwidth=5,bgcolor="white")


  # add grid lines
  gridlines <- lapply(seq(11.5,genes_n*12,by=12),vline)


  # plot dim
  plot_height <- 500 + 10*nrow(m)
  plot_width <- 100 + 6*ncol(m)
  # create plot
  p <- plotly::plot_ly(z=(m),type = "heatmap",
                       colorscale= colorScale,
                       colorbar = colorbar,
                       hoverinfo='text',text=conditions.text, width = plot_width, height = plot_height) %>%
    plotly::layout(yaxis = list(dtick = 1, ticktext = rownames(m), tickmode="array", tickvals = 0:(nrow(m)-1)),
                   xaxis = list(dtick = 1, ticktext = unique(colnames(m)), tickmode="array", tickvals = seq(6,12*genes_n,12)))

  # add k lines
  klines = geno_db_m[geno_db_m$K<lk_cutoff,]
  if(nrow(klines)>0){
    klines[, y:=match(SUBJECT,samples)] # row index
    klines[, yend:=y+0.5] # row index
    klines[, x:=(as.numeric(GENE_LOC)-1)*12] # col index
    klines[, xend:=x+1] # col index
    klines2 <- apply(klines, 1,function(x) kline(NR,NC,as.numeric(x["x"]),as.numeric(x["y"])))

    p <- p %>% plotly::layout(shapes = c(gridlines,unlist(klines2,recursive = F)))
  }

  # add text annotations
  ids_text <- grep('^[0-9]|Del|Unk',geno_db_m$text_bottom,invert = T)
  if(length(ids_text)>0){
    annot = geno_db_m[ids_text,]
    annot[, y:=(match(SUBJECT,samples)-1)]
    annot[, x:=((as.numeric(GENE_LOC)-1)*12+as.numeric(id)*(12/n)-1.5 )]
    p <- p %>%  plotly::add_annotations(x = annot$x,
                                        y = annot$y,
                                        text = annot$text,
                                        xref = 'x',
                                        yref = 'y', size = 0.025, showarrow = FALSE, font=list(color='black',size=0.025))
  }


  # save html file
  htmlwidgets::saveWidget(p, file=file, selfcontained = F)
}
