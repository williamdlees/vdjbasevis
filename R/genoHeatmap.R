########################################################################################################
#' Graphical output of alleles division by genotypes
#'
#' The \code{genoHeatmap} function generates a graphical output of the alleles per gene in multiple samples.
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
#' @param    color_y              named list of the colors for y axis labels. Defualt is NULL.
#' @param    file                 file path for rendering the plot to pdf. If non is supplied than the plot is retured as object. Defualt is NULL.
#' @return
#'
#' The following list is returned:
#' \itemize{
#'   \item \code{'p'}:        heat-map visualization of the genotype inference for multiple samples.
#'   \item \code{'width'}:    Optimal width value for rendering plot.
#'   \item \code{'height'}:   Optimal width value for rendering plot.
#' }
#'
#' When a file is supplied the graph is also rendered directly to pdf.
#'
#' @details
#'
#' A \code{data.frame} created by \code{inferGenotypeBaysian}.
#'
#'
#' @export
genoHeatmap <- function(geno_table, chain = c("IGH", "IGK", "IGL"), gene_sort = "position", removeIGH = TRUE, lk_cutoff = 1, mark_low_lk = TRUE, html = FALSE, n_line = 4, line_length=60, pseudo_genes = FALSE, ORF_genes = FALSE, file = NULL, color_y = NULL) {
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
  geno_db <- setDT(geno_db)[CJ("SUBJECT" = geno_db$SUBJECT, "GENE" = geno_db$GENE, unique=TRUE), on=c("SUBJECT", "GENE")]
  geno_db[is.na(geno_db$ALLELES) , c("ALLELES","K") := list("Unk", NA_integer_)]
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
  geno_db_m <- geno_db[, n:=  .N, by = c("SUBJECT", "GENE")][] # count number of alleles for group
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
  setkey(geno_db_m, "SUBJECT")

  # sort data for matrix
  geno_db_m[,"line":=12/geno_db_m$n]
  allele_code <- 1:length(allele_palette$AlleleCol)
  names(allele_code) <- gsub("\\^[0-9]+[-]","",allele_palette$AlleleCol)
  # sort the alleles in gene box
  geno_db_m[,"A_CODE":=allele_code[geno_db_m$ALLELES]+1]
  geno_db_m[grep("[0-9]_[0-9]",geno_db_m$ALLELES,perl = T), "A_CODE":=allele_code["NRA"]]
  setorderv(geno_db_m, c("SUBJECT","GENE_LOC","A_CODE"))

  # duplicate the data by 12 box to gene
  geno_db_m[,"id" := 1:.N, by = c("SUBJECT", "GENE")]
  geno_db_f = geno_db_m[,c("n_line" = 1:get("line")), by = c("SUBJECT", "GENE", "GENE_LOC", "ALLELES_G", "A_CODE", "text_bottom"), nomatch = 0]

  # transform allele codes to matrix, 12 box for each gene. each row is an individual
  m <- matrix(geno_db_f[[5]],ncol = 12*genes_n,byrow = T,dimnames = list(unique(geno_db_f[[1]]),geno_db_f[[2]][1:(12*genes_n)]))

  #### sort legend
  # fit width of legend column text by the longest allele
  longest_allele <- max(nchar(allele_palette$AlleleCol))*3+40
  # legend length
  leg_length <- next_divisor(longest_allele,genes_n*12)
  # create legend values for matrix
  leg <- function(x,allele_text_length) c(rep(x,4),rep(0,allele_text_length))
  seqs <- lapply(allele_code+1,leg, allele_text_length = leg_length-4)
  # add values for legend to complete the matrix, white boxes
  add = ceiling(length(unlist(seqs))/(genes_n*12))*(genes_n*12) - length(unlist(seqs))
  add = ifelse(add<0,0,add)
  # legend matrix
  m2 <- matrix(c(unlist(seqs),rep(0,add)), ncol = genes_n*12, byrow = T)

  # start and end values for plot parts
  start <- c(0.01)
  end <- c(0.15,0.13)
  size = 1
  short_reads_rows = 0
  ## add short read text annotation at the bottom
  if(nra){
  bottom_annot <- unique(grep("[0-9][0-9]_[0-9][0-9]", geno_db_m$text_bottom , value = T, perl = T))

  # Create text for annotating the short reads labels
  annot <- splitlines(bottom_annot, genes_n*12/(4*size))
  annot <- annot[!is.na(dplyr::na_if(annot,""))]
  short_reads_rows = length(annot)
  annot <- paste0(annot,collapse = '\n')
  # add the start and value for the third part
  end <- c(0.2, 0.22, 0.06)
  start <- c(0.03,0.01)
  }

  # set the height and width of plot
  height <- samples_n * 0.1 + 2 + nrow(m2)*0.2 + short_reads_rows*0.4 # number of samples, number of rows in legend, number of rows in bottom annotation
  width <- genes_n * 0.3 + 1.5 # numer of genes
  size_text = nrow(m)/(height*width) # text size for heatmap annoations
  size_text_leg = ncol(m2)/(width*longest_allele)+1 # text size for legend annotations

  if(!is.null(file)) pdf(file,onefile = F, width = width, height = height, family = "serif")

  # plot layout
  layout.matrix <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)

  graphics::layout(mat = layout.matrix,
         heights = c(3, 0.5, 0.3) # Heights of the three rows
  )
  # heatmap main plot
  par(mar=c(2,6,8,6))
  image(t(m),col = names(allele_palette$AlleleCol), breaks = 1:(length(allele_palette$AlleleCol)+1),axes=F)
  # add grid lines for genes
  grid(lwd=1,nx = genes_n,ny=0,col = "white",lty = 1)
  # add axis annotations
  axis(3,(0:(genes_n-1))/genes_n+6/(12*genes_n),names(gene_loc),las=3) # top
  axis(1,(0:(genes_n-1))/genes_n+6/(12*genes_n),names(gene_loc),las=3) # bottom

  # color y tick labels if supplied
  colors <- "black"
  if(!is.null(color_y)) colors <- color_y[rownames(m)]

  Map(axis, side=2, at=(0:(samples_n-1))/(samples_n-1), col.axis=colors, labels=rownames(m), lwd=0, las=1, cex.axis=0.8) #left
  axis(2,at=(0:(samples_n-1))/(samples_n-1),labels=FALSE)

  # draw lines for low lk values
  sub_geno = geno_db_m[geno_db_m$K<lk_cutoff,]
  NR = samples_n
  NC = genes_n*12
  apply(sub_geno, 1,function(x){

    I = which(x["SUBJECT"]==samples)-1    # row index
    J = (as.numeric(x["GENE_LOC"])-1)*12            # column index
    draw_segment(NR,NC,I,J,lwd=1,col="white")}
  )

  # ad text annotations
  ids_text <- grep('^[0-9]|Del|Unk',geno_db_m$text_bottom,invert = T)
  sub_geno = geno_db_m[ids_text,]
  NR = samples_n
  NC = genes_n*12
  apply(sub_geno, 1,function(x){

    I = which(x["SUBJECT"]==samples)-1    # row index
    J = (as.numeric(x["GENE_LOC"])-1)*12            # column index
    ALLELE =  as.numeric(x["id"])                   # allele index
    N_ALLELES = as.numeric(x["n"])                  # number of alleles
    TEXT =  x["text"]                   # text
    Write_text(NR,NC,I,J,ALLELE,N_ALLELES,TEXT,cex=size_text)}
  )

  # legend plot
  par(mar=c(0,6,4,6))
  image(t(m2),col = names(allele_palette$AlleleCol), breaks = 1:(length(allele_palette$AlleleCol)+1),axes=F)
  # add grid lines for each row
  grid(lwd=1,ny = nrow(m2), nx = 0,col = "black",lty = 2)
  # add text anotation for legend
  NR = nrow(m2)
  NC = genes_n*12
  names(allele_code) <- allele_palette$AlleleCol
  invisible(tapply(allele_code,names(allele_code),function(x){
    ii = which(m2==x+1,arr.ind = T)[1,]
    I = ii[[1]]-1              # row index
    J = ii[[2]]                # column index
    ALLELE =  1                # allele index
    N_ALLELES = 1              # number of alleles
    TEXT =  names(x)            # text
    STEP_X<-1/(NC-1)
    STEP_Y<-ifelse(1/(NR-1)==Inf,0,1/(NR-1))
    text(STEP_X*J+STEP_X*leg_length*0.5,
         STEP_Y*I,
         TEXT,cex = size_text_leg)
  }
  ))
  # bottom text for nra, only if exists
  if(nra){
    # add the text at the bottom
    par(mar=c(1,6,1,6))
    gplots::textplot(annot, halign = "center", cex = size)
  }
  p1.base <- recordPlot()
  invisible(dev.off())
  return(list(p = p1.base, width = width, height = height))
  if(!is.null(file)) dev.off()
  # embed the fonts to file
  if(!is.null(file)) embedFonts(file)
}
