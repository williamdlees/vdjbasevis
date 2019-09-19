#' @include vdjbaseVis.R
NULL

########################################################################################################
#' Sort data frame by genes
#'
#' \code{sortDFByGene} Sort the \code{data.frame} by the genes names or position. For sorting by gene names the \code{sortAlleles} function by TIgGER is used.
#' For sorting by position the defualt package gene location list is used.
#'
#' @param    DATA                 data frame to sort
#' @param    chain                the Ig chain: IGH,IGK,IGL. Default is IGH.
#' @param    method               the method for sorting the genes. If by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#'
#' @return   sorted \code{data.frame}
#' @export
sortDFByGene <- function(DATA, chain = c("IGH", "IGK", "IGL"), method = c("name", "position"), removeIGH = FALSE, geno = FALSE, peseudo_remove = F, ORF_remove = F) {
  if (missing(chain)) {
    chain = "IGH"
  }
  chain <- match.arg(chain)

  if (missing(method)) {
    method = "position"
  }
  method <- match.arg(method)

  GENE.loc.tmp <- GENE.loc[[chain]]
  vs <- grep("V",GENE.loc.tmp,value = T)
  ps <- grep("V",PSEUDO[[chain]],value = T)
  ps <- ps[!ps %fin% vs]
  orf <- unique(grep("OR|NL", DATA$GENE,value = T,perl = T))
  GENE.loc.tmp <- c(vs,ps,orf,grep("V",GENE.loc.tmp,value = T,invert = T,perl = T))


  if(!peseudo_remove){
    DATA <- DATA[!DATA$GENE %fin% PSEUDO[[chain]],]
    GENE.loc.tmp <- GENE.loc.tmp[!GENE.loc.tmp %fin% ps]
  }
  if(!ORF_remove){
    DATA <- DATA[!DATA$GENE %fin% orf,]
    GENE.loc.tmp <- GENE.loc.tmp[!GENE.loc.tmp %fin% orf]
  }

  names(GENE.loc.tmp) <- GENE.loc.tmp

  if (method == "name") {
    DATA$GENE <- factor(DATA$GENE, levels = rev(sortAlleles(unique(DATA$GENE), method = method)))
    if (removeIGH) {
      DATA$GENE <- gsub("IG[H|K|L]", "", DATA$GENE)
      DATA$GENE <- factor(DATA$GENE, levels = rev(sortAlleles(unique(DATA$GENE), method = method)))
      if(!geno) DATA$hapBy <- gsub("IG[H|K|L]", "", DATA$hapBy)
    }
  } else {
    DATA$GENE <- factor(DATA$GENE, levels = rev(GENE.loc.tmp))
    if (removeIGH) {
      GENE.loc.tmp <- gsub("IG[H|K|L]", "", GENE.loc.tmp)
      names(GENE.loc.tmp) <- GENE.loc.tmp
      DATA$GENE <- gsub("IG[H|K|L]", "", DATA$GENE, perl = T)
      DATA$GENE <- factor(DATA$GENE, levels = rev(GENE.loc.tmp))
      if(!geno) DATA$hapBy <- gsub("IG[H|K|L]", "", DATA$hapBy)
    }
  }

  return(DATA)
}

########################################################################################################
#' Creates the allele color palette for haplotype graphical output
#'
#' \code{allelePalette} Takes a list of the haplotype alleles and returns the allele color palette.
#'
#' @param   hap_alleles          a list of the haplotype alleles.
#'
#' @return   Haplotype allele color palette
#' @export
allelePalette <- function(hap_alleles, NRA = TRUE) {

  Alleles <- grep("[012]", unique(hap_alleles), value = T, perl = T)
  AlleleCol.tmp <- sort(unique(sapply(strsplit(Alleles, "_"), "[", 1)))
  tmp.col <- ALLELE_PALETTE[AlleleCol.tmp]


  novels <- grep("_", Alleles, value = T)
  if (length(novels) > 0) {
    novels.col <- ALLELE_PALETTE[sapply(strsplit(novels, "_"), "[", 1)]
    names(novels.col) <- novels
    alleles.comb <- c(tmp.col, novels.col)[order(names(c(tmp.col, novels.col)))]
  } else {
    alleles.comb <- c(tmp.col)[order(names(c(tmp.col)))]

  }

  AlleleCol <- names(c(alleles.comb, Unk = "#dedede", Del = "#6d6d6d", NR = "#000000", NRA = "#fbf7f5"))
  names(AlleleCol) <- c(alleles.comb, Unk = "#dedede", Del = "#6d6d6d", NR = "#000000", NRA = "#fbf7f5")
  rm_allele <- function(allele,alleles,AlleleCol){
    if(!allele %in% alleles){
      id <- which(allele == AlleleCol)
      return(AlleleCol[-id])
    }
    return(AlleleCol)
  }
  AlleleCol <- rm_allele("NR",hap_alleles,AlleleCol)
  AlleleCol <- rm_allele("Del",hap_alleles,AlleleCol)
  AlleleCol <- rm_allele("Unk",hap_alleles,AlleleCol)
  AlleleCol <- rm_allele("NRA",hap_alleles,AlleleCol)


  transper <- sapply(AlleleCol, function(x) {
    if (grepl("_", x)) {
      mom_allele <- strsplit(x, "_")[[1]][1]
      all_novel <- grep(paste0(mom_allele, "_"), AlleleCol, value = T)
      if (length(all_novel) == 1) {
        return(0.5)
      }
      if (length(all_novel) == 2) {
        m = which(all_novel == x)
        return(ifelse(m == 1, 0.6, 0.3))
      }
      if (length(all_novel) == 3) {
        m = which(all_novel == x)
        if (m == 1) {
          return(0.6)
        }
        return(ifelse(m == 2, 0.4, 0.2))
      }
      if (length(all_novel) > 9) {
        m = which(all_novel == x)
        if (m == 1) {
          return(1)
        }
        return(1 - m/20)
      }
      if (length(all_novel) > 3) {
        m = which(all_novel == x)
        if (m == 1) {
          return(0.85)
        }
        return(0.85 - m/10)
      }
    } else (1)
  })
  names(transper) <- AlleleCol

  # remove 'mother' allele if added (when there is no germline allele but there is a novel)

  special <- c("Unk", "Del", "NR", "NRA")[c("Unk", "Del", "NR", "NRA") %in% AlleleCol]

  AlleleCol <- AlleleCol[AlleleCol %in% c(sort(grep("[012]", unique(hap_alleles), value = T, perl = T)), special)]

  transper <- transper[names(transper) %in% AlleleCol]

  return(list(transper = transper, AlleleCol = AlleleCol))
}

########################################################################################################
# Get diagonal line for legend
#
getDigLegend <- function(color){

  return(ggplotGrob(ggplot(data.frame(x=c(1,2),y=c(3,4)), aes_string("x","y")) + geom_abline(aes_string(colour="color", intercept = 1, slope = 1), show.legend = T) +
                      scale_color_manual(values = c("white"), name = "", drop = FALSE) + guides(color = guide_legend(override.aes = list(size = 0.5), order = 2)) +
                      theme(legend.justification = "center", legend.key = element_rect(fill = "gray"), legend.position = "bottom")))
}

########################################################################################################
# Get a list of patterns, replacments, and a string
#
# Returns the replaced string
#
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern) != length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}

########################################################################################################
#' Creates the non reliable allele text annotation for plots
#'
#' \code{nonReliableAllelesText_V2} Takes the haplotype data frame
#'
#' @param   geno_table          a data frame of the haplotypes.
#'
#' @return   Non reliable alleles text data frame for plots annotation.
#' @export
nonReliableAllelesText_V2 <- function(non_reliable_alleles_text, size = 3, map = F) {

  id <- grepl("[0-9][0-9]_[0-9][0-9]", non_reliable_alleles_text$ALLELES)
  num_text <- sapply(1:length(unique(non_reliable_alleles_text$ALLELES[id])),function(i) paste0('[*',i,']'))
  names(num_text) <- unique(non_reliable_alleles_text$ALLELES[id])
  non_reliable_alleles_text$text <- ''
  non_reliable_alleles_text$text[id] <- num_text[non_reliable_alleles_text$ALLELES]
  non_reliable_alleles_text$text_bottom <- ''
  non_reliable_alleles_text$text_bottom[id] <- paste(num_text[non_reliable_alleles_text$ALLELES[id]],non_reliable_alleles_text$ALLELES[id])
  non_reliable_alleles_text$pos <- ifelse(non_reliable_alleles_text$freq == 1, 1,
                                          ifelse(non_reliable_alleles_text$freq == 2, 0.5,
                                                 ifelse(non_reliable_alleles_text$freq == 3, 0.33, 0.25)))
  non_reliable_alleles_text$size <- size
  if(!map){
    non_reliable_alleles_text <- non_reliable_alleles_text %>% ungroup() %>% group_by(GENE, SUBJECT, hapBy) %>%
      dplyr::mutate(pos = pos + ifelse(dplyr::row_number()==2,dplyr::row_number()-1.5,dplyr::row_number()-1))
    }else{
    non_reliable_alleles_text <- non_reliable_alleles_text %>% ungroup() %>% group_by(GENE, SUBJECT) %>%
      dplyr::mutate(pos = ifelse(n == 1, 0.5,
                          ifelse(n == 2, seq(0.25,1,by = 0.5)[1:max(dplyr::row_number())],
                                 ifelse(n == 3, seq(0.165,1,by = 0.33)[1:max(dplyr::row_number())],
                                        seq(0.125,1,by = 0.25)[1:max(dplyr::row_number())]))))
  }
  #non_reliable_alleles_text$NEW_ALLELES <- non_reliable_alleles_text$ALLELES
  non_reliable_alleles_text$ALLELES[id] <- "NRA"
  return(non_reliable_alleles_text)
}

########################################################################################################
#' Creates the novel allele text annotation for plots
#'
#' \code{novelAlleleAnnotation} Takes the haplotype data frame
#'
#' @param   novel_allele         data frame with the novel allele cordinates.
#'
#' @return   novel alleles text data frame for plots annotation.
#' @export
novelAlleleAnnotation <- function(novel_allele, new_label, size = 3) {
  if (nrow(novel_allele) != 0) {
    novel_allele$text <- sapply(new_label[novel_allele$ALLELES],function(s) strsplit(s,'-')[[1]][1])
    novel_allele$text_bottom <- paste(new_label[novel_allele$ALLELES],novel_allele$ALLELES)
    novel_allele$pos <- ifelse(novel_allele$freq == 1, 1,
                               ifelse(novel_allele$freq == 2, 0.5,
                                      ifelse(novel_allele$freq == 3, 0.33, 0.25)))
    novel_allele$size = size
    novel_allele <- novel_allele %>% dplyr::ungroup() %>% dplyr::group_by(GENE, SUBJECT) %>%
      dplyr::mutate(pos = ifelse(n == 1, 0.5,
                          ifelse(n == 2, seq(0.25,1,by = 0.5)[1:max(dplyr::row_number())],
                                 ifelse(n == 3, seq(0.165,1,by = 0.33)[1:max(dplyr::row_number())],
                                        seq(0.125,1,by = 0.25)[1:max(dplyr::row_number())]))))
    return(novel_allele)
  } else {
    return(setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("GENE", "ALLELES", "n", "freq", "text", "pos", "size")))
  }
}

########################################################################################################
#' Split lines of short reads assignments
#'
#' The \code{splitlines} function sliplits the text by line width.
#'
#' @param    bottom_annot         annotation text.
#' @param    line_width           the line width allowed.
#'
#' @return Seperated text to lines
#'
#' @export
splitlines<-function(bottom_annot,line_width=60){
  if(line_width<=max(sapply(bottom_annot,nchar))){
    line_width = max(sapply(bottom_annot,nchar))+2
    print(paste0("Set line width to ",line_width))
  }
  collapsed_annot<-paste(bottom_annot,collapse=", ")
  L<-nchar(collapsed_annot)
  if(L<line_width)return(collapsed_annot)
  vec_annot<-substring(collapsed_annot,1:L,1:L)
  ind<-grep("[",vec_annot,fixed=TRUE)
  aligned_text<-NULL
  #i<-line_width
  i_previous<-1
  while(i_previous<=(L-line_width)){
    temp<-findInterval(line_width*(length(aligned_text)+1)+1,ind)
    aligned_text<-c(aligned_text,substring(collapsed_annot,i_previous,ind[temp]-1))
    i_previous<-ind[temp]
  }
  aligned_text<-c(aligned_text,substring(collapsed_annot,i_previous,L))
  return(aligned_text)
}

########################################################################################################
#' Creates the novel allele text annotation for plots
#'
#' \code{novelAlleleAnnotation} Takes the haplotype data frame
#'
#' @param   novel_allele         data frame with the novel allele cordinates.
#'
#' @return   novel alleles text data frame for plots annotation.
#' @export
novelAlleleAnnotation_geno <- function(novel_allele, new_label, size = 3) {
  if (nrow(novel_allele) != 0) {
    novel_allele$text <- sapply(new_label[novel_allele$ALLELES],function(s) strsplit(s,'-')[[1]][1])
    novel_allele$text_bottom <- paste(new_label[novel_allele$ALLELES],novel_allele$ALLELES)
    novel_allele$pos <- ifelse(novel_allele$freq == 1, 1,
                               ifelse(novel_allele$freq == 2, 0.5,
                                      ifelse(novel_allele$freq == 3, 0.33, 0.25)))
    novel_allele$size = size
    novel_allele <- novel_allele %>% dplyr::ungroup() %>% dplyr::group_by(GENE, SUBJECT) %>%
      dplyr::mutate(pos = ifelse(n == 1, 0.5,
                          ifelse(n == 2, seq(0.25,1,by = 0.5)[1:max(dplyr::row_number())],
                                 ifelse(n == 3, seq(0.165,1,by = 0.33)[1:max(dplyr::row_number())],
                                        seq(0.125,1,by = 0.25)[1:max(dplyr::row_number())]))))
    return(novel_allele)
  } else {
    return(setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("GENE", "ALLELES", "n", "freq", "text", "pos", "size")))
  }
}

########################################################################################################
#' Write text annotations for heatmap graph
#'
#' \code{Write_text} takes values for plotting text on heatmap
#'
#' @param   NR           number of rows.
#' @param   NC           number of columns.
#' @param   I            row index.
#' @param   J            column index.
#' @param   ALLELE       allele index in the individual gene box.
#' @param   N_ALLELES    number of alleles for individual in gene box.
#' @param   TEXT         annotation text.
#'
#' @return   plotting text annotation.
#' @export
Write_text<-function(NR,NC,I,J,ALLELE,N_ALLELES,TEXT,...){
  STEP_X<-1/(NC-1)
  STEP_Y<-1/(NR-1)
  text(STEP_X*J-STEP_X/2+STEP_X*12/N_ALLELES*(ALLELE-1/2),
       STEP_Y*I,
       TEXT,...)
}
########################################################################################################
#' Draw lk value lines on heatmap
#'
#' \code{draw_segment} takes values for plotting text on heatmap
#'
#' @param   NR           number of rows.
#' @param   NC           number of columns.
#' @param   I            row index.
#' @param   J            column index.
#'
#' @return   plotting lk lines.
#' @export
draw_segment<-function(NR,NC,I,J,...){
  STEP_X<-1/(NC-1)
  STEP_Y<-1/(NR-1)
  points(c(STEP_X*(J-0.5),STEP_X*(J+1.5)),
         c(STEP_Y*(I),STEP_Y*(I+0.5)), type = "l",...)
  points(c(STEP_X*(J-0.5),STEP_X*(J+3.5)),
         c(STEP_Y*(I-0.5),STEP_Y*(I+0.5)), type = "l",...)
  points(c(STEP_X*(J+1.5),STEP_X*(J+5.5)),
         c(STEP_Y*(I-0.5),STEP_Y*(I+0.5)), type = "l",...)
  points(c(STEP_X*(J+3.5),STEP_X*(J+7.5)),
         c(STEP_Y*(I-0.5),STEP_Y*(I+0.5)), type = "l",...)
  points(c(STEP_X*(J+5.5),STEP_X*(J+9.5)),
         c(STEP_Y*(I-0.5),STEP_Y*(I+0.5)), type = "l",...)
  points(c(STEP_X*(J+7.5),STEP_X*(J+11.5)),
         c(STEP_Y*(I-0.5),STEP_Y*(I+0.5)), type = "l",...)
  points(c(STEP_X*(J+9.5),STEP_X*(J+11.5)),
         c(STEP_Y*(I-0.5),STEP_Y*(I)), type = "l",...)
}
########################################################################################################
# finds next dividor for an int number
next_divisor <- function(x,y){
  if(x>y) return(y)
  while(T){
    if(y%%x==0) return(x)
    x <- x+1
  }
}

########################################################################################################
#' Sort multiple genotypes
#'
#' The \code{sortGenotype} function sort multiple genotype tables.
#'
#'
#' @param    geno_table           genoytpe summary table. See details.
#' @param    chain                the IG chain: IGH,IGK,IGL. Default is IGH.
#' @param    gene_sort            if by 'name' the genes in the output are ordered lexicographically,
#' if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
#' @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#' @param    lk_cutoff            the lK cutoff value to be considerd low for texture layer. Defualt is lK<1.
#'
#' @return
#'
#' A \code{data.frame} of the genotypes.
#'
#' @details
#'
#' A \code{data.frame} created by \code{inferGenotypeBaysian}.
#'
#'
#' @export
sortGenotype <- function(geno_table, chain = c("IGH", "IGK", "IGL"), gene_sort = "position", removeIGH = TRUE, lk_cutoff = 1){
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
  return(geno_db_m)
}
