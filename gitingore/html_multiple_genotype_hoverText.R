# Or - Create Multiple Genotype Graphs html/pdf("T"/"F").

#### need packeges: ####
#install.packages(c("mltools","optparse"))
#for function "bin_data"
library('mltools')
#for "do.call" function
library('base')
# for "get_legend" function
library('cowplot')
library('ggplot2')
#library('caroline')
library(optparse)
library(magrittr) # need to run every time you start R and want to use %>%
library(stats)
library(dplyr)    # alternative, this also loads %>%
#colors
library(RColorBrewer)
library(wesanderson)
library(scales)
library('tigger')
library('plotly')
# dcd
library(stringr)
#drop NA
library(tidyr)

########## VDJbase server ##############

option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="excel file name"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="graph.pdf file name"),
  make_option(c("-s", "--sysdata_file"), type="character", default=NULL,
              help="sysdata file name"),
  make_option(c("-t", "--is_html"), type="character", default=NULL,
              help="type of file F - pdf T - html"),
  make_option(c("-p", "--with_pseudo"), type="character", default="F",
              help="With/out pseudo and ORF genes: F - without T - with"),
  make_option("--samp", type="character", default=NULL,
              help="Sample name")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_file)){
  stop("input file must be supplied", call.=FALSE)
}

if (is.null(opt$output_file)){
  stop("output reference file must be supplied", call.=FALSE)
}

if (is.null(opt$sysdata_file)){
  stop("sysdata file must be supplied", call.=FALSE)
}

if (is.null(opt$is_html)){
  stop("type of file must be supplied", call.=FALSE)
}


######### VDJbase server- loading data(use melt function) #############

# read genotype table
genotype_path<-opt$input_file
output_file<-opt$output_file
data<- read.delim(file= genotype_path ,header=TRUE,sep="\t",stringsAsFactors = F)

if (!is.null(opt$samp)) {
  data$SUBJECT <- opt$samp
}
# load the "sysdata"
#load("/home/aviv/Rscripts/sysdata.rda")
load(opt$sysdata_file)

html_output <- as.logical(opt$is_html)  # for pdf set "F"
pseudo_ORF_genes <- as.logical(opt$with_pseudo)
# html_output <- ifelse(html_output == "T", TRUE, FALSE)


########################################################################################################
# Creates the allele color palette for haplotype graphical output
#
# \code{allelePalette} Takes a list of the haplotype alleles and returns the allele color palette.
#
# @param   hap_alleles          a list of the haplotype alleles.
#
# @return   Haplotype allele color palette
########################################################################################################
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
# Creates the non reliable allele text annotation for plots
#
# \code{nonReliableAllelesText_V2} Takes the haplotype data frame
#
# @param   hap_table          a data frame of the haplotypes.
#
# @return   Non reliable alleles text data frame for plots annotation.
#
nonReliableAllelesText_V2 <- function(non_reliable_alleles_text, size = 3, map = F) {

  if (nrow(non_reliable_alleles_text) != 0) {
    num_text <- sapply(1:length(unique(non_reliable_alleles_text$ALLELES)),function(i) paste0('[*',i,']'))
    names(num_text) <- unique(non_reliable_alleles_text$ALLELES)
    non_reliable_alleles_text$text <- num_text[non_reliable_alleles_text$ALLELES]
    non_reliable_alleles_text$text_bottom <- paste(num_text[non_reliable_alleles_text$ALLELES],non_reliable_alleles_text$ALLELES)
    non_reliable_alleles_text$pos <- ifelse(non_reliable_alleles_text$freq == 1, 1,
                                            ifelse(non_reliable_alleles_text$freq == 2, 0.5,
                                                   ifelse(non_reliable_alleles_text$freq == 3, 0.33, 0.25)))
    non_reliable_alleles_text$size = size
    if(!map){
      non_reliable_alleles_text <- non_reliable_alleles_text %>% ungroup() %>% group_by(.data$GENE, .data$SUBJECT, .data$hapBy) %>%
        mutate(pos = .data$pos + ifelse(dplyr::row_number()==2,dplyr::row_number()-1.5,dplyr::row_number()-1))}
    else{
      non_reliable_alleles_text <- non_reliable_alleles_text %>% ungroup() %>% group_by(.data$GENE, .data$SUBJECT) %>%
        mutate(pos = ifelse(.data$n == 1, 0.5,
                            ifelse(.data$n == 2, seq(0.25,1,by = 0.5)[1:max(dplyr::row_number())],
                                   ifelse(.data$n == 3, seq(0.165,1,by = 0.33)[1:max(dplyr::row_number())],
                                          seq(0.125,1,by = 0.25)[1:max(dplyr::row_number())]))))
    }

    non_reliable_alleles_text$ALLELES[grep("[0-9][0-9]_[0-9][0-9]", non_reliable_alleles_text$ALLELES)] <- "NRA"
    return(non_reliable_alleles_text)
  } else {
    if(!map) return(setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("GENE", "ALLELES", "hapBy", "n", "freq", "text", "pos", "size")))
    else return(setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("GENE", "ALLELES", "n", "freq", "text", "pos", "size")))
  }
}

########################################################################################################
# Creates the novel allele text annotation for plots
#
# \code{novelAlleleAnnotation} Takes the haplotype data frame
#
# @param   novel_allele         data frame with the novel allele cordinates.
#
# @return   novel alleles text data frame for plots annotation.
#
novelAlleleAnnotation_geno <- function(novel_allele, new_label, size = 3) {
  if (nrow(novel_allele) != 0) {
    novel_allele$text <- sapply(new_label[novel_allele$ALLELES],function(s) strsplit(s,'-')[[1]][1])
    novel_allele$text_bottom <- paste(new_label[novel_allele$ALLELES],novel_allele$ALLELES)
    novel_allele$pos <- ifelse(novel_allele$freq == 1, 1,
                               ifelse(novel_allele$freq == 2, 0.5,
                                      ifelse(novel_allele$freq == 3, 0.33, 0.25)))
    novel_allele$size = size
    novel_allele <- novel_allele %>% ungroup() %>% group_by(.data$GENE, .data$SUBJECT) %>%
      mutate(pos = ifelse(.data$n == 1, 0.5,
                          ifelse(.data$n == 2, seq(0.25,1,by = 0.5)[1:max(dplyr::row_number())],
                                 ifelse(.data$n == 3, seq(0.165,1,by = 0.33)[1:max(dplyr::row_number())],
                                        seq(0.125,1,by = 0.25)[1:max(dplyr::row_number())]))))
    return(novel_allele)
  } else {
    return(setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("GENE", "ALLELES", "n", "freq", "text", "pos", "size")))
  }
}
########################################################################################################
# Split lines of short reads assignments
#
# The \code{splitlines} function sliplits the text by line width
#
# @param    bottom_annot         annotation text.
# @param    line_width           the line width allowed.
#
# @return
# Seperated text to lines
#
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
# Sort data frame by genes
#
# \code{sortDFByGene} Sort the \code{data.frame} by the genes names or position. For sorting by gene names the \code{sortAlleles} function by TIgGER is used.
# For sorting by position the defualt package gene location list is used.
#
# @param    DATA                 data frame to sort
# @param    chain                the Ig chain: IGH,IGK,IGL. Default is IGH.
# @param    method               the method for sorting the genes. If by 'name' the genes in the output are ordered lexicographically,
# if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.
# @param    removeIGH            if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.
#
# @return   sorted \code{data.frame}
#
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
  ps <- ps[!ps %in% vs]
  orf <- unique(grep("OR|NL", DATA$GENE,value = T))
  GENE.loc.tmp <- c(vs,ps,orf,grep("V",GENE.loc.tmp,value = T,invert = T))


  if(!peseudo_remove){
    DATA <- DATA[!(DATA$GENE %in% PSEUDO[[chain]]),]
    GENE.loc.tmp <- GENE.loc.tmp[!GENE.loc.tmp %in% ps]
  }
  if(!ORF_remove){
    DATA <- DATA[!grepl("OR|NL", DATA$GENE),]
    GENE.loc.tmp <- GENE.loc.tmp[!GENE.loc.tmp %in% orf]
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
      DATA$GENE <- gsub("IG[H|K|L]", "", DATA$GENE)
      DATA$GENE <- factor(DATA$GENE, levels = rev(GENE.loc.tmp))
      if(!geno) DATA$hapBy <- gsub("IG[H|K|L]", "", DATA$hapBy)
    }
  }

  return(DATA)
}

########################################################################################################
# Creates genotype graphical output division by subjects in pdf/html format
#
# \code{multipleGenoytpe} Takes a genoytpe table and returns genotype graphical output.
#
# @param   gen_table          a table of the genoytpes.
# @param   html               character "F" or "T". "T" for html and "F" for pdf output.
#
# @return   genotype graphical output
########################################################################################################
multipleGenoytpe <- function(gen_table, chain = "IGH", html = FALSE, removeIGH = TRUE, text_size = 14, gene_sort='position', facet_by='SUBJECT', pseudo_genes = F, ORF_genes=F){


      #~~~~~~~~~~~~~~~~~~~~~~~~~~ arrange genotype work file ~~~~~~~~~~~~~~~~~~~~~~~
      # change single number (1,2,3) to allele notation (01,02,03)
      i1 <- grepl("^[0-9]$", gen_table$GENOTYPED_ALLELES)
      gen_table$GENOTYPED_ALLELES[i1] <- paste0("0", gen_table$GENOTYPED_ALLELES[i1])
      # Replace Deletion annotation
      gen_table$GENOTYPED_ALLELES <- gsub("Deletion","Del",gen_table$GENOTYPED_ALLELES)

      # select columns
      genotype_all <- gen_table%>% select(SUBJECT,GENE,ALLELES,COUNTS,TOTAL,K_DIFF,GENOTYPED_ALLELES,Freq_by_Clone)

      # change columns names Genotyped_alleles to Alleles
      genotype_all <- genotype_all %>% rename(temp=GENOTYPED_ALLELES)
      genotype_all<- genotype_all %>% rename(GENOTYPED_ALLELES=ALLELES)
      genotype_all <- genotype_all %>% rename(ALLELES=temp)
      genotype<-genotype_all

      # # for filter subjects or genes!!!!!!!!!!
      genotype <- subset(genotype, SUBJECT %in% unique(genotype$SUBJECT)[1:15])
      # genotype <- subset(genotype,  GENE%in% unique(genotype_all$GENE)[] )
      # # end: for filter subjects or genes!!!!!!!!!!!!!!!

      color_pes_orf <- c()
      if(pseudo_genes){
        color_pes_orf <- c(grep("V",PSEUDO[[chain]],value = T),color_pes_orf)
      }
      if(ORF_genes){
        color_pes_orf <- c(unique(grep("OR|NL", genotype$GENE,value = T)),color_pes_orf)
      }


      genotype <- sortDFByGene(DATA = genotype, chain = chain, method = gene_sort, removeIGH = removeIGH, geno = T,
                              peseudo_remove = pseudo_genes, ORF_remove = ORF_genes)

      # # remove pseudo genes + sort genes by loci for genotype graph and k_diff graph
      # GENE.loc.tmp <- GENE.loc[[chain]]  # take just IGH genes from the gene loci
      # names(GENE.loc.tmp) <- GENE.loc.tmp
      # # remove IGH from genes
      # if(removeIGH){
      #   GENE.loc.tmp <- gsub('IG[H|K|L]','',GENE.loc.tmp)
      #   names(GENE.loc.tmp) <- GENE.loc.tmp
      #   genotype$GENE <- gsub('IG[H|K|L]','',genotype$GENE)
      #   genotype$GENE = factor(genotype$GENE, levels = rev(GENE.loc.tmp))
      # } else {
      #   names(GENE.loc.tmp) <- GENE.loc.tmp
      #   genotype$GENE = factor(genotype$GENE, levels = rev(GENE.loc.tmp))
      # }
      #

      # split all multi alleles to number of raws
      alleles = strsplit(as.character(genotype$ALLELES), ",")
      count_allele = strsplit(as.character(genotype$Freq_by_Clone), ",")
      geno2 = genotype
      r = 1
      for (g in 1:nrow(genotype)) {
        for (a in 1:length(alleles[[g]])) {
          geno2[r, ] = genotype[g, ]
          geno2[r, ]$ALLELES = alleles[[g]][a]
          geno2[r, ]$Freq_by_Clone = count_allele[[g]][a]
          r = r + 1
        }
      }

      # Omit rows containing specific ALLELES/K_DIFF/GENE column of NA
      geno2  <-geno2 %>% drop_na(ALLELES)
      geno2  <-geno2 %>% drop_na(K_DIFF)
      geno2  <-geno2 %>% drop_na(GENE)
      #geno2  <-geno2 %>% drop_na(Freq_by_Clone)
      genotype <- genotype %>% drop_na(GENE)

      #~~~~~~~~~~~~~~~~~~~~~~~~~~ save and group K values ~~~~~~~~~~~~~~~~~~~~~~~
      kval.df = subset(genotype,select=c(SUBJECT,GENE,GENOTYPED_ALLELES,K_DIFF,Freq_by_Clone)) # Kval per gene

      # genes <- unique(kval.df$GENE)
      # for(s in unique(kval.df$SUBJECT)){
      #   sub <- kval.df[kval.df$SUBJECT==s,]
      #   # get missing genes
      #   id <- which(!genes %in% sub$GENE)
      #   if(length(id) != 0) kval.df <- bind_rows(kval.df,data.frame(SUBJECT = s, GENE = genes[id], K_DIFF = 1000,
      #                                                                 GENOTYPED_ALLELES = "Unk",
      #                                                                 Freq_by_Clone = NA, stringsAsFactors = F))
      # }

      #kval.df = subset(geno2,select=c(SUBJECT,GENE,GENOTYPED_ALLELES,K_DIFF,Freq_by_Clone)) # Kval per gene
      # categorize to groups
      kval.df$K_GROUPED <- bin_data(kval.df$K, bins=c(0, 1,2,3,4,5,10,20,50,Inf), binType = "explicit")
      kval.df <- kval.df[!is.na(na_if(kval.df$GENOTYPED_ALLELES,"")),]
      #~~~~~~~~~~~~~~~~~~~~~~~~~ select columns ~~~~~~~~~~~~~~~~~~~~~~~
      geno2 <- geno2 %>% select(SUBJECT,GENE,ALLELES,K_DIFF,GENOTYPED_ALLELES,Freq_by_Clone)
      kval.df2 <- kval.df %>% rename(ALLELES=K_GROUPED)
      unique_k <- unique(kval.df$K_GROUPED)

      # # define ratio between Kdiff to allele appearance, 3/4 [KK 0102 0102 0102]
      #geno2_kval <- rbind(geno2,geno2,geno2,kval.df2)
      genes <- unique(geno2$GENE)
      for(s in unique(geno2$SUBJECT)){
        sub <- geno2[geno2$SUBJECT==s,]
        # get missing genes
        id <- which(!genes %in% sub$GENE)
        if(length(id) != 0) geno2 <- bind_rows(geno2,data.frame(SUBJECT = s, GENE = genes[id], ALLELES = 'Unk', K_DIFF = 1000,
                                                                GENOTYPED_ALLELES = "Unk",
                                                                Freq_by_Clone = NA, stringsAsFactors = F))
      }



      geno2 <- as.data.frame(geno2 %>% group_by(.data$SUBJECT, .data$GENE) %>% mutate(n = dplyr::n()))
      #geno2$freq <- ifelse(geno2$n == 1, 0.75, ifelse(geno2$n == 2, 0.375, ifelse(geno2$n == 3, 0.25, 0.1875)))
      geno2 <- geno2 %>% group_by(.data$GENE, .data$SUBJECT) %>%
        mutate(freq = ifelse(.data$n == 1, 0.75,
                              ifelse(.data$n == 2, rep(0.375,2),
                                     ifelse(.data$n == 3, rep(0.25,3),
                                            rep(0.1875,4)))))
      kval.df2$freq <- 0.25
      geno2$ALLELE_TEXT <- geno2$ALLELES
      dagger = "†"
      if (length(grep("[0-9][0-9]_[0-9][0-9]", geno2$ALLELES)) != 0) {
        non_reliable_alleles_text <- nonReliableAllelesText_V2(non_reliable_alleles_text = geno2[grep("[0-9][0-9]_[0-9][0-9]", geno2$ALLELES), ],
                                                               map = T, size = 2)
        geno2$ALLELES[grep("[0-9][0-9]_[0-9][0-9]", geno2$ALLELES)] <- "NRA"
        allele_palette <- allelePalette(geno2$ALLELES)
        non_reliable_alleles_text$ALLELES <- factor(non_reliable_alleles_text$ALLELES, levels = allele_palette$AlleleCol)

      } else {
        non_reliable_alleles_text <- c()
        allele_palette <- allelePalette(geno2$ALLELES)
      }

      if(any(grepl('^[0-9]+[_][A-Z][0-9]+[A-Z]',geno2$ALLELES))){
      #check novel allele count
      novel <- data.frame(Novel=grep('^[0-9]+[_][A-Z][0-9]+[A-Z]',geno2$ALLELES,value=T),
                          Base = sapply(grep('[A-Z][0-9]+[A-Z]',geno2$ALLELES,value=T),
                                        function(x) strsplit(x,'[_]')[[1]][1])) %>%
        distinct() %>% group_by(.data$Base) %>% mutate(n = dplyr::n())
      # if any base is larger then 6 change to dagger and number
      if(any(novel$n>6)){
        id <- grep('^[0-9]+[_][A-Z][0-9]+[A-Z]',names(allele_palette$transper))
        allele_palette$transper[id] <- 1
        new_allele <- paste0(dagger,1:length(id),'-',allele_palette$AlleleCol[id])
        names(new_allele) <- allele_palette$AlleleCol[id]
        novel_allele_text <- novelAlleleAnnotation_geno(geno2[grep(paste0("^",names(new_allele),collapse = "|"), geno2$ALLELES), ],
                                                   new_label = new_allele, size = 2)
        for(i in 1:length(new_allele)){
          allele <- names(new_allele)[i]
          geno2$ALLELES[grep(allele,geno2$ALLELES)] <- new_allele[i]
        }
        allele_palette$AlleleCol[id] <- new_allele
        names(allele_palette$transper)[id] <- new_allele


      }else{
        novel_allele_text <- c()
      }}else{novel_allele_text <- c()}
      # levels defenitions: allele before Kdiff, D1=levels of allele, D2= levels of Kdiff
      D2 <- levels(factor(kval.df$K_GROUPED))
      D1 <- levels(factor(geno2$ALLELES))

      geno2_kval <- bind_rows(geno2 %>% as.data.frame(), kval.df2)
      # set levels for split legends with titles
      geno2_kval$ALLELES = factor(geno2_kval$ALLELES, levels=c("Allele           ",allele_palette$AlleleCol," ","log(lK)",D2))

      #~~~~~~~~~~~~~~~~~~~~~~~~~ Create graph ~~~~~~~~~~~~~~~~~~~~~~~
      # create specific pallete
      #allele_palette <- allelePalette(D1)
      #colors.set <- c( names(allele_palette$AlleleCol),brewer.pal(length(D2), name="PuBu"))

      # cerate graph
      #full_plot <- ggplot(geno2_kval, aes(x = GENE, y = freq, fill = ALLELES,
      #geom_col(position = "fill") + coord_flip() for y=freq   # geom_bar(position = "fill") + coord_flip()+  without y
      geno2_kval$text2 <- paste0("Individual: ",geno2_kval$SUBJECT,
                                 '<br />',"Gene:", geno2_kval$GENE ,
                                 '<br />',"Allele/lK:", geno2_kval$ALLELES,
                                 '<br />',"Count:", geno2_kval$Freq_by_Clone)

      col_x <- ifelse(levels(geno2_kval$GENE) %in% gsub(chain,"",color_pes_orf), "brown", "black")

      full_plot <- ggplot(geno2_kval,aes(text=text2)) + geom_col(data = geno2_kval,mapping = aes(x = GENE, y=freq, fill = ALLELES),
                                       position = "fill")+
          theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(colour = col_x),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             text = element_text(size = 10), strip.background = element_blank(),
                             strip.text = element_text(face = "bold"),axis.text = element_text(colour = "black")) + coord_flip()+
        xlab("Gene") + ylab("") + scale_fill_hue(name = "Allele", h = c(0, 270), h.start = 10)+
        scale_fill_manual(values=c("white",alpha(names(allele_palette$AlleleCol), allele_palette$transper),
                                       #colors.set[length(colors.set)],  #needs to be the NR colors
                                       "white","white",brewer.pal(length(D2), name="PuBu")),
                            drop=FALSE) +
          guides(fill=guide_legend(ncol =2)) +
          theme(legend.position="right",
                legend.key = element_rect(fill=NA),
                legend.title=element_blank()) + facet_grid(paste0(".~", facet_by)) # split by SAMPLES
      short_reads = F
      if(is.data.frame(non_reliable_alleles_text) | is.data.frame(novel_allele_text)){
        unique_text <- bind_rows(non_reliable_alleles_text, novel_allele_text)
        unique_text$ALLELE_TEXT <- unique_text$ALLELES
        unique_text$ALLELE_TEXT[unique_text$ALLELE_TEXT=="NRA"] <- unique_text$text_bottom[unique_text$ALLELE_TEXT=="NRA"]

        for(i in 1:nrow(unique_text)){

          s = unique_text$SUBJECT[i]
          a =  gsub('\\[[*][0-9]+\\] ',"",unique_text$ALLELE_TEXT[i])
          g = unique_text$GENE[i]

          pos <- geno2_kval %>% filter(SUBJECT==s,GENE==g,ALLELE_TEXT==a) %>% select(freq,n) %>% unlist()

          if(!grepl("[*]",unique_text$ALLELE_TEXT[i])) pos = pos[1]+0.5*pos[2]
          else pos = pos[1]
          unique_text$pos_f[i] <- max(unique_text$pos[i],pos)
        }


        unique_text$text2 <- paste0("Individual: ",unique_text$SUBJECT,
                                    '<br />', "Gene: ", unique_text$GENE,
                                    '<br />',"Allele/lK:", unique_text$text_bottom,
                                    '<br />',"Count: ", unique_text$Freq_by_Clone)

      short_reads = F
      bottom_annot <- c()
      if (is.data.frame(non_reliable_alleles_text)) {

        full_plot <- full_plot + geom_text(data = unique_text[unique_text$ALLELES=="NRA",],
                           aes_string(label = "text", x = "GENE", y = "pos",text = "text2"), angle = 0,
                           size = unique_text$size[unique_text$ALLELES=="NRA"])

        bottom_annot <- unique(non_reliable_alleles_text$text_bottom)
        names(bottom_annot) <- unique(non_reliable_alleles_text$text)
        short_reads = T
      }


      ## multiple novel alleles annot
      if (is.data.frame(novel_allele_text)) {
        full_plot <- full_plot + geom_text(data = unique_text[unique_text$ALLELES!="NRA",],
                                           aes_string(label = "text", x = "GENE", y = "pos", text = "text2"), angle = 0,
                                           size = unique_text$size[unique_text$ALLELES!="NRA"])
      }

      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~ Save as pdf ~~~~~~~~~~~~~~~~~~~~~~~
       if(!html){


         # save lenghts to set pdf proportion
         samples_len <- length(unique(geno2$SUBJECT))
         genes_len <- length(unique(geno2$GENE))
         allele_len <- length(unique(geno2$ALLELES))
         K_len <- length(unique_k)
         height_pdf = genes_len*0.2+samples_len*0.2
         width_pdf = ((2*samples_len)+5)
         ratio = allele_len -genes_len
         ratio2 = (allele_len+K_len)/2 - (height_pdf)
         if(ratio>=0 ){height_pdf = height_pdf+(ratio*0.2) }
         if(ratio2>=0){height_pdf = height_pdf+(ratio2*0.2)}

         if(short_reads){
           bottom_annot <- unique(non_reliable_alleles_text$text_bottom)
           # Create text for annotating the short reads labels
           gp <- ggplotGrob(full_plot)

           # Add 'lab' grob to that row, under the plot panel
           layout = grep('panel',gp$layout$name)

           line_width = (gp$layout[layout[length(layout)],]$r - gp$layout[layout[1],]$l) * 5.42

           annot <- splitlines(bottom_annot,line_width = line_width)
           annot <- annot[!is.na(dplyr::na_if(annot,""))]
           lab <- grid::textGrob(paste0(annot,collapse = '\n')
                                 , x = unit(.1, "npc"), just = c("left"), gp = grid::gpar(fontsize = 8, col = "black"))
           # Add a row below the 2nd from the bottom
           gp <- gtable::gtable_add_rows(gp, unit(1, "grobheight", lab), -2)

           short_reads_rows = length(annot)

           gp <- gtable::gtable_add_grob(gp, lab, t = -2, l = gp$layout[layout[1],]$l, r = gp$layout[layout[length(layout)],]$r)

           full_plot <- plot_grid(gp)
         }


          # save pdf
         pdf(output_file,height = height_pdf, width = width_pdf )
         print(full_plot)
         dev.off()

       }else{
         #~~~~~~~~~~~~~~~~~~~~~~~~~ Save as html ~~~~~~~~~~~~~~~~~~~~~~~
         #p.l.c <- ggplotly(full_plot)
         p.l.c <- ggplotly(full_plot,tooltip = "text2") %>% layout(hovermode='closest', hoverdistance = 10,
                                                                   yaxis = list(title = paste0(c(rep("&nbsp;", 3), "Gene",rep("&nbsp;", 3),
                                                                                                 rep("\n&nbsp;", 1)), collapse = ""))) #%>% config(displayModeBar = F)
         p.l.c$height = length(unique(geno2_kval$GENE)) * 10
         p.l.c$width = ifelse(length(unique(geno2_kval$SUBJECT))>4, "150%", ifelse(length(unique(geno2_kval$SUBJECT))==1, 450, "100%"))

         for(i in grep(dagger,p.l.c$x$data)){
           p.l.c$x$data[[i]]$hoverinfo <- "skip"
         }
         p.l.c$x$layout$xaxis$ticktext <- sapply(p.l.c$x$layout$xaxis$ticktext, function(x){
           ifelse(x %in% gsub(chain,"",color_pes_orf), paste0("[",x,"]"),x)
         },USE.NAMES = F)

         htmlwidgets::saveWidget(p.l.c , file.path(normalizePath(dirname(output_file)),basename(output_file)),selfcontained = F) #VDJbase


       }




}


######################### Run multiGenotype fuction ##########################################
multipleGenoytpe(gen_table = data ,html = html_output, pseudo_genes = pseudo_ORF_genes)




