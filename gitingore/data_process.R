

# load('test_func.RData')
# load('R/sysdata.rda')
# chain = "IGH"; gene_sort = "position"; removeIGH = TRUE; lk_cutoff = 1; mark_low_lk = TRUE; html = FALSE; n_line = 4; line_length=120; pseudo_genes = FALSE; ORF_genes = FALSE; size = 3
samples <- unique(geno_table$SUBJECT)
geno_db <- geno_table %>% select(.data$SUBJECT,.data$GENE,.data$GENOTYPED_ALLELES,.data$K_DIFF)
names(geno_db)[3:4] <- c("ALLELES", "K")
geno_db$ALLELES <- gsub("Deletion","Del",geno_db$ALLELES)
geno_db <- splitstackshape::cSplit(geno_db, "ALLELES", sep = ",", direction = "long", fixed = T, type.convert = F)
geno_db <- setDT(geno_db)[CJ(SUBJECT = SUBJECT, GENE = GENE, unique=TRUE), on=.(SUBJECT, GENE)]
geno_db[is.na(ALLELES) , c("ALLELES","K") := list("Unk", NA_integer_)]


geno_db$K[grep("Del",geno_db$ALLELES)] <- NA_integer_

color_pes_orf <- c()
if(pseudo_genes){
  color_pes_orf <- c(grep("V",PSEUDO[[chain]],value = T),color_pes_orf)
}
if(ORF_genes){
  color_pes_orf <- c(unique(grep("OR|NL", geno_db$GENE,value = T)),color_pes_orf)
}

geno_db <- sortDFByGene(DATA = geno_db, chain = chain, method = gene_sort, removeIGH = removeIGH, geno = T,
                        peseudo_remove = pseudo_genes, ORF_remove = ORF_genes)

#geno_db <- geno_db[geno_db$GENE %in% unique(geno_db$GENE)[1:15]]
#geno_db <- geno_db[geno_db$SUBJECT %in% unique(geno_db$SUBJECT)[1:15]]
geno_db$GENE <- factor(geno_db$GENE, levels = gsub("IG[H|K|L]", "", GENE.loc[[chain]]))
gene_loc <- 1:length(unique(geno_db$GENE)[order(match(unique(geno_db$GENE), levels(geno_db$GENE)))])
names(gene_loc) <- unique(geno_db$GENE)[order(match(unique(geno_db$GENE), levels(geno_db$GENE)))]
geno_db$GENE_LOC <- gene_loc[as.character(geno_db$GENE)]

gene_num <- round(length(unique(geno_db$GENE))/3)

geno_db_texture <- c()

geno_db_texture = geno_db[K<lk_cutoff & !ALLELES %in% c("Unk", "Del", "NR")][][CJ(SUBJECT = SUBJECT, GENE = GENE, n_line = 1:n_line, unique = TRUE), on = .(SUBJECT, GENE), nomatch = 0]
if(mark_low_lk & nrow(geno_db_texture) != 0){
  geno_db_texture[,n_line:=NULL]
  geno_db_texture[, c("points", "yend", "x", "xend") := list(seq(0, 0.9, length.out = n_line),
                                                             seq(0, 0.9, length.out = n_line) + 0.1,
                                                             gene_loc[as.character(GENE)] - 0.49,
                                                             gene_loc[as.character(GENE)] + 0.49),
                  by = list(SUBJECT,GENE,ALLELES)]
}else{
  mark_low_lk <- F
}
#sort the heatmap for plotting
geno_db_m <- geno_db[, n:=  .N, by = list(SUBJECT,GENE)][]
geno_db_m[, freq:= 1/n]
geno_db_m[, title := "Genotype"]
geno_db_m$ALLELES_G <- geno_db_m$ALLELES #for grouping
geno_db_m$text <- ''
geno_db_m$text_bottom <- geno_db_m$ALLELES
# change ambiguous alleles call
id_nra <- grepl("[0-9][0-9]_[0-9][0-9]", geno_db_m$ALLELES)
nra <- F
if (any(id_nra)) {
  num_text <- paste0('[*',1:length(unique(geno_db_m$ALLELES[id_nra])),']')
  names(num_text) <- unique(geno_db_m$ALLELES[id_nra])
  geno_db_m$text[id_nra] <- num_text[geno_db_m$ALLELES[id_nra]]
  geno_db_m$text_bottom[id_nra] <- paste(num_text[geno_db_m$ALLELES[id_nra]],geno_db_m$ALLELES[id_nra])
  geno_db_m$ALLELES[id_nra] <- "NRA"
  nra <- T
}
allele_palette <- allelePalette(geno_db_m$ALLELES)

val_novel <- grep('^[0-9]+[_][A-Z][0-9]+[A-Z]',geno_db_m$ALLELES, value = T)
novel <- F
novel_allele_text <- c()
if(length(val_novel)!=0){
  #check novel allele count
  # novel <- data.frame(Novel=val_novel,
  #                     Base = sapply(val_novel,function(x) strsplit(x,'[_]')[[1]][1])) %>%
  #   distinct() %>% group_by(.data$Base) %>% mutate(n = dplyr::n())
  # if any base is larger then 6 change to dagger and number
  dagger = "§"
  # if(any(novel$n>6)){
    id <- grep('^[0-9]+[_][A-Z][0-9]+[A-Z]',names(allele_palette$transper))
    allele_palette$transper[id] <- 1
    code_allele <- paste0(dagger,1:length(id))
    names(code_allele) <-allele_palette$AlleleCol[id]
    new_allele <- paste0(dagger,1:length(id),'-',allele_palette$AlleleCol[id])
    names(new_allele) <-allele_palette$AlleleCol[id]
    #novel_allele_text <- novelAlleleAnnotation(geno_db_m[grep(paste0("^",names(new_allele),collapse = "|"), geno_db_m$ALLELES), ],
    #                                           new_label = new_allele)
    ids <- geno_db_m$ALLELES %fin% names(new_allele)
    rep <- new_allele[geno_db_m$ALLELES[ids]]
    rep2 <- code_allele[geno_db_m$ALLELES[ids]]
    geno_db_m[ids, c("ALLELES","text_bottom","text") := list(rep,rep,rep2)]

    allele_palette$AlleleCol[id] <- new_allele
    names(allele_palette$transper)[id] <- new_allele
    novel <- T
}
# }

geno_db_m$ALLELES <- factor(geno_db_m$ALLELES, levels = allele_palette$AlleleCol)
