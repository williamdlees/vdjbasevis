#' @include vdjbaseVis.R
#' @include internal_functions.R
NULL

########################################################################################################
#' Creates genotype graphical output division by subjects in pdf/html format
#'
#' \code{multipleGenoytpe} Takes a genoytpe table and returns genotype graphical output.
#'
#' @param   gen_table          a table of the genoytpes.
#' @param   html               character "F" or "T". "T" for html and "F" for pdf output.
#'
#' @return   genotype graphical output
#' @export
multipleGenoytpe <- function(gen_table, chain = "IGH", html = FALSE, removeIGH = TRUE, text_size = 14, gene_sort='position', facet_by='subject', pseudo_genes = F, ORF_genes=F){


      #~~~~~~~~~~~~~~~~~~~~~~~~~~ arrange genotype work file ~~~~~~~~~~~~~~~~~~~~~~~
      # change single number (1,2,3) to allele notation (01,02,03)
      i1 <- grepl("^[0-9]$", gen_table$GENOTYPED_ALLELES)
      gen_table$GENOTYPED_ALLELES[i1] <- paste0("0", gen_table$GENOTYPED_ALLELES[i1])
      # Replace Deletion annotation
      gen_table$GENOTYPED_ALLELES <- gsub("Deletion","Del",gen_table$GENOTYPED_ALLELES)

      # select columns
      genotype_all <- gen_table[,c("subject","gene","alleles","counts","total","k_diff","GENOTYPED_ALLELES","Freq_by_Clone")]

      # change columns names Genotyped_alleles to Alleles
      genotype_all <- genotype_all  %>% dplyr::rename(temp=GENOTYPED_ALLELES)
      genotype_all<- genotype_all %>% dplyr::rename(GENOTYPED_ALLELES=alleles)
      genotype_all <- genotype_all %>% dplyr::rename(alleles=temp)
      genotype<-genotype_all

      # # for filter subjects or genes!!!!!!!!!!
      #genotype <- subset(genotype, subject %in% unique(genotype$subject)[1:15])
      # genotype <- subset(genotype,  gene%in% unique(genotype_all$gene)[] )
      # # end: for filter subjects or genes!!!!!!!!!!!!!!!

      color_pes_orf <- c()
      if(pseudo_genes){
        color_pes_orf <- c(grep("V",PSEUDO[[chain]],value = T),color_pes_orf)
      }
      if(ORF_genes){
        color_pes_orf <- c(unique(grep("OR|NL", genotype$gene,value = T)),color_pes_orf)
      }

      if(length(chain)!=1){
        genotype <- c()
      for(ch in chain){
          genotype <- rbind(genotype, sortDFByGene(DATA = genotype, chain = ch, method = gene_sort, removeIGH = removeIGH, geno = T,
                              peseudo_remove = pseudo_genes, ORF_remove = ORF_genes))
      }
      }else{
        genotype <- sortDFByGene(DATA = genotype, chain = ch, method = gene_sort, removeIGH = removeIGH, geno = T,
                                                 peseudo_remove = pseudo_genes, ORF_remove = ORF_genes)
      }
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
      alleles = strsplit(as.character(genotype$alleles), ",")
      count_allele = strsplit(as.character(genotype$Freq_by_Clone), ",")
      geno2 = genotype
      r = 1
      for (g in 1:nrow(genotype)) {
        for (a in 1:length(alleles[[g]])) {
          geno2[r, ] = genotype[g, ]
          geno2[r, ]$alleles = alleles[[g]][a]
          geno2[r, ]$Freq_by_Clone = count_allele[[g]][a]
          r = r + 1
        }
      }

      # Omit rows containing specific alleles/k_diff/GENE column of NA
      geno2  <-geno2   %>% tidyr::drop_na(alleles)
      geno2  <-geno2   %>% tidyr::drop_na(k_diff)
      geno2  <-geno2   %>% tidyr::drop_na(gene)
      #geno2  <-geno2   %>% drop_na(Freq_by_Clone)
      genotype <- genotype   %>% tidyr::drop_na(gene)

      #~~~~~~~~~~~~~~~~~~~~~~~~~~ save and group K values ~~~~~~~~~~~~~~~~~~~~~~~
      kval.df = subset(genotype,select=c(subject,gene,GENOTYPED_ALLELES,k_diff,Freq_by_Clone)) # Kval per gene

      # genes <- unique(kval.df$GENE)
      # for(s in unique(kval.df$subject)){
      #   sub <- kval.df[kval.df$subject==s,]
      #   # get missing genes
      #   id <- which(!genes %in% sub$GENE)
      #   if(length(id) != 0) kval.df <- bind_rows(kval.df,data.frame(subject = s, GENE = genes[id], k_diff = 1000,
      #                                                                 GENOTYPED_ALLELES = "Unk",
      #                                                                 Freq_by_Clone = NA, stringsAsFactors = F))
      # }

      #kval.df = subset(geno2,select=c(subject,GENE,GENOTYPED_ALLELES,k_diff,Freq_by_Clone)) # Kval per gene
      # categorize to groups
      kval.df$K_GROUPED <- bin_data(kval.df$K, bins=c(0, 1,2,3,4,5,10,20,50,Inf), binType = "explicit")
      kval.df <- kval.df[!is.na(na_if(kval.df$GENOTYPED_ALLELES,"")),]
      #~~~~~~~~~~~~~~~~~~~~~~~~~ select columns ~~~~~~~~~~~~~~~~~~~~~~~
      geno2 <- geno2   %>% dplyr::select(subject,gene,alleles,k_diff,GENOTYPED_ALLELES,Freq_by_Clone)
      kval.df2 <- kval.df   %>% dplyr::rename(alleles=K_GROUPED)
      unique_k <- unique(kval.df$K_GROUPED)

      # # define ratio between Kdiff to allele appearance, 3/4 [KK 0102 0102 0102]
      #geno2_kval <- rbind(geno2,geno2,geno2,kval.df2)
      genes <- unique(geno2$gene)
      for(s in unique(geno2$subject)){
        sub <- geno2[geno2$subject==s,]
        # get missing genes
        id <- which(!genes %in% sub$gene)
        if(length(id) != 0) geno2 <- bind_rows(geno2,data.frame(subject = s, gene = genes[id], alleles = 'Unk', k_diff = 1000,
                                                                GENOTYPED_ALLELES = "Unk",
                                                                Freq_by_Clone = NA, stringsAsFactors = F))
      }



      geno2 <- as.data.frame(geno2  %>% dplyr::group_by(.data$subject, .data$gene)  %>% dplyr::mutate(n = dplyr::n()))
      #geno2$freq <- ifelse(geno2$n == 1, 0.75, ifelse(geno2$n == 2, 0.375, ifelse(geno2$n == 3, 0.25, 0.1875)))
      geno2 <- geno2  %>% dplyr::group_by(.data$gene, .data$subject)  %>%
        dplyr::mutate(freq = ifelse(.data$n == 1, 0.75,
                              ifelse(.data$n == 2, rep(0.375,2),
                                     ifelse(.data$n == 3, rep(0.25,3),
                                            rep(0.1875,4)))))
      kval.df2$freq <- 0.25
      geno2$ALLELE_TEXT <- geno2$alleles
      dagger = "\u005E"
      if (length(grep("[0-9][0-9]_[0-9][0-9]", geno2$alleles)) != 0) {
        non_reliable_alleles_text <- nonReliableAllelesText_V2(non_reliable_alleles_text = geno2[grep("[0-9][0-9]_[0-9][0-9]", geno2$alleles), ],
                                                               map = T, size = 2)
        geno2$alleles[grep("[0-9][0-9]_[0-9][0-9]", geno2$alleles)] <- "NRA"
        allele_palette <- allelePalette(geno2$alleles)
        non_reliable_alleles_text$alleles <- factor(non_reliable_alleles_text$alleles, levels = allele_palette$AlleleCol)

      } else {
        non_reliable_alleles_text <- c()
        allele_palette <- allelePalette(geno2$alleles)
      }

      if(any(grepl('^[0-9]+[_][A-Z][0-9]+[A-Z]',geno2$alleles))){
        id <- grep('^[0-9]+[_][A-Z][0-9]+[A-Z]',names(allele_palette$transper))
        allele_palette$transper[id] <- 1
        new_allele <- paste0(dagger,1:length(id),'-',allele_palette$AlleleCol[id])
        names(new_allele) <- allele_palette$AlleleCol[id]
        novel_allele_text <- novelAlleleAnnotation_geno(geno2[grep(paste0("^",names(new_allele),collapse = "|"), geno2$alleles), ],
                                                   new_label = new_allele, size = 2)
        for(i in 1:length(new_allele)){
          allele <- names(new_allele)[i]
          geno2$alleles[grep(allele,geno2$alleles)] <- new_allele[i]
        }
        allele_palette$AlleleCol[id] <- new_allele
        names(allele_palette$transper)[id] <- new_allele
      }else{novel_allele_text <- c()}
      # levels defenitions: allele before Kdiff, D1=levels of allele, D2= levels of Kdiff
      D2 <- levels(factor(kval.df$K_GROUPED))
      D1 <- levels(factor(geno2$alleles))

      geno2_kval <- bind_rows(geno2  %>% as.data.frame(), kval.df2)
      # set levels for split legends with titles
      geno2_kval$alleles = factor(geno2_kval$alleles, levels=c("Allele           ",allele_palette$AlleleCol," ","log(lK)",D2))

      #~~~~~~~~~~~~~~~~~~~~~~~~~ Create graph ~~~~~~~~~~~~~~~~~~~~~~~
      # create specific pallete
      #allele_palette <- allelePalette(D1)
      #colors.set <- c( names(allele_palette$AlleleCol),brewer.pal(length(D2), name="PuBu"))

      # cerate graph
      #full_plot <- ggplot(geno2_kval, aes(x = GENE, y = freq, fill = ALLELES,
      #geom_col(position = "fill") + coord_flip() for y=freq   # geom_bar(position = "fill") + coord_flip()+  without y
      geno2_kval$text2 <- paste0("Individual: ",geno2_kval$subject,
                                 '<br />',"Gene:", geno2_kval$gene ,
                                 '<br />',"Allele/lK:", geno2_kval$alleles,
                                 '<br />',"Count:", geno2_kval$Freq_by_Clone)

      col_x <- ifelse(levels(geno2_kval$gene) %in% gsub(chain,"",color_pes_orf), "brown", "black")

      full_plot <- ggplot(geno2_kval,aes(text=text2)) + geom_col(data = geno2_kval,mapping = aes(x = gene, y=freq, fill = alleles),
                                       position = "fill")+
          theme_bw() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(),axis.text.y = element_text(colour = col_x),
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             text = element_text(size = 10), strip.background = element_blank(),
                             strip.text = element_text(face = "bold"),axis.text = element_text(colour = "black")) + coord_flip()+
        xlab("Gene") + ylab("") + scale_fill_hue(name = "Allele", h = c(0, 270), h.start = 10)+
        scale_fill_manual(values=c("white",alpha(names(allele_palette$AlleleCol), allele_palette$transper),
                                       #colors.set[length(colors.set)],  #needs to be the NR colors
                                       "white","white",RColorBrewer::brewer.pal(length(D2), name="PuBu")),
                            drop=FALSE) +
          guides(fill=guide_legend(ncol =2)) +
          theme(legend.position="right",
                legend.key = element_rect(fill=NA),
                legend.title=element_blank()) + facet_grid(paste0(".~", facet_by)) # split by SAMPLES
      short_reads = F
      if(is.data.frame(non_reliable_alleles_text) | is.data.frame(novel_allele_text)){
        tmp <- geno2_kval %>% dplyr::group_by(subject,gene) %>% dplyr::arrange(desc(alleles)) %>% dplyr::mutate(order = 1:dplyr::n())
        unique_text <- bind_rows(non_reliable_alleles_text, novel_allele_text)
        unique_text$ALLELE_TEXT <- as.character(unique_text$alleles)
        if(length(unique_text$text_bottom[unique_text$ALLELE_TEXT=="NRA"])>0){
        unique_text$ALLELE_TEXT[unique_text$ALLELE_TEXT=="NRA"] <- unique_text$text_bottom[unique_text$ALLELE_TEXT=="NRA"]}

        unique_text <- unique_text %>% dplyr::left_join(tmp[,c('subject','gene','ALLELE_TEXT','order')], by = c('subject','gene','alleles'='ALLELE_TEXT'))
        unique_text$pos_f <- unique_text$freq * (unique_text$order-1) + unique_text$freq/2

        unique_text$text2 <- paste0("Individual: ",unique_text$subject,
                                    '<br />', "Gene: ", unique_text$gene,
                                    '<br />',"Allele/lK:", unique_text$text_bottom,
                                    '<br />',"Count: ", unique_text$Freq_by_Clone)

      short_reads = F
      bottom_annot <- c()
      if (is.data.frame(non_reliable_alleles_text)) {

        full_plot <- full_plot + geom_text(data = unique_text[unique_text$alleles=="NRA",],
                           aes_string(label = "text", x = "gene", y = "pos_f",text = "text2"), angle = 0,
                           size = unique_text$size[unique_text$alleles=="NRA"])

        bottom_annot <- unique(non_reliable_alleles_text$text_bottom)
        names(bottom_annot) <- unique(non_reliable_alleles_text$text)
        short_reads = T
      }


      ## multiple novel alleles annot
      if (is.data.frame(novel_allele_text)) {
        full_plot <- full_plot + geom_text(data = unique_text[unique_text$alleles!="NRA",],
                                           aes_string(label = "text", x = "gene", y = "pos_f", text = "text2"), angle = 0,
                                           size = unique_text$size[unique_text$alleles!="NRA"])
      }

      }

      #~~~~~~~~~~~~~~~~~~~~~~~~~ Save as pdf ~~~~~~~~~~~~~~~~~~~~~~~
       if(!html){


         # save lenghts to set pdf proportion
         samples_len <- length(unique(geno2$subject))
         genes_len <- length(unique(geno2$gene))
         allele_len <- length(unique(geno2$alleles))
         K_len <- length(unique_k)
         height_pdf = genes_len*0.2+samples_len*0.2
         width_pdf = ((2*samples_len)+10)
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


          return(list(full_plot,height_pdf,width_pdf))

       }

}





