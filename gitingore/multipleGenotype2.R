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
multipleGenoytpe <- function(geno_table, chain = "IGH", html = FALSE, removeIGH = TRUE, text_size = 14, gene_sort='position', facet_by='SUBJECT', pseudo_genes = F, ORF_genes=F){


      #~~~~~~~~~~~~~~~~~~~~~~~~~~ arrange genotype work file ~~~~~~~~~~~~~~~~~~~~~~~
      # change single number (1,2,3) to allele notation (01,02,03)
      i1 <- grepl("^[0-9]$", gen_table$GENOTYPED_ALLELES)
      gen_table$GENOTYPED_ALLELES[i1] <- paste0("0", gen_table$GENOTYPED_ALLELES[i1])

      # samples
      samples <- unique(geno_table$SUBJECT)
      # select columns
      geno_db <- geno_table %>% select(.data$SUBJECT,.data$GENE,.data$GENOTYPED_ALLELES,.data$K_DIFF,.data$Freq_by_Clone)
      # rename the columns
      names(geno_db)[3:4] <- c("ALLELES", "K")
      # correct deletion annotations
      geno_db$ALLELES <- gsub("Deletion","Del",geno_db$ALLELES)
      # set data.table and correct missing Unk annotations and K
      geno_db <- setDT(geno_db)[CJ(SUBJECT = SUBJECT, GENE = GENE, unique=TRUE), on=.(SUBJECT, GENE)]
      geno_db[is.na(ALLELES) , c("ALLELES","K") := list("Unk", NA_integer_)]
      # set K value for deleted genes
      geno_db$K[grep("Del",geno_db$ALLELES)] <- NA_integer_

      # bin K values and add them as alleles to data
      kval.df <- geno_db
      bins_k <- cut(as.numeric(kval.df$K[!is.na(kval.df$K)]), c(0, 1, 2, 3, 4, 5, 10, 20, 50, Inf), include.lowest = T, right = F)
      K_GROUPED <- gsub(",", ", ", levels(bins_k))
      kval.df$ALLELES[!is.na(kval.df$K)] <- K_GROUPED[bins_k]
      kval.df$ALLELES[is.na(kval.df$K)] <- NA_character_
      kval.df$ALLELES <- factor(kval.df$ALLELES, levels = c(NA_character_, K_GROUPED))

      # bind the two tables
      geno_db <- rbind(geno_db,kval.df)
      # expand row, one allele per row
      geno_db <- splitstackshape::cSplit(geno_db, "ALLELES", sep = ",", direction = "long", fixed = T, type.convert = F)

      # add pseudo genes and orf to color base
      color_pes_orf <- c()
      if(pseudo_genes){
        color_pes_orf <- c(grep("V",PSEUDO[[chain]],value = T),color_pes_orf)
      }
      if(ORF_genes){
        color_pes_orf <- c(unique(grep("OR|NL", genotype$GENE,value = T)),color_pes_orf)
      }

      # sort the data, remove pseudo and orf if needed
      geno_db <- sortDFByGene(DATA = geno_db, chain = chain, method = gene_sort, removeIGH = removeIGH, geno = T,
                              peseudo_remove = pseudo_genes, ORF_remove = ORF_genes)


      geno2 <- as.data.frame(geno2  %>% group_by(.data$SUBJECT, .data$GENE)  %>% mutate(n = dplyr::n()))
      #geno2$freq <- ifelse(geno2$n == 1, 0.75, ifelse(geno2$n == 2, 0.375, ifelse(geno2$n == 3, 0.25, 0.1875)))
      geno2 <- geno2  %>% group_by(.data$GENE, .data$SUBJECT)  %>%
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
                                        function(x) strsplit(x,'[_]')[[1]][1]))  %>%
        distinct()  %>% group_by(.data$Base)  %>% mutate(n = dplyr::n())
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

      geno2_kval <- bind_rows(geno2  %>% as.data.frame(), kval.df2)
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
        unique_text$ALLELE_TEXT <- as.character(unique_text$ALLELES)
        unique_text$ALLELE_TEXT[unique_text$ALLELE_TEXT=="NRA"] <- unique_text$text_bottom[unique_text$ALLELE_TEXT=="NRA"]
        unique_text$pos_f <- unique_text$pos
        for(i in 1:nrow(unique_text)){

          s = unique_text$SUBJECT[i]
          a =  gsub('\\[[*][0-9]+\\] ',"",unique_text$ALLELE_TEXT[i])
          g = unique_text$GENE[i]

          pos <- geno2_kval  %>% filter(SUBJECT==s,GENE==g,ALLELE_TEXT==a)  %>% select(freq,n)  %>% unlist()

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
         samples_len <- length(unique(geno_db$SUBJECT))
         genes_len <- length(unique(geno_db$GENE))
         allele_len <- length(unique(grep('^[0-9]',geno_db$ALLELES)))
         K_len <- length(unique(kval.df$ALLELES))
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


          return(list(full_plot,height_pdf,width_pdf))

       }else{
         #~~~~~~~~~~~~~~~~~~~~~~~~~ Save as html ~~~~~~~~~~~~~~~~~~~~~~~
         #p.l.c <- ggplotly(full_plot)
         p.l.c <- ggplotly(full_plot,tooltip = "text2")  %>% plotly::layout(hovermode='closest', hoverdistance = 10,
                                                                   yaxis = list(title = paste0(c(rep("&nbsp;", 3), "Gene",rep("&nbsp;", 3),
                                                                                                 rep("\n&nbsp;", 1)), collapse = ""))) # %>% config(displayModeBar = F)
         p.l.c$height = length(unique(geno2_kval$GENE)) * 10
         p.l.c$width = ifelse(length(unique(geno2_kval$SUBJECT))>4, "150%", ifelse(length(unique(geno2_kval$SUBJECT))==1, 450, "100%"))

         for(i in grep(dagger,p.l.c$x$data)){
           p.l.c$x$data[[i]]$hoverinfo <- "skip"
         }
         p.l.c$x$layout$xaxis$ticktext <- sapply(p.l.c$x$layout$xaxis$ticktext, function(x){
           ifelse(x %in% gsub(chain,"",color_pes_orf), paste0("[",x,"]"),x)
         },USE.NAMES = F)

         return(list(p.l.c,0,0))


       }




}





