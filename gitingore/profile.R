library(reshape2)
library(dplyr)
library(ggplot2)
library(alakazam)
library(tigger)
library(seqinr)
library(cowplot)
library(gtools)
library(ggsignif)
library(gtools)
library(ggdendro)
library(stringr)
library(gridExtra)
library(grid)
library(mltools)
library(data.table)
library(lattice)
library(fastmatch) #%fin%
library(splitstackshape)
Rprof(tmp <- tempfile())
non_reliable_alleles_text <- non_reliable_alleles_text %>% ungroup() %>% mutate(s = 0.5/.data$n, by = 1/.data$n) %>%
  group_by(.data$GENE, .data$SUBJECT) %>%
  mutate(pos2 = seq.int(s[1],1,by = by[1]))
Rprof()
summaryRprof(tmp)

source('proftable.R')
Rprof("profile1.out", line.profiling=TRUE)
source('check_code.R')
Rprof(NULL)
summaryRprof("profile1.out", lines = "show")

proftable("profile1.out", lines=10)



library(rbenchmark)
library(splitstackshape)

benchmark(
  "separate_rows" = {
    tidyr::separate_rows(geno_db, "ALLELES", sep = ",")
  },
  "concat.split.multiple" = {
    splitstackshape::cSplit(geno_db, "ALLELES", sep = ",", direction = "long", fixed = T, type.convert = F)
  },
  "cSplit" = {
    splitstackshape::cSplit(geno_db, "ALLELES", sep = ",", direction = "long")
  },
  replications = 10,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self"))



benchmark(
  "grepl" = {
    gt = ggplotGrob(ggplot() + geom_col(data = heatmap.df, mapping = aes_string(x = "GENE_LOC", y = "freq", fill = "ALLELES"),
                                        position = "fill", width = 0.95, na.rm = T) +
                      scale_fill_manual(values = alpha(names(allele_palette$AlleleCol), allele_palette$transper), name = "Alleles", drop = FALSE))
  },
  "fin" = {
    gt1 = ggplotGrob(p)
    },
  replications = 1,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self"))


library(microbenchmark)
check_for_equal_coefs <- function(values) {
  all.equal(values[[1]],values[[2]], check.attributes = FALSE)
}
microbenchmark(
  # "separate_rows" = {
  #   l <- tidyr::separate_rows(geno_db, "ALLELES", sep = ",")
  # },
  "concat.split.multiple" = {
    l <- splitstackshape::cSplit(geno_db, "ALLELES", sep = ",", direction = "long", fixed = T, type.convert = F)
    l <- as.data.frame(l)
  },
  "cSplit" = {
    l <- splitstackshape::cSplit(geno_db, "ALLELES", sep = ",", direction = "long")
    l <- as.data.frame(l)
    l$ALLELES <- as.character(l$ALLELES)
  },
  check = check_for_equal_coefs
)
