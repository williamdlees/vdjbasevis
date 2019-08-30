library(rabhit)
library(dplyr)
library(optparse)
library(plotly)
library(vdjbaseVis)
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
#load(opt$sysdata_file)

html_output <- as.logical(opt$is_html)  # for pdf set "F"
pseudo_ORF_genes <- as.logical(opt$with_pseudo)

p <- multipleGenoytpe(gen_table = data ,html = html_output, pseudo_genes = pseudo_ORF_genes)

if(html_output){
  htmlwidgets::saveWidget(p[[1]] , file.path(normalizePath(dirname(output_file)),basename(output_file)),selfcontained = F) #VDJbase
}else{
  # save pdf
  pdf(output_file,height = p[[2]], width = p[[3]] )
  print(p[[1]])
  dev.off()
}

