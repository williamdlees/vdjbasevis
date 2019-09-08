library('plyr')
library('magrittr')
library('ggpubr')
library(optparse)
library('purrr')
require('readxl')
library(stringr)
option_list = list(
  make_option(c("-i", "--input_file"), type="character", default=NULL,
              help="excel file name", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="graph.pdf file name", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input_file)){
  stop("input file must be supplied", call.=FALSE)
}

if (is.null(opt$output_file)){
  stop("output reference file must be supplied", call.=FALSE)
}

# read data

data_merge <- readxl::excel_sheets(path) %>% purrr::set_names() %>% purrr::map(read_excel, path = path)

# reshape list to dataframe
require(reshape2)
data_merge$id <- rownames(data_merge)

