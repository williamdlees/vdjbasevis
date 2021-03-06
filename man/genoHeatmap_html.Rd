% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genoHeatmap_html.R
\name{genoHeatmap_html}
\alias{genoHeatmap_html}
\title{Graphical output of alleles division by genotypes}
\usage{
genoHeatmap_html(
  geno_table,
  chain = c("IGH", "IGK", "IGL"),
  gene_sort = "position",
  removeIGH = TRUE,
  lk_cutoff = 1,
  mark_low_lk = TRUE,
  n_line = 4,
  line_length = 60,
  pseudo_genes = FALSE,
  ORF_genes = FALSE,
  file = file.path(normalizePath(tempdir()), "genotype_heatmap.html")
)
}
\arguments{
\item{geno_table}{genoytpe summary table. See details.}

\item{chain}{the IG chain: IGH,IGK,IGL. Default is IGH.}

\item{gene_sort}{if by 'name' the genes in the output are ordered lexicographically,
if by 'position' only functional genes are used and are ordered by their chromosomal location. Default is 'position'.}

\item{removeIGH}{if TRUE, 'IGH'\'IGK'\'IGL' prefix is removed from gene names.}

\item{lk_cutoff}{the lK cutoff value to be considerd low for texture layer. Defualt is lK<1.}

\item{mark_low_lk}{if TRUE, a texture is add for low lK values. Defualt is TRUE.}
}
\value{
An interactive heat-map visualization of the genotype inference for multiple samples.
}
\description{
The \code{genoHeatmap_html} function generates an interactive graphical output of the alleles per gene in multiple samples.
}
\details{
A \code{data.frame} created by \code{inferGenotypeBaysian}.
}
