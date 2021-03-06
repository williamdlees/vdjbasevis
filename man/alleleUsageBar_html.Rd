% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alleleUsageBar_html.R
\name{alleleUsageBar_html}
\alias{alleleUsageBar_html}
\title{Graphical output of allele usage in given population}
\usage{
alleleUsageBar_html(gene_segment, chain = c("IGH", "IGK", "IGL"))
}
\arguments{
\item{gene_segment}{a data frame of allele usage count in a given population. See details.}

\item{chain}{the IG chain: IGH,IGK,IGL. Default is IGH.}
}
\value{
An interactive stacked barplot visualization of the allele hetrouzygousity in a given population.
}
\description{
The \code{alleleUsageBar_html} function generates an interactive graphical output of sum of unique allele that appeared in the given population.
}
\details{
A \code{data.frame} with the following columns.
\itemize{
  \item \code{'GENE'}:     The gene call
  \item \code{'COUNT'}:    Number of unique alleles that appeared in the given population.
}
}
\examples{
 gene_segment <- data.frame(GENE = c("V1-2",'V3-3','D2-8','D3-16','J4','J6'), COUNT = c(2,9,2,2,1,2))
 alleleUsageBar_html(gene_segment)
}
