% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep.gene.lsn.data.R
\name{prep.gene.lsn.data}
\alias{prep.gene.lsn.data}
\title{Prepare Gene and Lesion Data for GRIN Analysis}
\usage{
prep.gene.lsn.data(lsn.data, gene.data, mess.freq = 10)
}
\arguments{
\item{lsn.data}{A \code{data.frame} containing lesion data in GRIN-compatible format. Must include the following five columns:
\describe{
\item{ID}{Unique patient identifier.}
\item{chrom}{Chromosome on which the lesion is located.}
\item{loc.start}{Start position of the lesion in base pairs.}
\item{loc.end}{End position of the lesion in base pairs.}
\item{lsn.type}{Type of lesion (e.g., gain, loss, mutation, fusion, etc...).}
}}

\item{gene.data}{A \code{data.frame} containing gene annotation data with the following four required columns:
\describe{
\item{gene}{Ensembl gene ID.}
\item{chrom}{Chromosome on which the gene is located.}
\item{loc.start}{Start position of the gene in base pairs.}
\item{loc.end}{End position of the gene in base pairs.}
}}

\item{mess.freq}{Integer specifying the frequency at which progress messages are displayed. Messages are printed every \code{mess.freq}-th lesion block processed (default is 10).}
}
\value{
A list with the following components:
\describe{
\item{lsn.data}{Original lesion data.}
\item{gene.data}{Original gene annotation data.}
\item{gene.lsn.data}{Combined and ordered data.frame of gene and lesion intervals. The \code{cty} column encodes position type: 1 = gene start, 2 = lesion start, 3 = lesion end, 4 = gene end.}
\item{gene.index}{Index data.frame indicating the start and end rows for each chromosome within \code{gene.lsn.data} for genes.}
\item{lsn.index}{Index data.frame indicating the start and end rows for each lesion (grouped by type, chromosome, and subject) within \code{gene.lsn.data}.}
}
}
\description{
Prepares and indexes gene and lesion data for downstream GRIN (Genomic Random Interval) analysis. This function merges and orders gene and lesion coordinates to support efficient computation of overlaps between genes and all different types of genomic lesions (structural or sequence lesions).
}
\details{
This function performs pre-processing by ordering and indexing both gene and lesion data. It combines gene and lesion coordinates into a unified structure, marking each with a specific code (\code{cty}) that identifies whether the row represents a gene or lesion. This merged data is then used in the \code{find.gene.lsn.overlaps()} function to detect gene-lesion overlaps.
}
\examples{
data(lesion_data)
data(hg38_gene_annotation)

# Prepare gene and lesion data for GRIN analysis:
prep.gene.lsn <- prep.gene.lsn.data(lesion_data, hg38_gene_annotation)
}
\references{
Pounds, S., et al. (2013). A genomic random interval model for statistical analysis of genomic lesion data.
Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
}
\seealso{
\code{\link{order.index.gene.data}}, \code{\link{order.index.lsn.data}}, \code{\link{find.gene.lsn.overlaps}}
}
\author{
Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
}
