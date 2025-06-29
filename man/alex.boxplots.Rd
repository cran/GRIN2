% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alex.boxplots.R
\name{alex.boxplots}
\alias{alex.boxplots}
\title{Generate Box Plots of Gene Expression by Lesion Groups}
\usage{
alex.boxplots(out.dir, alex.data, alex.kw.results, q, gene.annotation)
}
\arguments{
\item{out.dir}{Path to the directory where the resulting PDF file containing the box plots will be saved.}

\item{alex.data}{A list returned by the \code{alex.prep.lsn.expr} function, containing three components: \code{"row.mtch"}, \code{"alex.expr"} (gene expression matrix), and \code{"alex.lsn"} (lesion matrix). Rows of \code{alex.expr} and \code{alex.lsn} is ordered by Ensembl gene IDs; columns represent patient IDs.}

\item{alex.kw.results}{Output of the \code{KW.hit.express} function containing Kruskal Wallis test results.}

\item{q}{Q-value threshold. Only genes with q-values below this threshold will be included in the output plots.}

\item{gene.annotation}{A \code{data.frame} containing gene annotation data, either provided by the user or retrieved using \code{get.ensembl.annotation}. Must contain four columns: \code{"gene"} (Ensembl gene ID), \code{"chrom"} (chromosome), \code{"loc.start"} (start position), and \code{"loc.end"} (end position).}
}
\value{
A PDF file saved in \code{out.dir}, containing one box plot per page for each selected gene, showing gene expression distribution across lesion groups.
}
\description{
Generates box plots of gene expression levels stratified by lesion groups for a subset of genes selected based on a user-specified q-value threshold from Kruskal Wallis test results.
}
\examples{
data(expr_data)
data(lesion_data)
data(hg38_gene_annotation)

# Prepare expression and lesion data
alex.data <- alex.prep.lsn.expr(expr_data, lesion_data,
                                hg38_gene_annotation, min.expr = 5, min.pts.lsn = 5)

# Run Kruskal Wallis test
alex.kw.results <- KW.hit.express(alex.data, hg38_gene_annotation, min.grp.size = 5)

# Generate box plots for significant genes
dir.create(resultsFolder <- file.path(tempdir(), "temp.out"))
alex.boxplots(out.dir = resultsFolder,
              alex.data = alex.data,
              alex.kw.results = alex.kw.results,
              q = 1e-15,
              gene.annotation = hg38_gene_annotation)
unlink(resultsFolder, recursive = TRUE)
}
\references{
Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
}
\seealso{
\code{\link{alex.prep.lsn.expr}}, \code{\link{KW.hit.express}}
}
\author{
Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}
}
