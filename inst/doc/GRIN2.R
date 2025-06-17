## ----loadLibs, include = FALSE------------------------------------------------
# Required package:
library(GRIN2)

# Suggested packages: load if available
if (requireNamespace("knitr", quietly = TRUE)) {
  library(knitr)
}

if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)
}

#old<- options()
# Load example datasets from GRIN2 package
data(clin_data)
data(lesion_data)
data(expr_data)
data(hg38_gene_annotation)
data(hg38_chrom_size)
data(hg38_cytoband)
data(pathways)

# Set knitr options
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  echo = TRUE,
  error = TRUE,
  eval = TRUE,
  digits = 3,
  tidy = FALSE,
  background = "#FFFF00",
  fig.align = 'center',
  warning = FALSE,
  message = FALSE
)

# Reset options at end of vignette
oldOpt <- options(width = 55, digits = 3)
on.exit(options(oldOpt), add = TRUE)

# Helper function
getInfo <- function(what = "Suggests") {
  text <- packageDescription("GRIN2")[what][[1]]
  text <- gsub("\n", ", ", text, fixed = TRUE)
  text <- gsub(">=", "$\\\\ge$", text, fixed = TRUE)
  eachPkg <- strsplit(text, ", ", fixed = TRUE)[[1]]
  eachPkg <- gsub(",", "", eachPkg, fixed = TRUE)
  length(eachPkg)
}

## -----------------------------------------------------------------------------
data(clin_data)
data(lesion_data)
data(expr_data)

# Please make sure that coordinates in the lesion data file (loc.start and loc.end) are based on GRCh38 (hg38) genome assembly. Multiple tools that include the UCSC LiftOver tool (<https://genome.ucsc.edu/cgi-bin/hgLiftOver>) can be used for such conversions.

head(lesion_data)

# Specify a folder on your local machine to store the analysis results:
# resultsPath=tempdir()
# knitr::opts_knit$set(root.dir = normalizePath(path = resultsPath))


## ----message=FALSE, warning=FALSE---------------------------------------------

# A subset (417 genes) retrieved from ensembl BioMart will be used as an example gene annotation data file:
data(hg38_gene_annotation)

# To retrieve a complete annotation data for genes and regulatory features from ensembl BioMart, users can use get.ensembl.annotation function:
# hg38.ann <- get.ensembl.annotation("Human_GRCh38")

head(hg38_gene_annotation)


## ----message=FALSE, warning=FALSE---------------------------------------------

# chromosome size data:
head(hg38_chrom_size)

# chromosome size data can be also retrieved for GRCh38 (hg38) genome build from chr.info txt file available on UCSC genome browser:


## ----message=FALSE, warning=FALSE---------------------------------------------

# grin.stats is the function that can be used to excute the GRIN statistical framework:
grin.results=grin.stats(lesion_data, 
                        hg38_gene_annotation, 
                        hg38_chrom_size)

# In addition to annotated genes, users can run GRIN analysis on regulatory sequences from ensembl regulatory build and FANTOM5 project.


## -----------------------------------------------------------------------------
# Extract GRIN results table:
grin.table=grin.results$gene.hits
sorted.results <- grin.table[order(as.numeric(as.character(grin.table$p2.nsubj))),]

## -----------------------------------------------------------------------------
head(sorted.results[,c(7,11:14)])

## -----------------------------------------------------------------------------
head(sorted.results[,c(7,19:22)])

## -----------------------------------------------------------------------------
head(sorted.results[,c(7,27:30)])

## -----------------------------------------------------------------------------
head(sorted.results[,c(7,31:34)])

## -----------------------------------------------------------------------------
# write.grin.xlsx function return an excel file with multiple sheets that include GRIN results table, interpretation of each column in the results, and methods paragraph
# write.grin.xlsx(grin.results, "T-ALL_GRIN_result_annotated_genes.xlsx")

# To return the results table without other information (will be helpful in case of large lesion data files where the gene.lsn.data sheet will be > 1 million rows that halt the write.grin.xlsx function).
grin.res.table=grin.results$gene.hits


## ----genomewide.lesion.plot, fig.height = 7, fig.width = 8, fig.cap="Figure 1. Genome-wide lesion plot"----
genomewide.plot=genomewide.lsn.plot(grin.results, 
                                    max.log10q=50) 
# This function use the list of grin.results


## ----grin.stacked.barplot, fig.height = 10, fig.width = 10, fig.cap="Figure 2. stacked barplot with number of patients affected by different types of lesions in a list of genes of interest"----
# This barplot shows the number of patients affected by different types of lesions in a list of genes of interest:
count.genes=as.vector(c("CDKN2A", "NOTCH1", "CDKN2B", "TAL1", "FBXW7", "PTEN", "IRF8",
                        "NRAS", "BCL11B", "MYB", "LEF1","RB1", "MLLT3", "EZH2", "ETV6",
                        "CTCF", "JAK1", "KRAS", "RUNX1", "IKZF1", "KMT2A", "RPL11", "TCF7",
                        "WT1", "JAK2", "JAK3", "FLT3"))

# return the stacked barplot
grin.barplt(grin.results,
            count.genes)


## ----message=FALSE, warning=FALSE---------------------------------------------
# First identify the list of genes to be included in the oncoprint:
oncoprint.genes=as.vector(c("ENSG00000101307", "ENSG00000171862", "ENSG00000138795", 
                            "ENSG00000139083", "ENSG00000162434", "ENSG00000134371",
                            "ENSG00000118058", "ENSG00000171843", "ENSG00000139687",
                            "ENSG00000184674", "ENSG00000118513", "ENSG00000197888",
                            "ENSG00000111276", "ENSG00000258223", "ENSG00000187266",
                            "ENSG00000174473", "ENSG00000133433", "ENSG00000159216",
                            "ENSG00000107104", "ENSG00000099984", "ENSG00000078403",
                            "ENSG00000183150", "ENSG00000081059", "ENSG00000175354",
                            "ENSG00000164438"))

# Prepare a lesion matrix for the selected list of genes with each row as a gene and each column is a patient (this matrix is compatible with oncoPrint function in ComplexHeatmap package):
oncoprint.mtx=grin.oncoprint.mtx(grin.results, 
                                 oncoprint.genes)

head(oncoprint.mtx[,1:6])

# Use onco.print.props function to specify a height proportion for each lesion category:
onco.props<-onco.print.props(lesion.data,
                             hgt = c("gain"=5, "loss"=4, "mutation"=2, "fusion"=1))
column_title = "" # optional


# use oncoprint function from ComplexHeatmap library to plot the oncoprint:


## ----locus lesion plot, fig.keep='last', fig.height =8, fig.width = 9, fig.cap="Figure 3. lesion plot showing all deletions affecting CDKN2A locus on chromosome 9"----

# First call the data for "hg38_cytoband":
data(hg38_cytoband)

# Plots Showing a Selected Type of Lesion (Ex: deletions) Affecting a region of Interest:
cdkn2a.locus=lsn.transcripts.plot(grin.results, transTrack = FALSE,
                                  hg38.cytoband=hg38_cytoband, chrom=9,                                                      plot.start=19900000, plot.end=25600000,                                                    lesion.grp = "loss", spec.lsn.clr = "blue")

# This type of lesion plots can be also prepared for a certain gene or locus of interest with transcripts track added using the same lsn.transcripts.plot function.


## ----regional.lesion.plot, fig.keep='last', fig.height =10, fig.width = 8, fig.cap="Figure 4. Plots showing all different types of lesions affecting the whole chromosome"----

# Plots Showing Different Types of Lesions Affecting the whole chromosome:
chrom.plot=lsn.transcripts.plot(grin.results, transTrack = FALSE,
                                  hg38.cytoband=hg38_cytoband, chrom=9,                                                      plot.start=1, plot.end=141000000)


## ----message=FALSE, warning=FALSE---------------------------------------------
# Prepare gene and lesion data for later computations
# This lesion matrix has all lesion types that affect a single gene in one row. It can be used to run association analysis with expression data (part of alex.prep.lsn.expr function)

# First step is to prepare gene and lesion data for later computations
gene.lsn=prep.gene.lsn.data(lesion_data,
                            hg38_gene_annotation)
# Then determine lesions that overlap each gene (locus)
gene.lsn.overlap= find.gene.lsn.overlaps(gene.lsn)
# Finally, build the lesion matrix using prep.lsn.type.matrix function: 
gene.lsn.type.mtx=prep.lsn.type.matrix(gene.lsn.overlap,
                                       min.ngrp=5)
# prep.lsn.type.matrix function return each gene in a row, if the gene is affected by multiple types of lesions (for example gain AND mutations), entry will be denoted as "multiple" for this specific patient.
# min.ngrp can be used to specify the minimum number of patients with a lesion to be included in the final lesion matrix.

head(gene.lsn.type.mtx[,1:5])


## ----message=FALSE, warning=FALSE---------------------------------------------
# alex.prep.lsn.expr function prepare expression, lesion data and return the set of genes with both types of data available ordered by gene IDs in rows and patient IDs in columns:

alex.data=alex.prep.lsn.expr(expr_data,
                             lesion_data,
                             hg38_gene_annotation,
                             min.expr=1,
                             min.pts.lsn=5)

# ALEX ordered lesion data:
alex.lsn=alex.data$alex.lsn
head(alex.lsn[,1:5])

# ALEX ordered expression data:
alex.expr=alex.data$alex.expr
head(alex.expr[,1:5])

## ----message=FALSE, warning=FALSE---------------------------------------------
# KW.hit.express function runs Kruskal-Wallis test for association between lesion groups and expression level of the same corresponding gene:

alex.kw.results=KW.hit.express(alex.data,
                               hg38_gene_annotation,
                               min.grp.size=5)


## -----------------------------------------------------------------------------
# order the genes by the ones with most significant KW q-value:
sorted.kw <- alex.kw.results[order(as.numeric(as.character(alex.kw.results$q.KW))),]

## -----------------------------------------------------------------------------
head(sorted.kw[,c(6,7,11,12)])

## -----------------------------------------------------------------------------
head(sorted.kw[,c(13:18)])

## -----------------------------------------------------------------------------
head(sorted.kw[,c(19:24)])

## ----JAK2.waterfall.plot, fig.height =6, fig.width = 6, fig.cap="Figure 5. JAK2 Water-fall plot which  offers a side-by-side graphical representation of lesion and expression data for each patient"----
# waterfall plots allow a side-by-side representation of expression and lesion data of the gene of interest.
# First prepare expression and lesion data for waterfall plots:
WT1.waterfall.prep=alex.waterfall.prep(alex.data,
                                        alex.kw.results,
                                        "WT1",
                                        lesion_data)

# alex.waterfall.plot can be used to return the plot
WT1.waterfall.plot=alex.waterfall.plot(WT1.waterfall.prep,
                                        lesion_data)

## ----lesion.expression.pathway, fig.height =6, fig.width = 6, fig.cap="Figure 6. Ordered Lesion and Expression Data based on the Clustering Analysis on the pathway level (JAK/STAT pathway)"----

data("pathways")

# alex.pathway will return two panels figure of lesion and expression data of ordered subjects based on the computed lesions distance in all genes assigned to the pathway of interest:
alex.path=alex.pathway(alex.data,
                       lsn.data = lesion_data,
                       pathways = pathways, 
                       selected.pathway = "Jak_Pathway")

# To return ordered lesion and expression data of the genes assigned to the pathway of interest (same patients order in the plot):
alex.path[1:10,1:5]


## ----message=FALSE, warning=FALSE---------------------------------------------
# This type of lesion matrices with each gene affected by a certain type of lesion in a separate row is very helpful to run multiple levels of association analysis that include association between lesions and treatment outcomes.

# Users should first Prepare gene and lesion data and determine lesions that overlap each gene (locus):
gene.lsn=prep.gene.lsn.data(lesion_data,
                            hg38_gene_annotation)    
gene.lsn.overlap= find.gene.lsn.overlaps(gene.lsn)

# use prep.binary.lsn.mtx function to prepare the lesion binary matrix:
lsn.binary.mtx.atleast5=prep.binary.lsn.mtx(gene.lsn.overlap,
                                            min.ngrp=5)
# Each row is a lesion type that affect a certain gene for example NOTCH1_mutation (entry will be labelled as 1 if the patient is affected by by this type of lesion and 0 otherwise).
# min.ngrp can be used to specify the minimum number of patients with a lesion to be included in the final lesion matrix.

head(lsn.binary.mtx.atleast5[,1:5])


## ----message=FALSE, warning=FALSE---------------------------------------------
# Prepare Event-free Survival (EFS) and Overall Survival (OS) as survival objects:
clin_data$EFS <- Surv(clin_data$efs.time, clin_data$efs.censor)
clin_data$OS <- Surv(clin_data$os.time, clin_data$os.censor)

# List all clinical variables of interest to be included in the association analysis:
clinvars=c("MRD.binary", "EFS", "OS")

# Run association analysis between lesions and clinical variables:
assc.outcomes=grin.assoc.lsn.outcome(lsn.binary.mtx.atleast5,
                                     clin_data,
                                     hg38_gene_annotation,
                                     clinvars)

# Optional: Adjust for covariates using the 'covariate' argument

## -----------------------------------------------------------------------------
# order the genes by the ones with most significant KW q-value:
sorted.outcomes <- assc.outcomes[order(as.numeric(as.character(assc.outcomes$`MRD.binary.p-value`))),]

## -----------------------------------------------------------------------------
head(sorted.outcomes[1:7,c(6,11,14,15)])

## -----------------------------------------------------------------------------
head(sorted.outcomes[1:7,c(6, 16:19)])

## ----setup--------------------------------------------------------------------
library(GRIN2)

