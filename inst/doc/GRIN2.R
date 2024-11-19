## ----loadLibs, include = FALSE------------------------------------------------
library(GRIN2)
library(knitr)
library(ComplexHeatmap)
library(survival)
old<- options()
data(clin.data)
data(lesion.data)
data(expr.data)
data(hg19.gene.annotation)
data(kegg.ml.gsets)
data(hg19.chrom.size)
data(hg19_cytoband)
data(hg38_cytoband)
data(pathways)
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  echo=T,
  error=T,
  eval=T,
  digits = 3,
  tidy = FALSE,
  background = "#FFFF00",
  fig.align = 'center',
  warning = FALSE,
  message = FALSE
  )
oldOpt <-options(width = 55, digits = 3)
on.exit(options(oldOpt), add = TRUE)
getInfo <- function(what = "Suggests") {
  text <- packageDescription("GRIN2")[what][[1]]
  text <- gsub("\n", ", ", text, fixed = TRUE)
  text <- gsub(">=", "$\\\\ge$", text, fixed = TRUE)
  eachPkg <- strsplit(text, ", ", fixed = TRUE)[[1]]
  eachPkg <- gsub(",", "", eachPkg, fixed = TRUE)
  length(eachPkg)
}


## -----------------------------------------------------------------------------
data(clin.data)
data(lesion.data)
data(expr.data)

head(lesion.data)

# Specify a folder on your local machine to store the analysis results:
# resultsPath=tempdir()
# knitr::opts_knit$set(root.dir = normalizePath(path = resultsPath))


## ----message=FALSE, warning=FALSE---------------------------------------------
hg19.ann=get.ensembl.annotation("Human_GRCh37") 
# "Human_GRCh38" can be used instead of "Human_GRCh37" to retrieve data for hg38

# 1) Gene annotation data that include around 20,000 coding genes and 25,000 Non-coding processed transcripts such as lncRNAs, miRNAs, snRNA and snoRNAs:
gene.annotation=hg19.ann$gene.annotation

# 2)Annotation data for regulatory features retrieved from ensembl regulatory build that  include around 500,000 feauters (promoters, enhancer, TF and CTCF binding sites, etc...). Ensembl imports publicly available data from different large epigenomic consortia that  includes ENCODE, Roadmap Epigenomics and Blueprint (118 epigenome):
hg19.reg.annotation=hg19.ann$reg.annotation.predicted

# 3)Annotation data for experimentally validated regulatory features retrieved from FANTOM5  project:
hg19.reg.FANTOM=hg19.ann$reg.annotation.validated

# Instead of retrieving annotation data from Ensembl BioMart, users can use their own gene  annotation data files. File should has four required columns that include "gene" which is the ensembl ID of annotated genes to which the lesion data will be overlapped, "chrom" which is the chromosome on which the gene is located, "loc.start" which is the gene start position, and "loc.end" the gene end position. hg19.gene annotation will be used as an example gene annotation data file:
data(hg19.gene.annotation)

head(hg19.gene.annotation)


## ----message=FALSE, warning=FALSE---------------------------------------------
# To retrieve chromosome size data for GRCh37 (hg19) genome build from chr.info txt file  available on UCSC genome browser
hg19.chrom.size=get.chrom.length("Human_GRCh37")
# "Human_GRCh38" can be used to retrieve chrom size data for hg38

# Instead of retrieving chromosome size data from UCSC genome browser, users can use their own files that should has two required columns that include "chrom" with the chromosome number and "size" for the size of the chromosome in base pairs:
# data(hg19.chrom.size)

head(hg19.chrom.size)


## ----message=FALSE, warning=FALSE---------------------------------------------
# Users can run GRIN analysis by just specifying the genome.version in grin.stats function. 
# A) Gene annotation data will be directly retrieved from Ensembl BioMart for the specified  genome assembly using get.ensembl.annotation function and chromosome size data will be also retrieved from UCSC genome browser:
# grin.results=grin.stats(lesion.data, 
#                         genome.version="Human_GRCh37")
# "Human_GRCh38" can be used instead of "Human_GRCh37" for hg38 genome assembly

# Users can also use their own annotation and chromosome size data files to run GRIN analysis:
grin.results=grin.stats(lesion.data, 
                        hg19.gene.annotation, 
                        hg19.chrom.size)
# it takes around 2 minutes to map 6,887 lesions to around 57,000 annotated genes and return the GRIN results.

# B) To run GRIN for computationally predicted regulatory features from Ensembl regulatory build:
# First get a group of 500 regulatory features for an example run: 
hg19.reg.example=hg19.reg.annotation[396500:397000,]
# whole file with around 500,000 feature takes around 25 minutes to return the results:
# Run GRIN analysis:
grin.results.reg=grin.stats(lesion.data, 
                            hg19.reg.example, 
                            hg19.chrom.size)

# C) To run GRIN analysis for experimentally verified regulatory features from FANTOM5 project:
# First get a group of 500 FANTOM5 regulatory features for an example run:
hg19.fantom.example=hg19.reg.FANTOM[232500:233000,]

grin.results.fantom=grin.stats(lesion.data, 
                               hg19.fantom.example, 
                               hg19.chrom.size)


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
                                    max.log10q=150) 
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


## ----oncorprint.constellation.test.significant.genes, fig.height =6, fig.width = 13, fig.cap="Figure 3. OncoPrint for a selected group of genes significant in the constelaation test for the gene to be affected by at least three types of lesions (q3.nsubj<0.05)"----
# Use onco.print.props function to specify a hgt for each lesion category to show all lesions that might affect a certain patient. For example, if the same patient is affected by gain and mutation, only 25% of the oncoprint rectangle will be filled with the mutations green color and the rest will appear as the gain red color.
onco.props<-onco.print.props(lesion.data,
                             hgt = c("gain"=5, "loss"=4, "mutation"=2, "fusion"=1))
column_title = "" # optional

# use oncoprint function from ComplexHeatmap library to plot the oncoprint:
oncoPrint(oncoprint.mtx,
          alter_fun = onco.props$alter_func,
          col = onco.props$col,
          column_title = column_title,
          heatmap_legend_param = onco.props$heatmap_legend_param)


## -----------------------------------------------------------------------------
# First we should call the pathways data file:
data(pathways)
head(pathways)

# define a list of pathways of interest:
PI3K_Pathway=pathways[pathways$pathway=="PI3K_Pathway",]
PI3K_ensembl=as.vector(PI3K_Pathway$ensembl.id)
Bcell_Pathway=pathways[pathways$pathway=="Bcell_Pathway",]
Bcell_ensembl=as.vector(Bcell_Pathway$ensembl.id)
Jak_Pathway=pathways[pathways$pathway=="Jak_Pathway",]
Jak_ensembl=as.vector(Jak_Pathway$ensembl.id)
Ras_Pathway=pathways[pathways$pathway=="Ras_Pathway",]
Ras_ensembl=as.vector(Ras_Pathway$ensembl.id)

oncoprint.genes=c(PI3K_ensembl, Bcell_ensembl, Jak_ensembl, Ras_ensembl)

# prepare the oncoprint matrix:
oncoprint.mtx.path=grin.oncoprint.mtx(grin.results,
                                      oncoprint.genes)
Gene=as.data.frame(rownames(oncoprint.mtx.path))
colnames(Gene)="gene.name"
Gene$index=1:nrow(Gene)
merged.df=merge(Gene,pathways, by="gene.name", all.x=TRUE)
merged.df=merged.df[order(merged.df$index), ]

sel.pathways=factor(merged.df$pathway,
                    levels=c("PI3K_Pathway", "Jak_Pathway", "Ras_Pathway", "Bcell_Pathway"))


## ----oncorprint.pathways, fig.height =9, fig.width = 13, fig.cap="Figure 4. OncoPrint for genes annotated to a selected group of pathways"----
# Use onco.print.props function to specify a hgt for each lesion category to show all lesions that might affect a certain patient:
onco.props.path<-onco.print.props(lesion.data,
                                  hgt = c("gain"=5, "loss"=4, "mutation"=2, "fusion"=1))
column_title = "" # optional

# use oncoprint function from complexheatmap library to plot the oncoprint
oncoPrint(oncoprint.mtx.path,
          alter_fun = onco.props.path$alter_func,
          col = onco.props.path$col,
          column_title = column_title,
          heatmap_legend_param = onco.props.path$heatmap_legend_param,
          row_split=sel.pathways)


## ----gene.plot.wt1, fig.keep='last', fig.height =8, fig.width = 9, fig.cap="Figure 5. lesion plot showing all different types of lesions affecting WT1 gene with transcripts track directly retreived from Ensembl database"----
# First we need to call "hg19_cytoband" and "hg38_cytoband" before calling the plot function:
data(hg19_cytoband)
data(hg38_cytoband)

# lsn.transcripts.plot function can be used to generate a plot that shows all different types of lesions that affect a gene of interest with a transcripts track directly retrieved from Ensembl genome browser:
lsn.transcripts.plot(grin.results, 
                     genome="hg19",  
                     gene="WT1", 
                     hg19.cytoband=hg19_cytoband)

# for hg38 genome assembly:
# library(AnnotationHub)
# ah <- AnnotationHub()
# retrieve gene transcripts for human GRCh38 genome assembly from Ensembl (version 110):
# gtf.V110 <- ah[["AH113665"]]

#lsn.transcripts.plot(grin.results, 
#                     genome="hg38",  
#                     gene="WT1", 
#                     hg38.transcripts=gtf.V110, 
#                     hg38.cytoband=hg38_cytoband)

## ----regional.lesion.plot, fig.keep='last', fig.height =10, fig.width = 8, fig.cap="Figure 6. Regional lesion plot showing a specific type of lesion that affect a region of interest"----
# lsn.transcripts.plot function can be also used to generate a plot for lesions of a specific lesion group that span a certain locus of interest with transcripts track added:
lsn.transcripts.plot(grin.results, 
                     genome="hg19", 
                     hg19.cytoband=hg19_cytoband,
                     chrom=9,
                     plot.start=21800000,
                     plot.end=22200000,
                     lesion.grp = "loss",
                     spec.lsn.clr = "blue")

# for hg38 genome assembly:
# ah <- AnnotationHub()
# retrieve gene transcripts for human GRCh38 genome assembly from Ensembl (version 110):
# gtf.V110 <- ah[["AH113665"]]

#lsn.transcripts.plot(grin.results,
#                     genome="hg38",
#                     hg38.transcripts="gtf.v110", 
#                     hg38.cytoband=hg38_cytoband,
#                     chrom=9,
#                     plot.start=21800000,
#                     plot.end=22200000,
#                     lesion.grp = "loss",
#                     spec.lsn.clr = "blue")


## ----locus.plot, fig.keep='last', fig.height =7, fig.width = 8, fig.cap="Figure 7. Regional lesion plot showing a specific type of lesion that affect a region of interest"----
# lsn.transcripts.plot function can be used to generate a plot for all lesions of a specific lesion type that affect a locus or region of interest without adding transcripts track. This will allow plotting a larger locus of the chromosome such as a chromosome band.transTrack argument should be set as FALSE.

lsn.transcripts.plot(grin.results, 
                     genome="hg19",
                     transTrack = FALSE,
                     hg19.cytoband=hg19_cytoband,
                     chrom=9,
                     plot.start=19900000,
                     plot.end=25600000,
                     lesion.grp = "loss",
                     spec.lsn.clr = "blue")

# for hg38 genome assembly:
# ah <- AnnotationHub()
# retrieve gene transcripts for human GRCh38 genome assembly from Ensembl (version 110):
# gtf.V110 <- ah[["AH113665"]]

#lsn.transcripts.plot(grin.results,
#               genome="hg38",
#               transTrack = FALSE,
#               hg38.transcripts="gtf.v110", 
#               hg38.cytoband=hg38_cytoband,
#               chrom=9,
#               plot.start=19900000,
#               plot.end=25600000,
#               lesion.grp = "loss",
#               spec.lsn.clr = "blue")

## ----chromosome.plot, fig.keep='last', fig.height =9, fig.width = 8, fig.cap="Figure 8. Lesion plot showing different types of lesions that affect a chromosome of interest"----
# lsn.transcripts.plot function can be also used to generate a plot for all types of lesions that affect a chromosome of interest with plot.start=1 and plot.end is the chr size:
lsn.transcripts.plot(grin.results,
                     genome="hg19",
                     transTrack = FALSE,
                     hg19.cytoband=hg19_cytoband,
                     chrom=9,
                     plot.start=1,
                     plot.end=141000000)

# for hg38 genome assembly:
#ah <- AnnotationHub()
#gtf.V110 <- ah[["AH113665"]]

#lsn.transcripts.plot(grin.results,
#                     genome="hg38",
#                     transTrack = FALSE,
#                     hg38.transcripts="gtf.v110", 
#                     hg38.cytoband=hg38_cytoband,
#                     chrom=9,
#                     plot.start=1, 
#                     plot.end=141000000)

## ----reg.feature.plot, fig.keep='last', fig.height =7, fig.width = 8, fig.cap="Figure 9. A plot that shows all different types of lesions that affect a regulatory feature of interests in addition the feature GRIN statistics"----
# grin.stats.lsn.plot function can be used to generate plots that show all different types of lesions that affect a regulatory feature of interest in addition to the GRIN statistics. Plot does not include transcripts track that's typically not available for those features. 
# grin.stats.lsn.plot(grin.results.reg, 
#                     feature="ENSR00000105619")

# Same plot can be also prepared for regulatory features from the FANTOM5 project (for example: the NRAS promoter site affected by 18 mutations)
grin.stats.lsn.plot(grin.results.fantom,
                    feature="p6@NRAS,0.2452")


## ----message=FALSE, warning=FALSE---------------------------------------------
# Prepare gene and lesion data for later computations
# This lesion matrix has all lesion types that affect a single gene in one row. It can be used to run association analysis with expression data (part of alex.prep.lsn.expr function)

# First step is to prepare gene and lesion data for later computations
gene.lsn=prep.gene.lsn.data(lesion.data,
                            hg19.gene.annotation)
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

alex.data=alex.prep.lsn.expr(expr.data,
                             lesion.data,
                             hg19.gene.annotation,
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
                               hg19.gene.annotation,
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

## ----message=FALSE, warning=FALSE---------------------------------------------
# return boxplots for a list of top significant genes to the pre-specified results folder:
# alex.boxplots(out.dir=resultsPath,
#               alex.data, alex.kw.results,
#                1e-15, hg19.gene.annotation)


## ----JAK2.waterfall.plot, fig.height =6, fig.width = 6, fig.cap="Figure 10. JAK2 Water-fall plot which  offers a side-by-side graphical representation of lesion and expression data for each patient"----
# waterfall plots allow a side-by-side representation of expression and lesion data of the gene of interest.
# First prepare expression and lesion data for waterfall plots:
WT1.waterfall.prep=alex.waterfall.prep(alex.data,
                                        alex.kw.results,
                                        "WT1",
                                        lesion.data)

# alex.waterfall.plot can be used to return the plot
WT1.waterfall.plot=alex.waterfall.plot(WT1.waterfall.prep,
                                        lesion.data)

## ----message=FALSE, warning=FALSE---------------------------------------------
# To prepare Waterfall plots for top significant genes in the KW Results Table, users can use top.alex.waterfall.plots function by specifying a directory to store the plots, and minimum KW.q, for example:

# top.alex.waterfall.plots(out.dir=resultsPath, 
#                          alex.data,
#                          alex.kw.results,
#                          1e-15,
#                          lesion.data)

## ----lesion.expression.pathway, fig.height =6, fig.width = 6, fig.cap="Figure 11. Ordered Lesion and Expression Data based on the Clustering Analysis on the pathway level (JAK/STAT pathway)"----
# alex.pathway function will run association analysis between lesion and expression data for all genes in a specified pathway (example: JAK/STAT pathway).
# Function will return two panels figure of lesion and expression data of ordered subjects based on the computed lesions distance in all genes assigned to the pathway of interest:
alex.path=alex.pathway(alex.data,
                       lesion.data,
                       pathways, 
                       "Jak_Pathway")

# To return ordered lesion and expression data of the genes assigned to the pathway of interest (same patients order in the plot):
alex.path[1:10,1:5]


## ----message=FALSE, warning=FALSE---------------------------------------------
# This type of lesion matrices with each gene affected by a certain type of lesion in a separate row is very helpful to run multiple levels of association analysis that include association between lesions and treatment outcomes.

# Users should first Prepare gene and lesion data and determine lesions that overlap each gene (locus):
gene.lsn=prep.gene.lsn.data(lesion.data,
                            hg19.gene.annotation)    
gene.lsn.overlap= find.gene.lsn.overlaps(gene.lsn)

# use prep.binary.lsn.mtx function to prepare the lesion binary matrix:
lsn.binary.mtx.atleast5=prep.binary.lsn.mtx(gene.lsn.overlap,
                                            min.ngrp=5)
# Each row is a lesion type that affect a certain gene for example NOTCH1_mutation (entry will be labelled as 1 if the patient is affected by by this type of lesion and 0 otherwise).
# min.ngrp can be used to specify the minimum number of patients with a lesion to be included in the final lesion matrix.

head(lsn.binary.mtx.atleast5[,1:5])


## ----message=FALSE, warning=FALSE---------------------------------------------
# Prepare Event-free Survival (EFS) and Overall Survival (OS) as survival objects:
clin.data$EFS <- Surv(clin.data$efs.time, clin.data$efs.censor)
clin.data$OS <- Surv(clin.data$os.time, clin.data$os.censor)

# List all clinical variables of interest to be included in the association analysis:
clinvars=c("MRD.binary", "EFS", "OS")

# Run association analysis between lesions and clinical variables:
assc.outcomes=grin.assoc.lsn.outcome(lsn.binary.mtx.atleast5,
                                     clin.data,
                                     hg19.gene.annotation,
                                     clinvars)


# Run models adjusted for one or a group of covariates:
# assc.outcomes.adj=grin.assoc.lsn.outcome(lsn.binary.mtx.atleast5,
#                                           clin.data,
#                                           hg19.gene.annotation,
#                                           clinvars,
#                                           covariate="Sex")



## -----------------------------------------------------------------------------
# order the genes by the ones with most significant KW q-value:
sorted.outcomes <- assc.outcomes[order(as.numeric(as.character(assc.outcomes$`MRD.binary.p-value`))),]

## -----------------------------------------------------------------------------
head(sorted.outcomes[1:7,c(6,11,14,15)])

## -----------------------------------------------------------------------------
head(sorted.outcomes[1:7,c(6, 16:19)])

## ----message=FALSE, warning=FALSE---------------------------------------------
# This analysis is lesion type specific and covers the entire genome.It's meant to cover and asses the regions without any annotated genes or regulatory features. The first boundary for each chromosome will start from the first nucleotide base on the chromosome till the start position of the first lesion that affect the chromosome. Similarly, the last boundary will start from the end position of the last lesion that affect the chromosome till the last base on the chromosome.

# First extract data for gains and deletions from the lesion data file:
gain=lesion.data[lesion.data$lsn.type=="gain",]
loss=lesion.data[lesion.data$lsn.type=="loss",]

# Then use grin.lsn.boundaries function to return the lesion boundaries:
lsn.bound.gain=grin.lsn.boundaries(gain, hg19.chrom.size)
lsn.bound.loss=grin.lsn.boundaries(loss, hg19.chrom.size)

# It return a table of ordered boundaries based on the unique start and end positions of different lesions in a specific category on each chromosome.
head(lsn.bound.loss[,1:5])


## ----message=FALSE, warning=FALSE---------------------------------------------
grin.results.gain.bound=grin.stats(gain,
                                   lsn.bound.gain,
                                   hg19.chrom.size)

grin.results.loss.bound=grin.stats(loss,
                                   lsn.bound.loss,
                                   hg19.chrom.size)


## ----loss.lesion.boundaries.significance.plot, fig.height = 8, fig.width = 8, fig.cap="Figure 12. Genome-wide -log10q plot of loss lesion boundaries"----
# genomewide.log10q.plot function will return a genome-wide plot based on -log(10) q-value testing if each of the evaluated lesion boundaries is significantly affect by a deletions in our example:
genomewide.log10q.plot(grin.results.loss.bound,
                       lsn.grps=c("loss"),
                       lsn.colors=c("loss" = "blue"),
                       max.log10q = 50)


## ----loss.annotated.genes.significance.plot, fig.height = 8, fig.width = 8, fig.cap="Figure 13. Genome-wide -log10q plot for annotated genes affected by deletions"----
# genomewide.log10q.plot function can be also used to return genome-wide significance plot for annotated genes to be affected by a certain type of lesions.
# Here we should use GRIN results for annotated genes affected by loss instead of lesion boundaries. Users can notice that some regions mostly without annotated markers were only captured in the lesion boundaries analysis that cover the entire genome:
genomewide.log10q.plot(grin.results, 
                       lsn.grps=c("loss"),
                       lsn.colors=c("loss" = "blue"),
                       max.log10q = 50)

## ----setup--------------------------------------------------------------------
library(GRIN2)

