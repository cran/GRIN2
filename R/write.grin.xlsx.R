
#' Write GRIN Results
#'
#' @description
#' The function Write GRIN results to an excel file with multiple sheets that include GRIN results, lesion data, gene annotation data, chromosome size, gene-lesion overlap and methods paragraph.
#'
#' @param grin.result output results of the grin.stats function.
#' @param output.file output file name ".xlsx".
#'
#' @return
#' This function return an excel file with seven sheets that include:
#' \item{gene.hits}{data table of GRIN results that include gene annotation, number of subjects affected by each lesion type for example gain, loss, mutation, etc.., and number of hits affecting each locus. The GRIN results table will also include P and FDR adjusted q-values showing the probability of each locus of being affected by one or a constellation of multiple types of lesions.}
#' \item{gene.lsn.data}{each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene ensembl ID and "ID"" column has the patient ID.}
#' \item{lsn.data}{input lesion data}
#' \item{gene.data}{input gene annotation data}
#' \item{chr.size}{data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs.}
#' \item{interpretation}{provides some details about  the content of each sheet in the output excel file and interpretation of each column in the "gene.hits" GRIN results table.}
#' \item{method.paragraph}{include a paragraph that explains the GRIN model and cite some references.}
#'
#' @export
#'
#' @importFrom writexl write_xlsx
#'
#' @references
#' Pounds, Stan, et al. (2013) A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [grin.stats()]
#'
#' @examples
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(hg19.chrom.size)
#'
#' # to directly retreive gene annotation and chromosome size files from Ensembl BioMart database,
#' # UCSC genome browsers and run the GRIN analysis:
#' grin.results=grin.stats(lesion.data,
#'                          hg19.gene.annotation,
#'                          hg19.chrom.size)
#'
#' # Write GRIN results in to an excel sheet ".xlsx" using write.grin.xlsx function.

write.grin.xlsx=function(grin.result, # output of the grin.stats function
                         output.file) # output file name ".xlsx"
{
  old_opt <- options(stringsAsFactors = FALSE)
  on.exit(options(old_opt), add = TRUE)

  rpt.res=grin.result[c("gene.hits",
                        "gene.lsn.data",
                        "lsn.data",
                        "gene.data",
                        "chr.size")]

  ################################
  # Contents of the data sheets
  sheet.int=cbind.data.frame(sheet.name=c("gene.hits",
                                          "gene.lsn.data",
                                          "lsn.data",
                                          "gene.data",
                                          "chr.size"),
                             col.name="entire sheet",
                             meaning=c("GRIN statistical results",
                                       "gene-lesion overlaps",
                                       "input lesion data",
                                       "input gene location data",
                                       "input chromosome size data"))

  ######################################
  # Interpretation of gene hit stats
  gh.cols=colnames(rpt.res[["gene.hits"]])
  genehit.int=cbind.data.frame(sheet.name="gene.hits",
                               col.name=gh.cols,
                               meaning=gh.cols)
  rownames(genehit.int)=gh.cols
  rpt.clms=c("gene.row","gene","loc.start","loc.end")
  rpt.defs=c("gene data row index",
             "gene name",
             "locus of left edge of gene",
             "locus of right edge of gene")

  rpt.indx=rep(NA,length(rpt.clms))
  for (i in 1:length(rpt.clms))
    rpt.indx[i]=which(genehit.int$col.name==rpt.clms[i])

  genehit.int$meaning[rpt.indx]=rpt.defs

  nsubj.clms=which(substring(gh.cols,1,5)=="nsubj")
  genehit.int$meaning[nsubj.clms]=paste0("Number of subjects with a ",
                                         substring(gh.cols[nsubj.clms],7)," lesion ",
                                         genehit.int$meaning[nsubj.clms])

  p.nsubj.clms=which(substring(gh.cols,1,7)=="p.nsubj")
  genehit.int$meaning[p.nsubj.clms]=paste0("p-value for the number of subjects with a ",
                                           substring(gh.cols[p.nsubj.clms],9)," lesion ",
                                           genehit.int$meaning[p.nsubj.clms])

  q.nsubj.clms=which(substring(gh.cols,1,7)=="q.nsubj")
  genehit.int$meaning[q.nsubj.clms]=paste0("FDR estimate for the number of subjects with a ",
                                           substring(gh.cols[q.nsubj.clms],9)," lesion ",
                                           genehit.int$meaning[q.nsubj.clms])

  nhit.clms=which(substring(gh.cols,1,4)=="nhit")
  genehit.int$meaning[nhit.clms]=paste0("Number of ",
                                        substring(gh.cols[nhit.clms],6)," lesions ",
                                        genehit.int$meaning[nhit.clms])

  p.nhit.clms=which(substring(gh.cols,1,6)=="p.nhit")
  genehit.int$meaning[p.nhit.clms]=paste0("p-value for the number of ",
                                          substring(gh.cols[p.nhit.clms],8)," lesions ",
                                          genehit.int$meaning[p.nhit.clms])

  q.nhit.clms=which(substring(gh.cols,1,6)=="q.nhit")
  genehit.int$meaning[q.nhit.clms]=paste0("FDR estimate for the number of ",
                                          substring(gh.cols[q.nhit.clms],8)," lesions ",
                                          genehit.int$meaning[q.nhit.clms])


  p1.nsubj.clms=paste0("p",1:length(nsubj.clms),".nsubj")
  genehit.int[p1.nsubj.clms,"meaning"]=paste0("p-value for the number of subjects ",
                                              "with any ",1:length(nsubj.clms)," type(s) of lesion ",
                                              "overlapping the gene locus")

  q1.nsubj.clms=paste0("q",1:length(nsubj.clms),".nsubj")
  genehit.int[q1.nsubj.clms,"meaning"]=paste0("FDR estimate for the number of subjects ",
                                              "with any ",1:length(nsubj.clms)," type(s) of lesion ",
                                              "overlapping the gene locus")

  p1.nhit.clms=paste0("p",1:length(nhit.clms),".nhit")
  genehit.int[p1.nhit.clms,"meaning"]=paste0("p-value for the number of ",
                                             "any ",1:length(nhit.clms)," type(s) of lesion ",
                                             "overlapping the gene locus")

  q1.nhit.clms=paste0("q",1:length(nhit.clms),".nhit")
  genehit.int[q1.nhit.clms,"meaning"]=paste0("FDR estimate for the number of",
                                             " any ",1:length(nhit.clms)," type(s) of lesion ",
                                             "overlapping the gene locus")

  #######################################
  # Lesion data interpretation
  lsn.clms=colnames(rpt.res[["lsn.data"]])
  lsn.int=cbind.data.frame(sheet.name="lsn.data",
                           col.name=lsn.clms,
                           meaning=lsn.clms)

  rpt.clms=c("ID","chrom","loc.start","loc.end","lsn.type")
  int.clms=c(which(lsn.int[,"col.name"]=="ID"),
             which(lsn.int[,"col.name"]=="chrom"),
             which(lsn.int[,"col.name"]=="loc.start"),
             which(lsn.int[,"col.name"]=="loc.end"),
             which(lsn.int[,"col.name"]=="lsn.type"))
  lsn.int[int.clms,"meaning"]=c("Input Subject Identifier",
                                "Input Chromosome",
                                "Input Gene Locus Left Edge",
                                "Input Gene Locus Right Edge",
                                "Input Lesion Type")

  ############################################
  # Gene data interpretation
  gene.clms=colnames(rpt.res[["gene.data"]])
  gene.int=cbind.data.frame(sheet.name="gene.data",
                            col.name=gene.clms,
                            meaning="Echoed from Input")
  rpt.clms=c("gene","chrom","loc.start","loc.end",
             "glp.row.start","glp.row.end")
  int.clms=c(which(gene.int[,"col.name"]=="gene"),
             which(gene.int[,"col.name"]=="chrom"),
             which(gene.int[,"col.name"]=="loc.start"),
             which(gene.int[,"col.name"]=="loc.end"))
  gene.int[int.clms,"meaning"]=c("Input Gene Locus Name",
                                 "Input Gene Locus Chromosome",
                                 "Input Gene Locus Left Edge",
                                 "Input Gene Locus Right Edge")

  ###############################################
  # Chromosome size data interpretation
  chr.clms=colnames(rpt.res[["chr.size"]])
  chr.int=cbind.data.frame(sheet.name="chr.size",
                           col.name=chr.clms,
                           meaning="Echoed from Input")
  int.clms=c(which(chr.int[,"col.name"]=="chrom"),
             which(chr.int[,"col.name"]=="size"))
  chr.int[int.clms,"meaning"]=c("Input Chromosome",
                                "Input Chromosome Size")


  rpt.res$interpretation=rbind.data.frame(sheet.int,
                                          genehit.int,
                                          lsn.int,
                                          gene.int,
                                          chr.int)


  #########################################
  # Write the methods paragraph
  lsn.types=sort(unique(rpt.res[["lsn.data"]][,"lsn.type"]))
  mp=c("The genomic random interval model [ref 1] was used to evaluate ",
       "the statistical significance of the number of subjects with ",
       paste0("each type of lesion (",
              paste0(lsn.types,collapse=","),")"),
       "in each gene.  For each type of lesion, robust false discovery ",
       "estimates were computed from p-values using Storey's q-value [ref 2] ",
       "with the Pounds-Cheng estimator of the proportion of hypothesis ",
       "tests with a true null [ref 3].  Additionally, p-values for the ",
       paste0("number of subjects with any 1 to ",length(lsn.types)," types of lesions "),
       "were computed using the beta distribution derived for order statistics ",
       "of the uniform distribution [ref 4].",
       "",
       "REFERENCES",
       "[ref 1] Pounds S, et al (2013)  A genomic random interval model for statistical analysis of genomic lesion data.  Bioinformatics, 2013 Sep 1;29(17):2088-95 (PMID: 23842812).",
       "[ref 2] Storey J (2002).  A direct approach to false discovery rates.  Journal of the Royal Statistical Society Series B.  64(3): 479-498.  (doi.org/10.1111/1467-9868.00346).",
       "[ref 3] Pounds S and Cheng C (2005)  Robust estimation of the false discovery rate.  Bioinformatics 22(16): 1979-87.  (PMID: 16777905).  ",
       "[ref 4] Casella G and Berger RL (1990)  Statistical Inference.  Wadsworth & Brooks/Cole: Pacific Grove, California.  Example 5.5.1.")

  rpt.res$methods.paragraph=cbind.data.frame(methods.paragraph=mp)

  m=length(rpt.res)
  for (i in 1:m)
  {
    rpt.res[[i]]=as.data.frame(rpt.res[[i]])
  }

  writexl::write_xlsx(rpt.res,output.file)

  return(invisible())

}
