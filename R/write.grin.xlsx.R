
#' Write GRIN Results to Excel File
#'
#' @description
#' Writes GRIN results to an Excel file containing multiple sheets. The output includes the GRIN summary statistics, input data (lesion and gene annotation), chromosome sizes, gene-lesion overlaps, and explanatory metadata to help with the results interpretation.
#'
#' @param grin.result A list returned by the \code{\link{grin.stats}} function, containing GRIN analysis output.
#' @param output.file A character string specifying the name of the output Excel file. Must end with \code{".xlsx"}.
#'
#' @return
#' The function creates a multi-sheet Excel file at the specified location. The file contains the following sheets:
#' \itemize{
#'   \item \code{gene.hits}: The GRIN results table. Includes gene annotation, number of subjects affected by each lesion type (e.g., gain, loss, mutation), total lesion hits per gene, and associated p-values and FDR-adjusted q-values for the probability of lesion enrichment.
#'   \item \code{gene.lsn.data}: A table where each row corresponds to a lesion overlapping a specific gene. Columns include \code{"gene"} (Ensembl gene ID) and \code{"ID"} (patient identifier).
#'   \item \code{lsn.data}: The input lesion dataset used in the GRIN analysis.
#'   \item \code{gene.data}: The input gene annotation dataset, typically retrieved from Ensembl.
#'   \item \code{chr.size}: A table listing the size (in base pairs) of chromosomes 1:22, X, and Y.
#'   \item \code{interpretation}: A guide to understanding the structure and content of each sheet, with detailed descriptions of columns in the \code{gene.hits} results table.
#'   \item \code{method.paragraph}: A summary of the GRIN methodology, including relevant references for citation.
#' }
#'
#' @export
#'
#' @importFrom writexl write_xlsx
#'
#' @references
#' Pounds, S., et al. (2013). A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}
#'
#' @seealso \code{\link{grin.stats}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#' data(hg38_chrom_size)
#'
#' # Run GRIN analysis using lesion, gene, and chromosome size data:
#' grin.results <- grin.stats(lesion_data,
#'                            hg38_gene_annotation,
#'                            hg38_chrom_size)
#'
#' # Write results to an Excel file:
#' tmp_file <- file.path(tempdir(), "GRIN_Results.xlsx")
#' write.grin.xlsx(grin.results, output.file = tmp_file)
#' if (file.exists(tmp_file)) file.remove(tmp_file)

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
