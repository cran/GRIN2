
#' Find Probability of Locus Hit
#'
#' @description
#' Computes the probability that each genomic locus (e.g., gene or regulatory region) is affected by one or more types of genomic lesions. This function estimates statistical significance for lesion enrichment using a convolution of independent but non-identical Bernoulli distributions.
#'
#' @param hit.cnt A list returned by the `count.hits()` function, containing the number of subjects and hits affecting each locus by lesion type.
#' @param chr.size A `data.frame` containing chromosome sizes for all 22 autosomes and the X and Y chromosomes. It must include two columns: \code{"chrom"} for chromosome number, and \code{"size"} for chromosome lengths in base pairs.
#'
#' @details
#'
#' This function estimates a p-value for each locus based on the probability of observing the observed number of lesions (or more) by chance, under a model where lesion events are treated as independent Bernoulli trials.
#'
#' For each lesion type, the model considers heterogeneity in lesion probability across loci based on their genomic context (e.g., locus size, chromosome size). These probabilities are then combined using a convolution of Bernoulli distributions to estimate the likelihood of observing the actual hit counts.
#'
#' In addition, the function calculates:
#' - **FDR-adjusted q-values** using the method of Pounds and Cheng (2006), which estimates the proportion of true null hypotheses.
#' - **p- and q-values for multi-lesion constellation hits**, i.e., the probability that a locus is affected by one (\code{p1}), two (\code{p2}), or more types of lesions simultaneously.
#'
#' @return
#' A list with the following components:
#' \item{gene.hits}{A `data.frame` containing GRIN statistical results. Includes gene annotations, the number of subjects and hits by lesion type, and the computed p-values and FDR-adjusted q-values for lesion enrichment across one or more lesion types.}
#' \item{lsn.data}{Original input lesion data.}
#' \item{gene.data}{Original input gene annotation data.}
#' \item{gene.lsn.data}{A `data.frame` in which each row corresponds to a gene overlapped by a specific lesion. Includes columns for Ensembl gene ID (\code{gene}) and patient/sample ID (\code{ID}).}
#' \item{chr.size}{Chromosome size information used in the computation.}
#' \item{gene.index}{A `data.frame` indexing rows in `gene.lsn.data` corresponding to each chromosome.}
#' \item{lsn.index}{A `data.frame` indexing rows in `gene.lsn.data` corresponding to each lesion.}
#'
#' @export
#'
#' @importFrom stats pbeta p.adjust
#'
#' @references
#' Pounds, S. et al. (2013). A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author
#' Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}

#'
#' @seealso \code{\link{prep.gene.lsn.data}}, \code{\link{find.gene.lsn.overlaps}}, \code{\link{count.hits}}
#'
#' @examples
#' data(lesion_data)
#' data(hg38_gene_annotation)
#' data(hg38_chrom_size)
#'
#' # 1) Prepare gene and lesion data:
#' prep.gene.lsn <- prep.gene.lsn.data(lesion_data, hg38_gene_annotation)
#'
#' # 2) Identify overlapping gene-lesion events:
#' gene.lsn.overlap <- find.gene.lsn.overlaps(prep.gene.lsn)
#'
#' # 3) Count number of subjects and lesions affecting each gene:
#' count.subj.hits <- count.hits(gene.lsn.overlap)
#'
#' # 4) Compute p- and q-values for lesion enrichment per gene:
#' hits.prob <- prob.hits(count.subj.hits, hg38_chrom_size)

prob.hits=function(hit.cnt, # Output results of the count.hits function
                   chr.size=NULL) # A data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs. data.frame should has two columns "chrom" with the chromosome number and "size" for the size of the chromosome in base pairs
{
  lsn.data=hit.cnt$lsn.data
  num.lsn=unique(lsn.data$lsn.type)
  ###################
  # order and index gene.lsn.data
  ord=order(hit.cnt$gene.lsn.data$lsn.type,
            hit.cnt$gene.lsn.data$lsn.chrom,
            hit.cnt$gene.lsn.data$gene.row,
            hit.cnt$gene.lsn.data$ID)
  hit.cnt$gene.lsn.data=hit.cnt$gene.lsn.data[ord,]
  m=nrow(hit.cnt$gene.lsn.data)

  new.sect=which((hit.cnt$gene.lsn.data$gene.chrom[-1]!=hit.cnt$gene.lsn.data$gene.chrom[-m])|
                   (hit.cnt$gene.lsn.data$lsn.type[-1]!=hit.cnt$gene.lsn.data$lsn.type[-m])|
                   (hit.cnt$gene.lsn.data$gene.row[-1]!=hit.cnt$gene.lsn.data$gene.row[-m]))
  sect.start=c(1,new.sect+1)
  sect.end=c(new.sect,m)
  gene.lsn.index=cbind.data.frame(lsn.type=hit.cnt$gene.lsn.data$lsn.type[sect.start],
                                  chrom=hit.cnt$gene.lsn.data$gene.chrom[sect.start],
                                  gene.row=hit.cnt$gene.lsn.data$gene.row[sect.start],
                                  row.start=sect.start,
                                  row.end=sect.end,
                                  n.lsns=sect.end-sect.start+1)
  k=nrow(gene.lsn.index)

  new.chr=which(gene.lsn.index$chrom[-1]!=gene.lsn.index$chrom[-k])
  chr.start=c(1,new.chr+1)
  chr.end=c(new.chr,k)
  gene.lsn.chr.index=cbind.data.frame(lsn.type=gene.lsn.index$lsn.type[chr.start],
                                      chrom=gene.lsn.index$chrom[chr.start],
                                      row.start=chr.start,
                                      row.end=chr.end,
                                      n.rows=chr.end-chr.start+1)

  nr.li=nrow(hit.cnt$lsn.index)
  new.chr=which((hit.cnt$lsn.index$lsn.type[-1]!=hit.cnt$lsn.index$lsn.type[-nr.li])|
                  (hit.cnt$lsn.index$chrom[-1]!=hit.cnt$lsn.index$chrom[-nr.li]))
  chr.start=c(1,new.chr+1)
  chr.end=c(new.chr,nr.li)
  lsn.chr.index=cbind.data.frame(lsn.type=hit.cnt$lsn.index$lsn.type[chr.start],
                                 chrom=hit.cnt$lsn.index$chrom[chr.start],
                                 row.start=chr.start,
                                 row.end=chr.end)

  b=nrow(gene.lsn.chr.index)
  g=nrow(hit.cnt$nhit.mtx)
  nlt=ncol(hit.cnt$nhit.mtx)

  p.nsubj=p.nhit=matrix(1,g,nlt)
  colnames(p.nsubj)=colnames(p.nhit)=colnames(hit.cnt$nhit.mtx)

  for (i in 1:b)
  {
    # find rows for affected genes
    gli.start.row=gene.lsn.chr.index$row.start[i]
    gli.end.row=gene.lsn.chr.index$row.end[i]
    gld.start.row=gene.lsn.index$row.start[gli.start.row]
    gld.end.row=gene.lsn.index$row.end[gli.end.row]
    gld.rows=gld.start.row:gld.end.row
    gene.rows=unique(hit.cnt$gene.lsn.data$gene.row[gld.rows])
    n.genes=length(gene.rows)

    # find rows for lesions of this type on this chromosomes
    lsn.chr.mtch=which((lsn.chr.index$lsn.type==gene.lsn.chr.index$lsn.type[i])&
                         (lsn.chr.index$chrom==gene.lsn.chr.index$chrom[i]))
    if (length(lsn.chr.mtch)>0)
    {
      lsn.index.start.row=lsn.chr.index$row.start[lsn.chr.mtch]
      lsn.index.end.row=lsn.chr.index$row.end[lsn.chr.mtch]
      lsn.start.row=hit.cnt$lsn.index$row.start[lsn.index.start.row]
      lsn.end.row=hit.cnt$lsn.index$row.end[lsn.index.end.row]
      lsn.rows=lsn.start.row:lsn.end.row
      n.lsns=length(lsn.rows)
      lsn.type=hit.cnt$lsn.data$lsn.type[lsn.start.row]

      message(paste0("Computing p-values for ",
                     n.genes," gene(s) on chromosome ",
                     gene.lsn.chr.index$chrom[i],
                     " affected by ",
                     n.lsns," ",lsn.type,
                     " (data block ",i," of ",b,"): ",date()))

      # find chromosome size
      chr.mtch=which(gene.lsn.chr.index$chrom[i]==chr.size$chrom)
      chrom.size=chr.size$size[chr.mtch]

      # obtain gene sizes, lesion sizes, and gene hit probabilities
      lsn.size=hit.cnt$lsn.data$loc.end[lsn.rows]-hit.cnt$lsn.data$loc.start[lsn.rows]+1
      gene.size=hit.cnt$gene.data$loc.end[gene.rows]-hit.cnt$gene.data$loc.start[gene.rows]+1

      log.pr=log(rep(lsn.size,each=n.genes)+rep(gene.size,times=n.lsns))-log(chrom.size)
      pr.gene.hit=matrix(exp(log.pr),n.genes,n.lsns)
      pr.gene.hit[pr.gene.hit>1]=1

      lsn.subj.IDs=hit.cnt$lsn.data$ID[lsn.rows]
      pr.subj=row.prob.subj.hit(pr.gene.hit,lsn.subj.IDs)

      #max.nsubj=max(hit.cnt$nsubj.mtx[gene.rows,lsn.type])
      #max.nhit=max(hit.cnt$nhit.mtx[gene.rows,lsn.type])

      #pr.nhit=row.bern.conv(pr.gene.hit,max.nhit) # in original GRIN library
      #pr.nsubj=row.bern.conv(pr.subj,max.nsubj)   # in original GRIN library

      for (j in 1:n.genes)
      {
        nsubj=hit.cnt$nsubj.mtx[gene.rows[j],lsn.type]
        nhit=hit.cnt$nhit.mtx[gene.rows[j],lsn.type]
        p.nsubj[gene.rows[j],lsn.type]=rpbc(nsubj,pr.subj[j,])
        p.nhit[gene.rows[j],lsn.type]=rpbc(nhit,pr.gene.hit[j,])
      }
    }

  }

  rownames(p.nhit)=rownames(hit.cnt$nhit.mtx)
  rownames(p.nsubj)=rownames(hit.cnt$nsubj.mtx)
  colnames(hit.cnt$nhit.mtx)=paste0("nhit.",colnames(hit.cnt$nhit.mtx))
  colnames(hit.cnt$nsubj.mtx)=paste0("nsubj.",colnames(hit.cnt$nsubj.mtx))
  colnames(p.nhit)=paste0("p.",colnames(hit.cnt$nhit.mtx))
  colnames(p.nsubj)=paste0("p.",colnames(hit.cnt$nsubj.mtx))

  # Compute q-values
  message(paste0("Computing q-values: ",date()))
  q.nhit=p.nhit
  q.nsubj=p.nsubj
  for (i in 1:ncol(q.nhit))
  {
    pi.hat=min(1,2*mean(p.nhit[,i],na.rm=T))
    q.nhit[,i]=pi.hat*stats::p.adjust(p.nhit[,i],method="fdr")
    pi.hat=min(1,2*mean(p.nsubj[,i],na.rm=T))
    q.nsubj[,i]=pi.hat*stats::p.adjust(p.nsubj[,i],method="fdr")
  }

  colnames(q.nhit)=paste0("q.",colnames(hit.cnt$nhit.mtx))
  colnames(q.nsubj)=paste0("q.",colnames(hit.cnt$nsubj.mtx))

  gd.clms=setdiff(colnames(hit.cnt$gene.data),c("glp.row.start","glp.row.end"))
  lsn.clms=setdiff(colnames(hit.cnt$lsn.data),c("glp.row.start","glp.row.end"))

  gd.clms=c("gene.row",setdiff(gd.clms,"gene.row"))
  lsn.clms=c("lsn.row",setdiff(lsn.clms,"lsn.row"))

  if (length(num.lsn) >1) {
    # Now get ordered p-values
    message(paste0("Computing p-values for number of lesion types affecting genes: ",date()))
    p.ord.nhit=p.order(p.nhit)
    colnames(p.ord.nhit)=paste0("p",1:ncol(p.nhit),".nhit")

    p.ord.nsubj=p.order(p.nsubj)
    colnames(p.ord.nsubj)=paste0("p",1:ncol(p.nsubj),".nsubj")

    # q-values of ordered p-values
    q.ord.nhit=p.ord.nhit
    q.ord.nsubj=p.ord.nsubj
    message(paste0("Computing q-values for number of lesion types affecting genes: ",date()))
    for (i in 1:ncol(p.ord.nhit))
    {
      pi.hat=min(1,2*mean(p.ord.nhit[,i],na.rm=T))
      q.ord.nhit[,i]=pi.hat*stats::p.adjust(p.ord.nhit[,i],method="fdr")
      pi.hat=min(1,2*mean(p.ord.nsubj[,i],na.rm=T))
      q.ord.nsubj[,i]=pi.hat*stats::p.adjust(p.ord.nsubj[,i],method="fdr")
    }
    colnames(q.ord.nsubj)=paste0("q",1:ncol(p.nsubj),".nsubj")
    colnames(q.ord.nhit)=paste0("q",1:ncol(p.nhit),".nhit")

    gene.res=cbind.data.frame(hit.cnt$gene.data[,gd.clms],
                              hit.cnt$nsubj.mtx,
                              p.nsubj,
                              q.nsubj,
                              p.ord.nsubj,
                              q.ord.nsubj,
                              hit.cnt$nhit.mtx,
                              p.nhit,
                              q.nhit,
                              p.ord.nhit,
                              q.ord.nhit)
  }

  # if lesion data has only one type of lesion, skip the constellation test for locus to be affected by multiple types of lesions
  else if (length(num.lsn)==1) {

    gene.res=cbind.data.frame(hit.cnt$gene.data[,gd.clms],
                              hit.cnt$nsubj.mtx,
                              p.nsubj,
                              q.nsubj,
                              hit.cnt$nhit.mtx,
                              p.nhit,
                              q.nhit)
  }

  res=list(gene.hits=gene.res, # A data table of GRIN results that includes gene annotation, number of subjects and number of hits affecting each locus. in addition, p and FDR adjusted q-values showing the probability of each locus being affected by one or a constellation of multiple types of lesions are also included in the GRIN results.
           lsn.data=hit.cnt$lsn.data[,lsn.clms], # Input lesion data
           gene.data=hit.cnt$gene.data[,gd.clms], # Input gene data
           gene.lsn.data=hit.cnt$gene.lsn.data, # Each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene and ID column has the patient ID
           chr.size=chr.size, # A data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs
           gene.index=hit.cnt$gene.index, # Data.frame that shows row start and row end for each chromosome in the gene.lsn.data table
           lsn.index=hit.cnt$lsn.index) # Data.frame that shows row start and row end for each lesion in the gene.lsn.data table

  return(res)
}
