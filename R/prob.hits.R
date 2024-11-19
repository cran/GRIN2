
#' Find Probablity of Locus Hit
#'
#' @description
#' The function evaluates the probability of a locus to be affected by one or a constellation of multiple types of lesions.
#'
#' @param hit.cnt output results of the count.hits function with number of subjects and number of hits affecting each locus.
#' @param chr.size data.frame with the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs. The data.frame should has two columns "chrom" with the chromosome number and "size" for the size of the chromosome in base pairs.
#'
#' @details
#' The function computes p-value for the probability of each locus (gene or regulatory feature) to be affected by different types of lesions based on a convolution of independent but non-identical Bernoulli distributions to determine whether a certain locus has an abundance of lesions that is statistically significant.In addition, FDR-adjusted q value is computed for each locus based on Pounds & Cheng (2006) estimator of the proportion of tests with a true null (pi.hat). The function also evaluates if a certain locus is affected by a constellation of multiple types of lesions and computes a p and adjusted q values for the locus to be affected by one type of lesions (p1), two types of lesions (p2), etc...
#'
#' @return
#' A list with the following components:
#' \item{gene.hits}{data table of GRIN results that include gene annotation, number of subjects affected by each lesion type for example gain, loss, mutation, etc.., and number of hits affecting each locus. The GRIN results table will also include P and FDR adjusted q-values showing the probability of each locus of being affected by one or a constellation of multiple types of lesions.}
#' \item{lsn.data}{input lesion data}
#' \item{gene.data}{input gene annotation data}
#' \item{gene.lsn.data}{each row represent a gene overlapped by a certain lesion. Column "gene" shows the overlapped gene ensembl ID and "ID"" column has the patient ID.}
#' \item{chr.size}{data table showing the size of the 22 autosomes, in addition to X and Y chromosomes in base pairs.}
#' \item{gene.index}{data.frame with overlapped gene-lesion data rows that belong to each chromosome in the gene.lsn.data table.}
#' \item{lsn.index}{data.frame that shows the overlapped gene-lesion data rows taht belong to each lesion in the gene.lsn.data table.}
#'
#' @export
#'
#' @importFrom stats pbeta p.adjust
#'
#' @references
#' Pounds, Stan, et al. (2013) A genomic random interval model for statistical analysis of genomic lesion data.
#'
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [prep.gene.lsn.data()], [find.gene.lsn.overlaps()], [count.hits()]
#'
#' @examples
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(hg19.chrom.size)
#'
#' # prepare gene and lesion data for later computations:
#' prep.gene.lsn=prep.gene.lsn.data(lesion.data,
#'                                  hg19.gene.annotation)
#'
#' # determine lesions that overlap each gene (locus):
#' gene.lsn.overlap=find.gene.lsn.overlaps(prep.gene.lsn)
#'
#' # count number of subjects affected by different types of lesions and number of hits that affect
#' # each locus:
#' count.subj.hits=count.hits(gene.lsn.overlap)
#'
#' # compute the probability of each locus to be affected by one or a constellation of multiple
#' # types of lesion
#' hits.prob=prob.hits(count.subj.hits, hg19.chrom.size)

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

      max.nsubj=max(hit.cnt$nsubj.mtx[gene.rows,lsn.type])
      max.nhit=max(hit.cnt$nhit.mtx[gene.rows,lsn.type])

      pr.nhit=row.bern.conv(pr.gene.hit,max.nhit)
      pr.nsubj=row.bern.conv(pr.subj,max.nsubj)

      for (j in 1:n.genes)
      {
        nsubj=hit.cnt$nsubj.mtx[gene.rows[j],lsn.type]
        nhit=hit.cnt$nhit.mtx[gene.rows[j],lsn.type]
        p.nsubj[gene.rows[j],lsn.type]=sum(pr.nsubj[j,(nsubj+1):(max.nsubj+1)])
        p.nhit[gene.rows[j],lsn.type]=sum(pr.nhit[j,(nhit+1):(max.nhit+1)])
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
