
# p.order: Function compute ordered p-values for each gene based on the probability that the gene is affected by each type of lesions individually
# and return a constellation p-values for the gene to be affected by at least one type of lesions (p1), at least two types of lesions (p2), etc...

p.order=function(P) # matrix of lesion hit probabilities, rows for genes, columns for lesion types
{
  k=ncol(P)
  p.mtx=apply(P,1,sort,na.last=T)
  p.mtx=t(p.mtx)

  res=p.mtx
  for (i in 1:k)
  {
    res[,i]=stats::pbeta(p.mtx[,i],i,k-i+1)
  }
  return(res)
}

