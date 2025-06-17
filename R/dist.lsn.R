
# Computes pairwise distances between subjects based on lesion profiles.

dist.lsn=function(lsn.mtx) # subjects in columns, genes in rows
{
  n=ncol(lsn.mtx)
  res=matrix(NA,n,n)
  colnames(res)=rownames(res)=colnames(lsn.mtx)
  diag(res)=0
  for (i in 1:(n-1))
    for (j in (i+1):n)
    {
      res[i,j]=res[j,i]=sum(lsn.mtx[,i]!=lsn.mtx[,j])
    }
  return(res)
}
