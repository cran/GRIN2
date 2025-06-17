
# row.stats.by.group: Function uses stat.by.group function to compute stats for each row of the
# expression and lesion data matrices and return a data.frame of expression level by lesion groups for
# each specified test.

row.stats.by.group=function(X,   # expression data
                            G,   # lesion data
                            stat, # stats to be computed
                            ...)

{
  if (any(dim(X)!=dim(G)))
    stop("X and G are of incompatible dimensions.  X and G must have the same dimensions.")

  if (any(colnames(X)!=colnames(G)))
    stop("X and G must have column names matched.")

  all.grps=sort(unique(G))
  all.grps=sort(unique(all.grps))

  k=length(all.grps)

  m=nrow(X)

  res=matrix(NA,m,k)
  colnames(res)=all.grps
  rownames(res)=rownames(X)

  for (i in 1:m)
  {
    res[i,]=stat.by.group(X[i,],G[i,],stat,all.grps,...)
  }

  return(res)
}
