
# row.prob.subj.hit: Compute the probability that a subject has a hit for each gene (part of the prob.hits function).

row.prob.subj.hit=function(P,   # matrix of lesion hit probabilities, rows for genes, columns for lesion types
                           IDs) # vector of subject IDs for each lesion, length must equal ncol(P)
{
  if (length(IDs)!=ncol(P))
    stop("length(IDs) must equal ncol(P).")

  g=nrow(P)
  l=ncol(P)

  ord=order(IDs)
  IDs=IDs[ord]
  P=matrix(P[,ord],g,l)
  l=length(IDs)

  new.ID=which(IDs[-1]!=IDs[-l])
  ID.start=c(1,new.ID+1)
  ID.end=c(new.ID,l)

  n=length(ID.start)
  m=nrow(P)
  Pr=matrix(NA,m,n)
  for (i in 1:n)
  {
    pr.mtx=matrix(P[,ID.start[i]:ID.end[i]],m,ID.end[i]-ID.start[i]+1)
    Pr[,i]=1-exp(rowSums(log(1-pr.mtx)))
  }

  return(Pr)
}
