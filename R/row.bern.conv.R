
# row.bern.conv: This function Compute a convolution of Bernoullis for each row of
# a Bernoulli success probability matrix (part of the prob.hits function).

row.bern.conv=function(P,          # matrix of lesion hit probabilities, rows for genes, columns for lesion types
                       max.x=NULL) # Maximum number of subjects or maximum number of hits

{
  m=nrow(P)
  n=ncol(P)
  if (is.null(max.x))
    max.x=(ncol(P))
  Pr=matrix(0,m,max.x+1)
  Pr[,1]=1

  for (i in 1:n)
  {
    P1=Pr*P[,i]
    P0=Pr*(1-P[,i])
    Pr0=P0
    Pr1=cbind(0,P1)
    Pr1[,max.x+1]=Pr1[,max.x+1]+Pr1[,max.x+2]
    Pr=Pr0+Pr1[,-(max.x+2)]
  }
  rs.Pr=rowSums(Pr)
  Pr=Pr/rs.Pr
  return(Pr)
}
