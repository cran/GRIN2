
# Compute the probability that a series of independent Bernoulli trials
# yields up to x successes

pbc=function(p,      # p: vector of Bernoulli success probabilities
             x=NULL) # x: scalar number of Bernoulli successes observed

{
  n=length(p)
  if (is.null(x)) x=n
  if (x<=0) return(1)

  k=min(n,x)
  pr=rep(0,k+1); pr[1]=1
  for (i in 1:n)
  {
    j=min(i,k)
    pr0=pr[1:(j+1)]*(1-p[i])
    pr1=pr[1:(j+1)]*p[i]
    pr1[j]=pr1[j]+pr1[j+1]
    pr1[-1]=pr1[-(j+1)]
    pr1[1]=0
    pr[1:(j+1)]=pr0+pr1
  }
  return(pr)
}
