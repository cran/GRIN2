
# compute the probability that a series of independent Bernoulli trials
# yields x or more successes.

rpbc=function(x,  # x: scalar number of Bernoulli successes observed
              p)  # p: vector of Bernoulli success probabilities

{
  n=length(p)
  if (x>n) return(0)
  if (x==n) return(exp(sum(log(p))))
  if (x<=0) return(1)
  res=pbc(p,x)
  return(res[x+1])
}
