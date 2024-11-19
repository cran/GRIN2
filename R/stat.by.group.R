
# stat.by.group: Function compute summary statistics such as expression mean, median and standard deviation for each lesion group based on the specified test.

stat.by.group=function(x,              # vector of expression data
                       g,              # vector of group labels corresponding to x (lesion groups)
                       stat,           # stats to be computed for the expression data by lesion groups (mean, median or standard deviation)
                       all.grps=NULL,  # vector of all the possible group labels
                       ...)            # additional arguments to stat
  
{
  if (is.null(all.grps))               # if all.grps isn't specified, define it from g
    all.grps=sort(unique(g))
  
  all.grps=sort(unique(all.grps))      # ensure grps don't appear twice and are in the same order for each application
  
  res=rep(NA,length(all.grps))         # initialize the result vector
  names(res)=all.grps                  # name the elements of the result vector
  
  for (i in all.grps)                  # loop over the groups
  {
    in.grp=(g%in%i)
    if (any(in.grp))
      res[i]=stat(x[in.grp],...)
  }
  
  return(res)
  
}
