
# one.KW.pvalue: Function Computes the KW p-value for one row of the lesion data matrix
# paired with one row of the expression data matrix

one.KW.pvalue=function(one.row.mtch,      # one row of the row.mtch matrix (output of the alex.prep.lsn.expr function)
                       expr.mtx,          # expression data matrix (alex.expr, output of the alex.prep.lsn.expr function)
                       hit.grps,          # lesion matrix (alex.lsn, output of the alex.prep.lsn.expr function)
                       min.grp.size=NULL) # minimum group size to perform the test (there should be at least two groups with number of patients > min.grp.size)

{
  # extract data for the test
  expr.row=one.row.mtch["expr.row"]
  hit.row=one.row.mtch["hit.row"]
  y=expr.mtx[expr.row,]
  grps=hit.grps[hit.row,]

  # identify data to include in the test
  grps=define.grps(grps,min.grp.size)
  exc=(grps=="EXCLUDE")
  grps=grps[!exc]
  y=y[!exc]

  # return NA if there is no test to perform
  uniq.grps=unique(grps)
  if (length(uniq.grps)<2) return(NA)

  # perform the KW test
  kw.res=stats::kruskal.test(y~grps)

  # return the KW p-value
  return(kw.res$p.value)
}

