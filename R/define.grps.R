
# define.grps: Function defines the lesion groups to be included in the Kruskal-Wallis test.
# There should be at least two groups with the minimum number of subjects with lesions as specified in the min.grp.size, otherwise the test will return NA).
# The function exclude lesion groups that might have a very small number of patients that might exaggerates the KW test results.

define.grps=function(grps,               # Lesion groups
                     min.grp.size=NULL)  # minimum number of patients in the lesion group to be included in the KW test.
{
  grp.tbl=table(grps)
  small.grp=(grp.tbl<min.grp.size)
  exclude.grp=grps%in%(names(grp.tbl)[small.grp])
  grps[exclude.grp]="EXCLUDE"
  return(grps)
}
