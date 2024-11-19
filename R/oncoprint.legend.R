
# oncoprint.legend: function specify name and colors of the oncoprint legend lesion groups

oncoprint.legend = function(lsn.type){
  res <- list()
  res$title <- 'Lesion Category'
  res$at <- lsn.type
  res$labels <- sub('^(\\w?)', '\\U\\1', lsn.type, perl=T)
  return(res)
}

