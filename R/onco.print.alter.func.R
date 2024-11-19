
# onco.print.alter.func: The function define some important characteristics of the oncoprint such as
# background color, rectangles height and width.

onco.print.alter.func=function(type.wdth.hgt.clr,  # data.frame with columns for lesion type, box size, box color
                               bg.clr="#CCCCCC",   # background color
                               bg.hgt=2,           # background height
                               bg.wdth=0.00005)    # background width

{
  res = NULL
  size = nrow(type.wdth.hgt.clr)
  bg.code=paste0("background = function(x,y,w,h){",
                 "grid::grid.rect(x,y,w-unit(",bg.wdth,',"pt"),',
                 "h-unit(",bg.hgt,',"pt"),',
                 "gp = grid::gpar(fill = '",bg.clr,"',col=NA))}")
  names(bg.code) <- "background"
  Rcode=paste0(type.wdth.hgt.clr[,"type"]," = function(x,y,w,h){",
               "grid::grid.rect(x,y,w-unit(",type.wdth.hgt.clr[,"wdth"],',"pt"),',
               "h-unit(h*",(1-(type.wdth.hgt.clr[,"hgt"]*(1/(size+1)))),',"pt"),',
               "gp = grid::gpar(fill ='",type.wdth.hgt.clr[,"clr"],"',col=NA))}")

  names(Rcode) = type.wdth.hgt.clr[,"type"]
  Rcode=c(bg.code,Rcode)
  return(Rcode)
}
