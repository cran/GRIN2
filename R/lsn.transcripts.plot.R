
#' Lesions Gene Transcripts Plot
#'
#' @description
#' Function prepare a plot with all types of lesions that spans either a gene or a region of interest.
#'
#' @param grin.res GRIN results (output of the grin.stats function).
#' @param genome either "hg19" or "hg38" genome assemblies can be specified based on the genome assembly that has been used to prepare the lesion data.
#' @param gene Gene name of interest.
#' @param transTrack In case of plots that span large genomic region such as a chromosome band or the whole chromosome, this argument should be specified as 'FALSE' to exclude the transcripts track from the plot.
#' @param lsn.clrs Lesion colors for the regional gene plot (If not provided by the user, colors will be automatically assigned using default.grin.colors function).
#' @param chrom chromosome number (should be only specified in the locus plots where plot.start and plot.end for the locus of interest are specified).
#' @param plot.start start position of the locus of interest.
#' @param plot.end end position of the locus of interest.
#' @param lesion.grp lesion group of interest (should be only specified in locus plots when chrom, plot.start, plot.end are specified).
#' @param spec.lsn.clr color assigned to the lesion of interest (should be specified when chrom, plot.start, plot.end and lesion.grp are specified).
#' @param extend.left specified number will be used to manually align the left side of the gene transcripts track directly retrieved from ensembl database with the gene lesions track if needed.
#' @param extend.right specified number will be used to manually align the right side of the gene transcripts track directly retrieved from ensembl database with the gene lesions track if needed.
#' @param expand Controls ratio of the gene locus (start and end position) to the whole plot with default value = 0.0005 (setting expand=0 will only plot the gene locus from the start to the end position without any of the upstream or downstream regions of the gene).
#' @param hg38.transcripts transcripts data retrieved from annotation hub for hg38 version 110 (should be only specified if genome="hg38").
#' @param hg19.cytoband hg19 chromosome bands start and end data in base pair (should be only specified if genome="hg19").
#' @param hg38.cytoband hg38 chromosome bands start and end data in base pair (should be only specified if genome="hg38").
#'
#' @details
#' Function return a plot with all lesions that affect either a gene or a region of interest. Top panel of the regional gene plot has the transcripts track with all transcripts annotated to the gene of interest directly retrieved from ensembl database. The middle panel will has all different types of lesions affecting the gene color coded according to the figure legend. Lower panel of the plot has all the GRIN statistics of the gene that include number of subjects affected by each type of lesions, -log10 p, and –log10q values showing if the gene is significantly affected by the corresponding lesion category. If a certain locus is specified, only transcripts track and the lesion panel will be returned (GRIN results panel will not be added to the plot). In case of plots that span large genomic region such as a chromosome band or the whole chromosome and by specifying transTrack=FALSE, transcripts track will not be added to the plot as well.
#'
#' @return
#' Function will return either a gene plot with the transcripts track, lesions panel and GRIN statistic for the gene of interest, a plot with all lesions and transcripts aligned to a certain locus of interest if chrom, plot.start and plot.end were specified or a plot with all lesions affecting a region of interest without the transcripts track.
#'
#' @export
#'
#' @importFrom graphics rect legend text segments
#' @importFrom gridGraphics grid.echo
#' @importFrom grid grid.grab editGrob grid.draw popViewport grid.newpage pushViewport viewport gpar grid.rect
#' @importFrom Gviz plotTracks IdeogramTrack GenomeAxisTrack GeneRegionTrack
#' @importFrom ensembldb getGeneRegionTrackForGviz
#' @importFrom EnsDb.Hsapiens.v75 EnsDb.Hsapiens.v75
#' @importFrom GenomeInfoDb seqlevelsStyle
#'
#' @references
#' Cao, X., Elsayed, A. H., & Pounds, S. B. (2023). Statistical Methods Inspired by Challenges in Pediatric Cancer Multi-omics.
#'
#' @author {Abdelrahman Elsayed \email{abdelrahman.elsayed@stjude.org} and Stanley Pounds \email{stanley.pounds@stjude.org}}
#'
#' @seealso [grin.stats()]
#'
#' @examples
#' \donttest{
#' data(lesion.data)
#' data(hg19.gene.annotation)
#' data(hg19.chrom.size)
#' data(hg19_cytoband)
#' data(hg38_cytoband)
#'
#' # run GRIN analysis using grin.stats function
#' grin.results=grin.stats(lesion.data,
#'                         hg19.gene.annotation,
#'                         hg19.chrom.size)
#'
#' # Plots showing different types of lesions affecting a gene of interest with a transcripts
#' # track that show all the gene transcripts retrieved from Ensembl (hg19 genome assembly):
#' WT1.gene.plot=lsn.transcripts.plot(grin.results, genome="hg19", gene="WT1",
#'                                     hg19.cytoband=hg19_cytoband)
#'
#' # Plots showing different types of lesions affecting a region of interest with a transcripts
#' # track added to the plot:
#' locus.plot=lsn.transcripts.plot(grin.results, genome="hg19", hg19.cytoband=hg19_cytoband,
#'                                 chrom=9, plot.start=21800000, plot.end=22200000,
#'                                 lesion.grp = "loss", spec.lsn.clr = "blue")
#'
#' # Plots Showing Different Types of Lesions Affecting a region of Interest without plotting the
#' # transcripts track (this will allow plotting a larger locus of the chromosome such as a
#' # chromosome band (should specify transTrack = FALSE):
#' noTranscripts=lsn.transcripts.plot(grin.results, genome="hg19", transTrack = FALSE,
#'                                    hg19.cytoband=hg19_cytoband, chrom=9, plot.start=19900000,
#'                                    plot.end=25600000, lesion.grp = "loss", spec.lsn.clr = "blue")
#'
#'  # Plots Showing Different Types of Lesions Affecting the whole chromosome:
#'  chrom.plot=lsn.transcripts.plot(grin.results, genome="hg19", transTrack = FALSE,
#'                                  hg19.cytoband=hg19_cytoband, chrom=9, plot.start=1,
#'                                   plot.end=141000000)
#'
#' # for GRCh38 (hg38) genome assembly, users should first call the AnnotationHub() web resource then
#' # specify ah[["AH113665"]] to retrieve the human hg38 gene transcripts. This formal class
#' # EnsDb object should be called afterwards in the 'hg38.transcripts' argument to return gene and
#' # regional plots.
#' }

lsn.transcripts.plot=function(grin.res,          # GRIN results (output of the grin.stats function)
                              genome,            # genome assembly (hg19 or hg38)
                              gene=NULL,         # gene name (should be only specified in the regional gene plots)
                              transTrack=TRUE,   # if specified as 'FALSE', transcripts track will not be added to the plot
                              lsn.clrs=NULL,     # Specified colors per lesion types (gene plots when gene name is specified). If not specified, colors will be automatically assigned using default.grin.colors function
                              chrom=NULL,        # chromosome number (should be only specified in the locus plots where plot.start and plot.end for the locus of interest are specified)
                              plot.start=NULL,   # start position of the locus of interest
                              plot.end=NULL,     # end position of the locus of interest
                              lesion.grp=NULL,   # lesion group of interest (should be only specified in locus plots when chrom, plot.start, plot.end are specified)
                              spec.lsn.clr=NULL, # color of the lesion of interest (should be specified when chrom, plot.start, plot.end and lesion.grp are specified)
                              extend.left=NULL,  # specified number will be used to manually align the left side of the gene transcripts track directly retrieved from ensembl database with the gene lesions track
                              extend.right=NULL, # specified number will be used to manually align the right side of the gene transcripts track directly retrieved from ensembl database with the gene lesions track
                              expand=0.0005,     # Controls ratio of the gene locus (start and end position) to the whole plot with default value = 0.0005 (setting expand=0 will only plot the gene locus from the start to the end position without any of the upstream or downstream regions of the gene)
                              hg38.transcripts=NULL, # transcripts data retrieved from annotation hub for hg38 version 110 (should be only specified if genome="hg38")
                              hg19.cytoband=NULL,    # hg19 chromosome bands start and end data retrieved from UCSC genome browser (should be only specified if genome="hg19")
                              hg38.cytoband=NULL)    # hg38 chromosome bands start and end data retrieved from UCSC genome browser (should be only specified if genome="hg38")

{
  # regional gene plot (gene name should be specified)
  if (transTrack==TRUE) {

    if (is.character(gene))
    {
      # Find the requested gene
      gene.data=grin.res[["gene.data"]]
      gene.mtch=which(gene.data[,"gene.name"]==gene)

      if (length(gene)!=1)
        stop("Exactly one gene must be specified!")
      if (length(gene.mtch)==0)
        stop(paste0(gene," not found in gene.data."))

      if (length(gene.mtch)>1)
        stop(paste0("Multiple matches of ",gene," found in gene data."))

      gene.chr=gene.data[gene.mtch,"chrom"]
      gene.start=gene.data[gene.mtch,"loc.start"]
      gene.end=gene.data[gene.mtch,"loc.end"]
      gene.size=(gene.end-gene.start)+1

      # Find lesions on the same chromosome as the gene
      lsn.dset=order.index.lsn.data(grin.res[["lsn.data"]])
      lsn.data=lsn.dset$lsn.data

      lsn.types=unique(lsn.data$lsn.type)
      if (is.null(lsn.clrs))
        lsn.clrs=default.grin.colors(lsn.types)

      lsn.index=lsn.dset$lsn.index
      lsn.ind.mtch=which(lsn.index$chrom==gene.chr)
      if (length(lsn.ind.mtch)==0)
        stop(paste0("No lesions overlap ",gene,"."))

      lsn.chr.rows=NULL
      for (i in lsn.ind.mtch)
      {
        blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
        lsn.chr.rows=c(lsn.chr.rows,blk.rows)
      }
      lsn.chr.rows=unlist(lsn.chr.rows)

      lsn.chr.data=lsn.data[lsn.chr.rows,]

      if (any(lsn.chr.data$chrom!=gene.chr))
        stop(paste0("Error in finding lesions on same chromosome as gene ",gene,"."))

      # Find lesions at overlap the gene
      ov.rows=which((lsn.chr.data$loc.start<=gene.end)&(lsn.chr.data$loc.end>=gene.start))
      if (length(ov.rows)==0)
        stop(paste0("No lesions overlap gene ",gene,"."))

      lsn.gene=lsn.chr.data[ov.rows,]

      lsn.gene$size=lsn.gene$loc.end-
        lsn.gene$loc.start+1
      lsn.gene=lsn.gene[order(lsn.gene$lsn.type, lsn.gene$size), ]

      # define plotting data
      x.start=gene.start-expand*gene.size
      x.end=gene.end+expand*gene.size

      lsn.gene$index=1:nrow(lsn.gene)
      lsn.gene$subj.num=as.numeric(as.factor(lsn.gene$index))
      lsn.gene$lsn.clr=lsn.clrs[lsn.gene$lsn.type]
      lsn.gene$type.num=as.numeric(as.factor(lsn.gene$lsn.type))


      n.type=max(lsn.gene$type.num)
      n.subj=max(lsn.gene$subj.num)

      lsn.gene$y0=-lsn.gene$subj.num+(lsn.gene$type.num-1)/n.type
      lsn.gene$y1=-lsn.gene$subj.num+lsn.gene$type.num/n.type

      lsn.gene$x0=pmax(lsn.gene$loc.start,x.start)
      lsn.gene$x1=pmin(lsn.gene$loc.end,x.end)

      # based on how many characters in the gene name, the right side of the lesion data plot should be shifted so that the two panels will be aligned
      if (nchar(gene)==3)
      {
        plot(c(x.start-0.055*(x.end-x.start),x.end+0.12*(x.end-x.start)),
             c(+0.0000001,-1.75)*n.subj,type="n",
             main="",
             xlab="",
             ylab="",axes=F)
      }
      else if (nchar(gene)==4) {
        plot(c(x.start-0.055*(x.end-x.start),x.end+0.13*(x.end-x.start)),
             c(+0.0000001,-1.75)*n.subj,type="n",
             main="",
             xlab="",
             ylab="",axes=F)
      }

      else if (nchar(gene)==5) {
        plot(c(x.start-0.055*(x.end-x.start),x.end+0.145*(x.end-x.start)),
             c(+0.0000001,-1.75)*n.subj,type="n",
             main="",
             xlab="",
             ylab="",axes=F)
      }

      else if (nchar(gene)==6) {
        plot(c(x.start-0.055*(x.end-x.start),x.end+0.1625*(x.end-x.start)),
             c(+0.0000001,-1.75)*n.subj,type="n",
             main="",
             xlab="",
             ylab="",axes=F)
      }
      else if (nchar(gene)==7) {
        plot(c(x.start-0.055*(x.end-x.start),x.end+0.175*(x.end-x.start)),
             c(+0.0000001,-1.75)*n.subj,type="n",
             main="",
             xlab="",
             ylab="",axes=F)
      }
      else if (nchar(gene)==8) {
        plot(c(x.start-0.055*(x.end-x.start),x.end+0.185*(x.end-x.start)),
             c(+0.0000001,-1.75)*n.subj,type="n",
             main="",
             xlab="",
             ylab="",axes=F)
      }
      else if (nchar(gene)>8) {
        plot(c(x.start-0.055*(x.end-x.start),x.end+0.21*(x.end-x.start)),
             c(+0.0000001,-1.75)*n.subj,type="n",
             main="",
             xlab="",
             ylab="",axes=F)
      }

      graphics::rect(x.start,
                     -(1:n.subj),
                     x.end,
                     -(1:n.subj)+1,
                     col=c("snow","gainsboro")[1+(1:n.subj)%%2],
                     border=NA)

      graphics::segments(c(gene.start,gene.end),
                         rep(-n.subj,2),
                         c(gene.start,gene.end),
                         rep(0,2),
                         col="darkgray")

      graphics::rect(lsn.gene$x0,
                     lsn.gene$y0,
                     lsn.gene$x1,
                     lsn.gene$y1,
                     col=lsn.gene$lsn.clr,
                     border=lsn.gene$lsn.clr)

      graphics::segments(c(gene.start,gene.end),
                         rep(-n.subj,2),
                         c(gene.start,gene.end),
                         rep(0,2),
                         col="darkgray",lty=2)

      graphics::text(c(gene.start,gene.end),
                     -n.subj,
                     c(gene.start,gene.end),
                     pos=1, cex=0.75)


      lgd=graphics::legend((x.start+x.end)/2,-1.10*n.subj,
                           fill=lsn.clrs,
                           legend=names(lsn.clrs),
                           ncol=length(lsn.clrs),
                           xjust=0.4,border=NA,
                           cex=0.7,bty="n")

      graphics::text(lgd$text$x[1]-0.08*diff(range(lgd$text$x)),
                     -c(1.25,1.35,1.45)*n.subj,
                     c("n","-log10p","-log10q"),pos=2, cex=0.7, font=2)

      gene.stats=grin.res[["gene.hits"]]
      stat.mtch=which(gene.stats$gene.name==gene)
      gene.stats=gene.stats[stat.mtch,]

      graphics::text(lgd$text$x,-1.25*n.subj,
                     gene.stats[,paste0("nsubj.",names(lsn.clrs))],
                     cex=0.7)
      graphics::text(lgd$text$x,-1.35*n.subj,
                     round(-log10(gene.stats[,paste0("p.nsubj.",names(lsn.clrs))]),2),
                     cex=0.7)
      graphics::text(lgd$text$x,-1.45*n.subj,
                     round(-log10(gene.stats[,paste0("q.nsubj.",names(lsn.clrs))]),2),
                     cex=0.7)

      lgd2=graphics::legend((x.start+x.end)/2,-1.5*n.subj,
                            legend=paste0("const", 1:length(lsn.types), ".typ"),
                            ncol=length(lsn.types),
                            xjust=0.42,border=NA,
                            cex=0.7,bty="n")

      graphics::text(lgd2$text$x[1]-0.0005*diff(range(lgd2$text$x)),
                     -c(1.65,1.75)*n.subj,
                     c("const.p", "const.q"),pos=2, cex=0.7, font=2)

      graphics::text(lgd2$text$x,-1.65*n.subj,
                     round(-log10(gene.stats[,paste0("p",1:length(lsn.types), ".nsubj")]),2),
                     cex=0.7, pos=4)
      graphics::text(lgd2$text$x,-1.75*n.subj,
                     round(-log10(gene.stats[,paste0("q",1:length(lsn.types), ".nsubj")]),2),
                     cex=0.7, pos=4)


      gridGraphics::grid.echo()
      lesion.plt <- grid::grid.grab()
      lesion.plt <- grid::editGrob(lesion.plt, gp=grid::gpar(fontsize=12))

      ## add the track of gene transcripts directly retrived from ensembl database based on the genome assembly
      if (genome=="hg19")
      {
        edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75      # same version of hg19 genome assembly used to retrieve annotation data in the get.ensembl.annotation function to make sure that coordinates will be consistent between gene annotation and the transcripts track files
        GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
        txs <- ensembldb::getGeneRegionTrackForGviz(edb, filter = ~ genename == gene)

        ## Define the individual tracks:
        ## - Ideogram
        hg19_cytoband=hg19.cytoband
        ideo_track <- Gviz::IdeogramTrack(genome = "hg19", chromosome = gene.chr, bands = hg19_cytoband)
        ## - Genome axis
        gaxis_track <- Gviz::GenomeAxisTrack()
        ## - Transcripts
        gene_track <- Gviz::GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                                            name = "Transcripts", geneSymbol = TRUE, size = 0.5)

        grid::grid.newpage()

        # specify plotting regions for gene lesions and transcript tracks
        grid::pushViewport(grid::viewport(height=0.78, width=1.2, y=0, just="bottom"))
        grid::grid.draw(lesion.plt)
        grid::popViewport(1)

        grid::pushViewport(grid::viewport(height=0.33, width=0.915, y=1, just="top"))
        Gviz::plotTracks(list(ideo_track, gaxis_track, gene_track), add=TRUE, sizes = c(1,2,6))
        grid::popViewport(1)

      }

      if (genome=="hg38")
      {
        trans.hg38=hg38.transcripts
        edb <- trans.hg38    # same version of hg38 genome assembly used to retrieve annotation data in the get.ensembl.annotation function to make sure that coordinates will be consistent between gene annotation and the transcripts track files
        GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
        txs <- ensembldb::getGeneRegionTrackForGviz(edb, filter = ~ genename == gene)

        ## Define the individual tracks:
        ## - Ideogram
        hg38_cytoband=hg38.cytoband
        ideo_track <- Gviz::IdeogramTrack(genome = "hg38", chromosome = gene.chr, bands = hg38_cytoband)
        ## - Genome axis
        gaxis_track <- Gviz::GenomeAxisTrack()
        ## - Transcripts
        gene_track <- Gviz::GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                                            name = "Transcripts", geneSymbol = TRUE, size = 0.5)

        grid::grid.newpage()

        # specify plotting regions for gene lesions and transcript tracks
        grid::pushViewport(grid::viewport(height=0.78, width=1.2, y=0, just="bottom"))
        grid::grid.draw(lesion.plt)
        grid::popViewport(1)

        grid::pushViewport(grid::viewport(height=0.33, width=0.915, y=1, just="top"))
        Gviz::plotTracks(list(ideo_track, gaxis_track, gene_track), add=TRUE, sizes = c(1,2,6))
        grid::popViewport(1)

      }
    }
    else {

      # if we want to specify a locus of interest instead of a gene
      locus.chr=chrom
      locus.start=plot.start
      locus.end=plot.end
      locus.size=(locus.end-locus.start)+1

      # Find lesions on the specified region of the chromosome
      lesion.data=grin.res[["lsn.data"]]
      lsn.type=lesion.grp
      lsn.data=lesion.data[lesion.data$lsn.type==lsn.type,]
      lsn.data=lsn.data[lsn.data$chrom==locus.chr,]
      lsn.dset=order.index.lsn.data(lsn.data)
      lsn.clr=spec.lsn.clr

      lsn.index=lsn.dset$lsn.index
      lsn.ind.mtch=which(lsn.index$chrom==locus.chr)
      if (length(lsn.ind.mtch)==0)
        stop(paste0("No lesions overlap"))

      lsn.chr.rows=NULL
      for (i in lsn.ind.mtch)
      {
        blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
        lsn.chr.rows=c(lsn.chr.rows,blk.rows)
      }
      lsn.chr.rows=unlist(lsn.chr.rows)

      lsn.chr.data=lsn.data[lsn.chr.rows,]
      lsn.chr.data=lsn.chr.data[lsn.chr.data$chrom==locus.chr,]

      if (any(lsn.chr.data$chrom!=locus.chr))
        stop(paste0("Error in finding lesions on same chromosome"))

      # Find lesions at overlap the locus
      ov.rows=which((lsn.chr.data$loc.start<=locus.end)&(lsn.chr.data$loc.end>=locus.start))
      if (length(ov.rows)==0)
        stop(paste0("No lesions overlap locus"))

      lsn.locus=lsn.chr.data[ov.rows,]

      lsn.locus$size=lsn.locus$loc.end-lsn.locus$loc.start+1
      lsn.locus=lsn.locus[order(lsn.locus$lsn.type, lsn.locus$size), ]

      # define plotting data
      x.start=locus.start-expand*locus.size
      x.end=locus.end+expand*locus.size

      lsn.locus$index=1:nrow(lsn.locus)
      lsn.locus$subj.num=as.numeric(as.factor(lsn.locus$index))
      lsn.locus$lsn.clr=lsn.clr
      lsn.locus$type.num=as.numeric(as.factor(lsn.locus$lsn.type))
      n.type=max(lsn.locus$type.num)
      n.subj=max(lsn.locus$subj.num)

      lsn.locus$y0=-lsn.locus$subj.num+(lsn.locus$type.num-1)/n.type
      lsn.locus$y1=-lsn.locus$subj.num+lsn.locus$type.num/n.type

      lsn.locus$x0=pmax(lsn.locus$loc.start,x.start)
      lsn.locus$x1=pmin(lsn.locus$loc.end,x.end)

      plot(c(x.start-0.05*(x.end-x.start),x.end+0.08*(x.end-x.start)),
           c(+0.1,-1.1)*n.subj,type="n",
           main="",
           xlab="",
           ylab="",axes=F)

      graphics::rect(x.start,
                     -(1:n.subj),
                     x.end,
                     -(1:n.subj)+1,
                     col=c("snow","gainsboro")[1+(1:n.subj)%%2],
                     border=NA)

      graphics::segments(c(locus.start,locus.end),
                         rep(-n.subj,2),
                         c(locus.start,locus.end),
                         rep(0,2),
                         col="darkgray")

      graphics::rect(lsn.locus$x0,
                     lsn.locus$y0,
                     lsn.locus$x1,
                     lsn.locus$y1,
                     col=lsn.locus$lsn.clr,
                     border=lsn.locus$lsn.clr)

      graphics::segments(c(locus.start,locus.end),
                         rep(-n.subj,2),
                         c(locus.start,locus.end),
                         rep(0,2),
                         col="darkgray",lty=2)

      graphics::text(c(locus.start,locus.end),
                     -n.subj,
                     c(locus.start,locus.end),
                     pos=1, cex=0.85)
      graphics::text((locus.start+locus.end)/2,
                     0,paste0(lsn.type, " _ ", "chr",locus.chr, ": ", locus.start, " - ", locus.end), pos=3,cex=1)

      gridGraphics::grid.echo()
      lesion.plt <- grid::grid.grab()
      lesion.plt <- grid::editGrob(lesion.plt, gp=grid::gpar(fontsize=12))

      ## add the transcripts track for all genes on the specified region  based on the genome assembly
      if (genome=="hg19")
      {
        edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
        GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
        txs <- ensembldb::getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                                    end=plot.end)

        ## Define the individual tracks:
        ## - Ideogram
        hg19_cytoband=hg19.cytoband
        ideo_track <- Gviz::IdeogramTrack(genome = "hg19", chromosome = locus.chr, bands = hg19_cytoband)
        ## - Genome axis
        gaxis_track <- Gviz::GenomeAxisTrack()
        ## - Transcripts
        gene_track <- Gviz::GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                                            name = "Transcripts", geneSymbol = TRUE, size = 0.5)

        grid::grid.newpage()

        # specify plotting regions for the locus lesions and transcript tracks
        grid::pushViewport(grid::viewport(height=0.45, width=1.2, y=0, just="bottom"))
        grid::grid.draw(lesion.plt)
        grid::popViewport(1)

        grid::pushViewport(grid::viewport(height=0.6, width=0.915, y=1, just="top"))
        Gviz::plotTracks(list(ideo_track, gaxis_track, gene_track), add=TRUE)
        grid::popViewport(1)

      }

      if (genome=="hg38")
      {
        trans.hg38=hg38.transcripts
        edb <- trans.hg38
        GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
        txs <- ensembldb::getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                                    end=plot.end)

        ## Define the individual tracks:
        ## - Ideogram
        hg38_cytoband=hg38.cytoband
        ideo_track <- Gviz::IdeogramTrack(genome = "hg38", chromosome = locus.chr, bands = hg38_cytoband)
        ## - Genome axis
        gaxis_track <- Gviz::GenomeAxisTrack()
        ## - Transcripts
        gene_track <- Gviz::GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                                            name = "Transcripts", geneSymbol = TRUE, size = 0.5)

        grid::grid.newpage()

        grid::pushViewport(grid::viewport(height=0.45, width=1.2, y=0, just="bottom"))
        grid::grid.draw(lesion.plt)
        grid::popViewport(1)

        grid::pushViewport(grid::viewport(height=0.6, width=0.915, y=1, just="top"))
        Gviz::plotTracks(list(ideo_track, gaxis_track, gene_track), add=TRUE)
        grid::popViewport(1)
      }
    }

  }

  else
  {
    if (is.character(lesion.grp)) {
      # To specify a locus of interest for the plotting purpose
      locus.chr=chrom
      locus.start=plot.start
      locus.end=plot.end
      locus.size=(locus.end-locus.start)+1

      # Find lesions on the specified region of the chromosome
      lesion.data=grin.res[["lsn.data"]]
      lsn.type=lesion.grp
      lsn.data=lesion.data[lesion.data$lsn.type==lsn.type,]
      lsn.data=lsn.data[lsn.data$chrom==locus.chr,]
      lsn.dset=order.index.lsn.data(lsn.data)
      lsn.clr=spec.lsn.clr

      lsn.index=lsn.dset$lsn.index
      lsn.ind.mtch=which(lsn.index$chrom==locus.chr)
      if (length(lsn.ind.mtch)==0)
        stop(paste0("No lesions overlap"))

      lsn.chr.rows=NULL
      for (i in lsn.ind.mtch)
      {
        blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
        lsn.chr.rows=c(lsn.chr.rows,blk.rows)
      }
      lsn.chr.rows=unlist(lsn.chr.rows)

      lsn.chr.data=lsn.data[lsn.chr.rows,]
      lsn.chr.data=lsn.chr.data[lsn.chr.data$chrom==locus.chr,]

      if (any(lsn.chr.data$chrom!=locus.chr))
        stop(paste0("Error in finding lesions on same chromosome"))

      # Find lesions that overlap the locus
      ov.rows=which((lsn.chr.data$loc.start<=locus.end)&(lsn.chr.data$loc.end>=locus.start))
      if (length(ov.rows)==0)
        stop(paste0("No lesions overlap locus"))

      lsn.locus=lsn.chr.data[ov.rows,]

      lsn.locus$size=lsn.locus$loc.end-lsn.locus$loc.start+1
      lsn.locus=lsn.locus[order(lsn.locus$lsn.type, lsn.locus$size), ]

      # define plotting data
      x.start=locus.start-expand*locus.size
      x.end=locus.end+expand*locus.size

      lsn.locus$index=1:nrow(lsn.locus)
      lsn.locus$subj.num=as.numeric(as.factor(lsn.locus$index))
      lsn.locus$lsn.clr=lsn.clr
      lsn.locus$type.num=as.numeric(as.factor(lsn.locus$lsn.type))
      n.type=max(lsn.locus$type.num)
      n.subj=max(lsn.locus$subj.num)

      lsn.locus$y0=-lsn.locus$subj.num+(lsn.locus$type.num-1)/n.type
      lsn.locus$y1=-lsn.locus$subj.num+lsn.locus$type.num/n.type

      lsn.locus$x0=pmax(lsn.locus$loc.start,x.start)
      lsn.locus$x1=pmin(lsn.locus$loc.end,x.end)

      plot(c(x.start-0.02*(x.end-x.start),x.end+0.1*(x.end-x.start)),
           c(+0.1,-1.1)*n.subj,type="n",
           main="",
           xlab="",
           ylab="",axes=F)

      graphics::rect(x.start,
                     -(1:n.subj),
                     x.end,
                     -(1:n.subj)+1,
                     col=c("snow","gainsboro")[1+(1:n.subj)%%2],
                     border=NA)

      graphics::segments(c(locus.start,locus.end),
                         rep(-n.subj,2),
                         c(locus.start,locus.end),
                         rep(0,2),
                         col="darkgray")

      graphics::rect(lsn.locus$x0,
                     lsn.locus$y0,
                     lsn.locus$x1,
                     lsn.locus$y1,
                     col=lsn.locus$lsn.clr,
                     border=lsn.locus$lsn.clr)

      graphics::segments(c(locus.start,locus.end),
                         rep(-n.subj,2),
                         c(locus.start,locus.end),
                         rep(0,2),
                         col="darkgray",lty=2)

      graphics::text(c(locus.start,locus.end),
                     -n.subj,
                     c(locus.start,locus.end),
                     pos=1, cex=0.85)
      graphics::text((locus.start+locus.end)/2,
                     0,paste0(lsn.type, " _ ", "chr",locus.chr, ": ", locus.start, " - ", locus.end), pos=3,cex=1)

      gridGraphics::grid.echo()
      lesion.plt <- grid::grid.grab()
      lesion.plt <- grid::editGrob(lesion.plt, gp=grid::gpar(fontsize=12))

      ## add the transcripts track for all genes on the specified region  based on the genome assembly
      if (genome=="hg19")
      {
        edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
        GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
        txs <- ensembldb::getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                                    end=plot.end)

        ## Define the individual tracks:
        ## - Ideogram
        hg19_cytoband=hg19.cytoband
        ideo_track <- Gviz::IdeogramTrack(genome = "hg19", chromosome = locus.chr,
                                          bands = hg19_cytoband, from=plot.start, to=plot.end)

        grid::grid.newpage()

        # specify plotting regions for the locus lesions and transcript tracks
        grid::pushViewport(grid::viewport(height=1, width=1.25, y=0, just="bottom"))
        grid::grid.draw(lesion.plt)
        grid::popViewport(1)

        grid::pushViewport(grid::viewport(height=0.15, width=0.88, y=1, just="top"))
        Gviz::plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                         showBandId = TRUE, cex.bands = 0.4)
        grid::popViewport(1)

      }

      if (genome=="hg38")
      {
        trans.hg38=hg38.transcripts
        edb <- trans.hg38
        GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
        txs <- ensembldb::getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                                    end=plot.end)

        ## Define the individual tracks:
        ## - Ideogram
        hg38_cytoband=hg38.cytoband
        ideo_track <- Gviz::IdeogramTrack(genome = "hg38", chromosome = locus.chr,
                                          bands = hg38_cytoband, from=plot.start, to=plot.end)

        grid::grid.newpage()

        grid::pushViewport(grid::viewport(height=1, width=1.25, y=0, just="bottom"))
        grid::grid.draw(lesion.plt)
        grid::popViewport(1)

        grid::pushViewport(grid::viewport(height=0.15, width=0.88, y=1, just="top"))
        Gviz::plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                         showBandId = TRUE, cex.bands = 0.4)
        grid::popViewport(1)
      }
    }
    else if (is.null(lesion.grp)) {
      # To specify a locus of interest for the plotting purpose
      locus.chr=chrom
      locus.start=plot.start
      locus.end=plot.end
      locus.size=(locus.end-locus.start)+1

      # Find lesions on the same chromosome as the gene
      lsn.dset=order.index.lsn.data(grin.res[["lsn.data"]])
      lsn.data=lsn.dset$lsn.data

      lsn.types=unique(lsn.data$lsn.type)
      if (is.null(lsn.clrs))
        lsn.clrs=default.grin.colors(lsn.types)

      lsn.index=lsn.dset$lsn.index
      lsn.ind.mtch=which(lsn.index$chrom==locus.chr)
      if (length(lsn.ind.mtch)==0)
        stop(paste0("No lesions overlap ",chrom,"."))

      lsn.chr.rows=NULL
      for (i in lsn.ind.mtch)
      {
        blk.rows=(lsn.index$row.start[i]:lsn.index$row.end[i])
        lsn.chr.rows=c(lsn.chr.rows,blk.rows)
      }
      lsn.chr.rows=unlist(lsn.chr.rows)

      lsn.chr.data=lsn.data[lsn.chr.rows,]

      if (any(lsn.chr.data$chrom!=locus.chr))
        stop(paste0("Error in finding lesions on same chromosome."))

      # Find lesions at overlap the locus
      ov.rows=which((lsn.chr.data$loc.start<=locus.end)&(lsn.chr.data$loc.end>=locus.start))
      if (length(ov.rows)==0)
        stop(paste0("No lesions overlap the chromosome."))

      lsn.locus=lsn.chr.data[ov.rows,]

      lsn.locus$size=lsn.locus$loc.end-lsn.locus$loc.start+1
      lsn.locus=lsn.locus[order(lsn.locus$lsn.type, lsn.locus$size), ]

      # define plotting data
      x.start=locus.start-expand*locus.size
      x.end=locus.end+expand*locus.size

      lsn.locus$index=1:nrow(lsn.locus)
      lsn.locus$subj.num=as.numeric(as.factor(lsn.locus$index))
      lsn.locus$lsn.clr=lsn.clrs[lsn.locus$lsn.type]
      lsn.locus$type.num=as.numeric(as.factor(lsn.locus$lsn.type))
      n.type=max(lsn.locus$type.num)
      n.subj=max(lsn.locus$subj.num)

      lsn.locus$y0=-lsn.locus$subj.num+(lsn.locus$type.num-1)/n.type
      lsn.locus$y1=-lsn.locus$subj.num+lsn.locus$type.num/n.type

      lsn.locus$x0=pmax(lsn.locus$loc.start,x.start)
      lsn.locus$x1=pmin(lsn.locus$loc.end,x.end)

      plot(c(x.start-0.02*(x.end-x.start),x.end+0.1*(x.end-x.start)),
           c(+0.1,-1.15)*n.subj,type="n",
           main="",
           xlab="",
           ylab="",axes=F)

      graphics::rect(x.start,
                     -(1:n.subj),
                     x.end,
                     -(1:n.subj)+1,
                     col=c("snow","gainsboro")[1+(1:n.subj)%%2],
                     border=NA)

      graphics::segments(c(locus.start,locus.end),
                         rep(-n.subj,2),
                         c(locus.start,locus.end),
                         rep(0,2),
                         col="darkgray")

      graphics::rect(lsn.locus$x0,
                     lsn.locus$y0,
                     lsn.locus$x1,
                     lsn.locus$y1,
                     col=lsn.locus$lsn.clr,
                     border=lsn.locus$lsn.clr)

      graphics::segments(c(locus.start,locus.end),
                         rep(-n.subj,2),
                         c(locus.start,locus.end),
                         rep(0,2),
                         col="darkgray",lty=2)

      graphics::text(c(locus.start,locus.end),
                     -n.subj,
                     c(locus.start,locus.end),
                     pos=1, cex=0.85)
      graphics::text((locus.start+locus.end)/2,
                     0,paste0("chr",locus.chr, ": ", locus.start, " - ", locus.end), pos=3,cex=1)

      lgd=graphics::legend((x.start+x.end)/2,-1.1*n.subj,
                           fill=lsn.clrs,
                           legend=names(lsn.clrs),
                           ncol=length(lsn.clrs),
                           xjust=0.45,border=NA,
                           cex=0.75,bty="n")

      gridGraphics::grid.echo()
      lesion.plt <- grid::grid.grab()
      lesion.plt <- grid::editGrob(lesion.plt, gp=grid::gpar(fontsize=12))

      ## add the transcripts track for all genes on the specified region  based on the genome assembly
      if (genome=="hg19")
      {
        edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
        GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
        txs <- ensembldb::getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                                    end=plot.end)

        ## Define the individual tracks:
        ## - Ideogram
        hg19_cytoband=hg19.cytoband
        ideo_track <- Gviz::IdeogramTrack(genome = "hg19", chromosome = locus.chr,
                                          bands = hg19_cytoband, from=plot.start, to=plot.end)

        grid::grid.newpage()

        # specify plotting regions for the locus lesions and transcript tracks
        grid::pushViewport(grid::viewport(height=1, width=1.25, y=0, just="bottom"))
        grid::grid.draw(lesion.plt)
        grid::popViewport(1)

        grid::pushViewport(grid::viewport(height=0.15, width=0.88, y=1, just="top"))
        Gviz::plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                         showBandId = TRUE, cex.bands = 0.4)
        grid::popViewport(1)

      }

      if (genome=="hg38")
      {
        trans.hg38=hg38.transcripts
        edb <- trans.hg38
        GenomeInfoDb::seqlevelsStyle(edb) <- "UCSC"
        txs <- ensembldb::getGeneRegionTrackForGviz(edb, chromosome=locus.chr, start=plot.start,
                                                    end=plot.end)

        ## Define the individual tracks:
        ## - Ideogram
        hg38_cytoband=hg38.cytoband
        ideo_track <- Gviz::IdeogramTrack(genome = "hg38", chromosome = locus.chr,
                                          bands = hg38_cytoband, from=plot.start, to=plot.end)

        grid::grid.newpage()

        grid::pushViewport(grid::viewport(height=1, width=1.25, y=0, just="bottom"))
        grid::grid.draw(lesion.plt)
        grid::popViewport(1)

        grid::pushViewport(grid::viewport(height=0.15, width=0.88, y=1, just="top"))
        Gviz::plotTracks(ideo_track, from = plot.start, to = plot.end, add=TRUE, showId = FALSE,
                         showBandId = TRUE, cex.bands = 0.4)
        grid::popViewport(1)
      }
    }
  }
}
