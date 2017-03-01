#' QQ plot of all sites from DMS Analysis
#'
#' @param DMS output dataframe produced by runlimma \code{DMS}
#' @param methdiff difference in methylation \code{methdiff}
#' @param pval pvalue cutoff \code{pval}
#' 
#' @return QQplot single site analysis output
#' 
#' @keywords limma, QQplot 
#'
#' @export
#' 
#' @examples
#' Circos(DMS,methdiff=10,pval=0.05,
#' hyper.col="springgreen3",hypo.col="mediumvioletred")




Circos <- function(DMS, methdiff = 5,  pval = 0.01, 
hyper.col = "springgreen3", hypo.col = "mediumvioletred") {
	require(BSgenome)			#v1.30.0
	require("BSgenome.Hsapiens.UCSC.hg19")
	chr.len = seqlengths(Hsapiens)  # get chromosome lengths
	#remove M and random chromosomes
	chr.len = chr.len[grep("_|X|Y|M", names(chr.len), invert = T)]
	chrom.length=chr.len
    require(methylKit)
    require(GenomicRanges)		#v1.14.4
    require(ggbio)				#v1.10.14

	   # chrom.length
    myIdeo <- GRanges(seqnames = names(chrom.length), 
		ranges = IRanges(start = 1, width = chrom.length))
    seqlevels(myIdeo) = names(chrom.length)
    seqlengths(myIdeo) = (chrom.length)


    hyper = subset(DMS, DMS$mean.meth.diff>0 & DMS$mean.meth.diff>methdiff)
	hyper = subset(hyper, hyper$P.Value<=pval)
	
	hypo =  subset(DMS, DMS$mean.meth.diff<0& DMS$mean.meth.diff< -methdiff)
	hypo = subset(hypo, hypo$P.Value<=pval)
	
    g.per = as(hyper, "GRanges")
	seqlevels(g.per, force=TRUE) <- seqlevels(myIdeo)
	seqlengths(g.per)=(chrom.length)
	
	
    g.po = as(hypo, "GRanges")
 	seqlevels(g.po, force=TRUE) <- seqlevels(myIdeo)
	seqlengths(g.po)=(chrom.length)
	
    values(g.po)$id = "hypo"
    values(g.per)$id = "hyper"
       p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", 
            radius = 42, trackWidth = 2.6)
	#trackWidth adjusts thickness of grey bars
	#radius adjusts size of circle

        p <- p + layout_circle(c(g.po, g.per), geom = "point", 
                 size = 3.0, aes(x = midpoint, 
				 #this size adjusts the size of the points plotted
            y = mean.meth.diff, color = id), radius = 35, trackWidth = 28) +              
            scale_colour_manual(values = c(hyper.col, hypo.col))
        p + layout_circle(myIdeo, geom = "text", aes(label = seqnames), 
            vjust = 0, radius =55, trackWidth = 2,size=7) +theme(legend.text = element_text(colour="black", size = 20, face = "bold"))+ theme(legend.title=element_blank())
	#adjust this radius to move chrom labels further away (bigger number)
#dev.off()
}
