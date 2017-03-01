#' Karyogram:
#'
#' @param methylKit.obj object with chr,start,end positions \code{methylKit.obj}
#' @param cols color vector \code{cols}
#' 
#' @return Karyogram with CpGs plotted across the genome (chr1-22)
#' 
#' @keywords methylBase, methylKit, Karyogram, karyo
#'
#' @export
#' 
#' @examples
#' karyogram<-karyo(methylBase,cols="limegreen")
#' plot(karyogram)

karyo<-function(methylKit.obj,cols){
	require(ggbio)
	require(GenomicRanges)
	require(methylKit)
	require(BSgenome)
	require("BSgenome.Hsapiens.UCSC.hg19")
	
	#get chromosome lenghts for Hg19
	chr.length=seqlengths(Hsapiens)#getchromosomelengths
	#remove M and random chromosomes
	chr.length=chr.length[grep("_|M",names(chr.length),invert=T)]
	
	myKaryo<-GRanges(seqnames=names(chr.length),ranges=IRanges(start=1,
			width=chr.length))
		seqlevels(myKaryo)=names(chr.length)
		seqlengths(myKaryo)=(chr.length)

		g.sites=as(methylKit.obj,"GRanges")
		seqlevels(g.sites,force=TRUE)<-seqlevels(myKaryo)
		seqlengths(g.sites)=(chr.length)

		
		p<-ggplot()+layout_karyogram(myKaryo,cytoband=FALSE)+ 
		theme(axis.text = element_text(size = rel(4)))	+ 
		theme(axis.title = element_blank())
		p+layout_karyogram(c(g.sites),col=cols)
}