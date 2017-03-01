#' Pie Chart: GeneAnnotation
#'
#' @param methylKit.obj methylBase object \code{methylKit.obj}
#' @param gene.obj GRaangesList w/ promoter,exon,intron,and TSS \code{gene.obj}
#' 
#' @return Pie Chart for genic annotation
#' 
#' @keywords methylBase, methylKit, PieChart, annotation 
#'
#' @export
#' 
#' @examples
#' pie<-GenePie(methylBase,gene.obj)
#' plot(pie)

GenePie<-function(methylKit.obj,gene.obj){
	require(methylKit)
	#annotate CpGs wit hpromoter/exon/intron using annotation data
	g.sites=as(methylKit.obj,"GRanges")
	diffAnn=annotate.WithGenicParts(g.sites,gene.obj)
	slices <- diffAnn@precedence
	lbls <- names(diffAnn@precedence)
	pct <- round(slices/sum(slices)*100)  
	lbls <- paste(lbls, "\n",pct,"%", sep="") 
	pie(slices, labels=lbls, col=rainbow(4),radius=1.0)
}