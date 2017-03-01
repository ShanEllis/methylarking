#' Pie Chart: GeneAnnotation
#'
#' @param methylKit.obj obj w/ chr, start, end positions \code{methylKit.obj}
#' @param gene.obj GRaangesList w/ promoter,exon,intron,and TSS \code{gene.obj}
#' 
#' @return Pie Chart for genic annotation
#' 
#' @keywords methylBase, methylKit, PieChart, annotation 
#'
#' @export
#' 
#' @examples
#' GenePie(methylKit.obj,gene.obj)
GenePie<-function(methylKit.obj,gene.obj){
	require(methylKit)
	#annotate CpGs with promoter/exon/intron using annotation data
	g.sites=as(methylKit.obj,"GRanges")
	diffAnn=annotate.WithGenicParts(g.sites,gene.obj)
	slices <- diffAnn@precedence
	lbls <- names(diffAnn@precedence)
	pct <- round(slices/sum(slices)*100)  
	lbls <- paste(lbls, "\n",pct,"%", sep="") 
	pie(slices, labels=lbls, col=rainbow(length(diffAnn@precedence)),radius=1.0)
	Gene_annot<-slices
	return(Gene_annot)
}


#' Pie Chart: CGI Annotation
#'
#' @param methylKit.obj obj with chr,start,end positions \\code{methylKit.obj}
#' @param cpg.obj SimpleGenomicRangesList w/ CGI info \code{cpg.obj}
#' 
#' @return Pie Chart for CGI annotation
#' 
#' @keywords methylBase, methylKit, PieChart,CGI, annotation 
#'
#' @export
#' 
#' @examples
#' CpGPie(methylKit.obj,cpg.obj)
CpGPie<-function(methylKit.obj,cpg.obj){
	require(methylKit)
	#annotate CpGs with CGI info using annotation data
	g.sites=as(methylKit.obj,"GRanges")
	diffCpGann=annotate.WithFeature.Flank(g.sites,cpg.obj$CpGi,cpg.obj$shores,
	feature.name="CpGi",flank.name="shores")
	slices <- diffCpGann@precedence
	lbls <- names(diffCpGann@precedence)
	pct <- round(slices/sum(slices)*100)  
	lbls <- paste(lbls, "\n",pct,"%", sep="") 
	pie(slices, labels=lbls, col=terrain.colors(length(diffCpGann@precedence)),
	radius=1.0)
	CpG_annot<-slices
	return(CpG_annot)
}

