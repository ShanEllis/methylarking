#' Plot number of sites by chromosome
#'
#' @param DMS output dataframe produced by runlimma \code{DMS}
#' @param pval filenames to attach \code{pval}
#' 
#' @return QQplot single site analysis output
#' 
#' @keywords limma, QQplot 
#'
#' @export
#' 
#' @examples
#' sitesperchr(DMS,pval=pvalue)

sitesperchr <- function(DMS,pval){
	png(paste("sitesperchr_", pval, ".jpeg",sep=""),width=1500,height=500,
	type="cairo")
	chrOrder<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
	"chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18",
	"chr19","chr20","chr21","chr22","chrX","chrY","chrM")
	chrs<-factor(DMS$chr,levels=chrOrder, ordered=TRUE)
	barplot(t(summary(chrs)),main="Number of CpG sites per chromosome")
	dev.off()
}