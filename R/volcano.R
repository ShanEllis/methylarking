#' Volcano plot of all sites from DMS Analysis
#'
#' @param DMS output dataframe produced by runlimma \code{DMS}
#' 
#' @return volcano plot for single site analysis output
#' 
#' @keywords limma, volcano 
#'
#' @export
#' 
#' @examples
#' volcano(DMS,"test")

volcano <- function(DMS,filename){
	jpeg(paste("volcano_", filename, ".jpeg",sep=""), res=600,width = 16, 
		height = 16,unit="cm",type="cairo")
	plot(-log10(DMS$P.Value)~DMS$mean.meth.diff,
		xlab="Mean Methylation Difference",
		ylab=expression(-log[10]*P),
		main="Volcano Plot")
	dev.off()	
}
