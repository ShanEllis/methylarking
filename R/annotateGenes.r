#' Annotate genomic position using bumphunter's matchGenes
#'
#' @param DMS obj w/ chr, start, end positions \code{DMS}
#' 
#' @return annotated data frame; see bumphunter documentation for further 
#' explanation of output
#' 
#' @keywords annotate, bumphunter, matchGenes 
#'
#' @export
#' 
#' @examples
#' DMS_annotated=annotateGenes(DMS)

annotateGenes <-function(DMS){
	require("bumphunter")
	#remove X,Y,M
	DMS_noXY <- subset(DMS, DMS$chr %in% names(summary(DMS$chr))[1:22])
	DMS_annotated <- matchGenes(DMS_noXY)
	DMS_annotated = cbind(DMS_noXY, DMS_annotated)
	return(DMS_annotated)
}