#' Mean Methylation Difference -- by case-control
#'
#' @param methyl.mat methylation matrix \code{methyl.mat}
#' @param treatment vector with case control information \code{treatment}

#' 
#' @return mean methylation differences between cases and controls
#' 
#' @keywords methylBase, methylKit, methyl.diff, mean methylation difference
#'
#' @export
#' 
#' @examples
#' #' mean.meth.diff=apply(methyl.obj,1,mean.meth)

mean.meth <- function(X,treatment){
	meth=na.omit(X)
	treat=treatment[!is.na(X)]
	means=by(meth,treat,mean)
	return(means[2]-means[1])
}