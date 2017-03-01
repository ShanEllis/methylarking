#' Boxplot of single sites from output from DMS Analysis
#'
#' @param DMS output dataframe produced by runlimma \code{DMS}
#' @param methyl.mat matrix with methylation values \code{methyl.mat}
#' @param treatment vector with case control status \code{treatment}
#' @param n number of sites you'd like to output \code{n}
#' 
#' @return plot single site analysis output
#' 
#' @keywords limma, boxplot, 
#'
#' @export
#' 
#' @examples
#' boxplot_DMS(DMS,methyl.mat,treatment,filename,n=15)


boxplot_DMS<-function(DMS,methyl.mat,treatment,filename,n=15){
	png(paste("Boxplot_", filename, ".png",sep=""),height=1000,width=1400,
		type=c("cairo"))

	par(mfrow=c(3,ceiling(n/3)))
	for (i in 1:n){
	boxplot(as.numeric(methyl.mat[rownames(DMS)[i],])~as.factor(treatment),
		pch=19,col=c("grey48","red"),main=paste("","","",sep="\n"),
		ylab="",xaxt="n")
		axis(side=1,at=c(1,2),labels=c("controls","cases"),lty=NULL)
		#Calculate mean
		means<-by(na.omit(methyl.mat[rownames(DMS)[i],]),
			treatment[!is.na(methyl.mat[rownames(DMS)[i],])],mean)
		points(c(1,2),means,pch=21,cex=1.5,bg="blue")
		par(new=T)
		stripchart(as.numeric(methyl.mat[rownames(DMS)[i],])
		~as.factor(treatment),col="grey48",
		main=paste(rownames(DMS)[i],paste("log10pval: ",
		round(-log10(DMS$P.Value[i]),2),sep=""),
		paste("mean.diff: ",round(DMS$mean.meth.diff[i],2),sep=""),sep="\n"),
		vertical=TRUE,ylab="residuals",xaxt="n")
	}
dev.off()
}