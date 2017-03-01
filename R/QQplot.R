#' QQ plot of all sites from DMS Analysis broken down by Mean Methylation
#'
#' @param DMS output dataframe produced by runlimma \code{DMS}
#' @param filename filenames to attach \code{filename}
#' 
#' @return QQplot single site analysis output
#' 
#' @keywords limma, QQplot 
#'
#' @export
#' 
#' @examples
#' QQplot(DMS,"test")

#QQ function
obs <- function(pval){
observed_data <-sort(as.numeric(paste(pval)))
lobs_data = -(log10(observed_data))
return(lobs_data)
}
exps <- function(pval){
observed_data <-  sort(as.numeric(paste(pval)))
lobs_data <- -(log10(observed_data))
expected_data <- c(1:length(observed_data)) 
lexp_data <- -(log10(expected_data / (length(expected_data)+1)))
return(lexp_data)
}
QQplot_sub <- function(DMS,filename){
	cols=rainbow(10)
	lbls=c("[0-10]","(10-20]","(20-30]","(30-40]","(40-50]","(50-60]","(60-70]",
		"(70-80]","(80-90]","(90-100]")
	j=1
	jpeg(paste("QQPlot_meanmeth_", filename, ".jpeg",sep=""), res=600,
		width = 16, height = 16, unit="cm",type="cairo")
	par(oma=c(1,2.5,0,0),mgp=c(2.5,1,0) )
	plot(c(0,10), c(0,10), col="black", lwd=1, type="l", 
		xlab=expression(Expected -log[10]*(pval)), 
		ylab=expression(Observed -log[10]* (pval)), cex.lab=1.7, xlim=c(0,10), 
		ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l",main="QQ plot", 
		cex.main=1.5)
	for(i in c(0,10,20,30,40,50,60,70,80,90)){
		if( i==0|i==5|i==90|i==95){
			a<-subset(DMS,DMS$AveExpr > i & DMS$AveExpr <= i+5 )
			points(exps(a$P.Value), obs(a$P.Value), pch=23, cex=.5, 
				col=cols[j], bg=cols[j])
			j=j+1
			
		}
		else{
			a<-subset(DMS,DMS$AveExpr > i & DMS$AveExpr <= i+10 )
			points(exps(a$P.Value), obs(a$P.Value), pch=23, cex=.5, 
				col=cols[j], bg=cols[j])
			j=j+1
		}
		legend(5.5,4,c(lbls), cex=1.0, col=cols, pch=18,ncol=2)

	}
	dev.off()
}

#' QQ plot of all sites from DMS Analysis
#'
#' @param DMS output dataframe produced by runlimma \code{DMS}
#' @param filename filenames to attach \code{filename}
#' @param cols color of points \code{cols}

#' 
#' @return QQplot single site analysis output
#' 
#' @keywords limma, QQplot 
#'
#' @export
#' 
#' @examples
#' QQplot(DMS,"test",cols="red")
QQplot <- function(DMS,filename,cols="dodgerblue"){
	require('GenABEL')
	lbls=paste(filename,"\n","l=",round(estlambda(DMS$P.Value)$estimate,2),sep="")
	jpeg(paste("QQPlot_pval",pvalue,"_", filename, ".jpeg",sep=""), res=600,
		width = 16, height = 16,unit="cm",type="cairo")
	par(oma=c(1,2.5,0,0),mgp=c(2.5,1,0) )
	plot(c(0,10), c(0,10), col="black", lwd=1, type="l", 
		xlab=expression(Expected -log[10]*(pval)), 
		ylab=expression(Observed -log[10]* (pval)), 
		cex.lab=1.0, xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", 
		bty="l",main="QQ plot", cex.main=1.5)
	points(exps(DMS$P.Value), obs(DMS$P.Value), pch=23, cex=.5, col=cols, 
		bg=cols)
	legend(5.5,4,c(lbls), cex=1.0, col=cols, pch=18,ncol=2,box.col=NA)
	dev.off()
}
