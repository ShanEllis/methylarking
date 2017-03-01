#' Manhattan plot of all sites from DMS Analysis
#'
#' @param DMS output dataframe produced by runlimma \code{DMS}
#' @param filename filenames to attach \code{filename}
#' @param sig level of genome-wide significance \code{sig}
#' 
#' @return plot single site analysis output
#' 
#' @keywords limma, manhattan, 
#'
#' @export
#' 
#' @examples
#' manhattan(DMS,"test")

manhattan<-function(DMS, filename, sig=NULL){
	require(methylKit)

	#if necessary, get rid of sex chromosomes
	data<-subset(DMS, DMS$chr!="chrX")
	data<-subset(data, (data$chr!="chrY"))
	data<-subset(data, (data$chr!="chrM"))
	
	data$chr <- gsub("chr","",data$chr)
	data$chr <- as.numeric(data$chr)
	data2 <- data[order(data$chr),]
	data=data2
	data2$meth.diff <- 
	cex.val <- 0.8 - ((abs(data[,"mean.meth.diff"]) < 5 )*.5) + 
	((abs(data[,"mean.meth.diff"]) > 20)*.7)	
	data$cex.val<-cex.val

	# library(GenABEL)
	# z.sq<-(data$Effect / data$StdErr)^2
	# lambda<-estlambda(data$P.Value)
	# cat(paste("Lambda =",lambda,"\n"))
	##working from metal output###
	if (is.null(sig)==TRUE){
		ymin1=round( max(-log10(data$P.Value))+1)
	}
	else{
		ymin1=sig+1
	}
	title=c() 
	jpeg(paste("manhattan_", filename, ".jpeg",sep=""),res=400,width = 37, 
	height = 12,units="cm",type="cairo")
	chr <- c(1:22)
	#Summary statistics
	data$position<-round(data$start,digits=0)
	#print(summary(data))
	print(table(data$chr))
	par(mar=c(5,5,4,2))
	phy.max<-tapply(data$start, data$chr,max,na.rm=T)
	cumlen=0
	for(i in chr){
		cat(paste("Now working on chromosome ",i,"\r"))
		data[data$chr==i,"loc"]<-data[data$chr==i,"position"]+cumlen
		cumlen<-cumlen+phy.max[i]
	}
	phy.med<-tapply(data$loc,data$chr,median,na.rm=T)[chr]
	print(phy.med)
	data$mlgpval<- -log(data[,"P.Value"], base=10)
	plot(data[,"loc"],data[,"mlgpval"],type="n",yaxt="n",xaxt="n",
		xlab="chromosome",
		ylab=expression(-log[10]*P),main=title,
		xlim=c(0,max(data$loc,na.rm=T)),cex.lab=1.5,ylim=c(0,ymin1))
	col=rep(c("black","gray48"),13)
	axis(side=2, at=seq(from=0,to=ymin1,by=1), labels=seq(from=0,to=ymin1,by=1),
		tick=T,cex.axis=1.0,las=1)
	axis(side=1, at=phy.med[c(1:22)], labels=chr[c(1:22)],
		tick=T,cex.axis=1.2,las=1)
	#axis(side=1, at=phy.med[c(20:23)], labels=chr[c(20:23)],
	#	tick=T,cex.axis=1,las=1)
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
		col = "white")
	for(i in chr){
		cat(paste("Now working on chromosome ",i,"\r"))
		#if(which(chr==i) %in% seq(2,30,2)) col="blue" else col="red"
		points(data[data$chr==i,"loc"],data[data$chr==i,"mlgpval"],
			col=col[i],pch=20,cex=data[data$chr==i,"cex.val"])
		#,cex = data[data$chr==i,"cex.val"]
	}
	# add in line for significance
	if (is.null(sig)==FALSE){
		abline(h=sig,lty="dotted",lwd=2,col="chartreuse4")
	}
	
	dev.off()
}