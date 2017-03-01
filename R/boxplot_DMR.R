#' Boxplot of regions from output from DMR Analysis
#'
#' @param annotated DMR table output from runbumphunter  \code{DMR}
#' @param methyl.mat matrix with methylation values \code{methyl.mat}
#' @param treatment vector with case control status \code{treatment}
#' @param n number of sites you'd like to output \code{n}
#' 
#' @return plot single site analysis output
#' 
#' @keywords bumphunter, boxplot, DMR
#'
#' @export
#' 
#' @examples
#' boxplot_DMR(DMR,methyl.mat,treatment,filename,n=15)

boxplot_DMR<-function(DMR,methyl.mat,treatment,filename,n=15){
	require(bumphunter)
	require(ggbio)
	require(ggplot2)
	
	DMR <- subset(DMR, DMR$L>1)
	png(paste("Boxplot_bumphunter_",filename, ".png",sep=""),height=1400,
		width=1400,	type=c("cairo"))
	par(mfrow=c(3,ceiling(n/3)))
	#chr
	chr<-grep("chr",unlist(strsplit(rownames(methyl.mat),"[.]")),val=T)
	#position
	pos<-as.numeric(grep("chr",unlist(strsplit(rownames(methyl.mat),"[.]")),
		val=T,invert=TRUE))
	plots <- list()  
	for(i in 1:n){
		indexstart<-DMR$indexStart[i]
		indexend<-DMR$indexEnd[i]
		grabBeta<- methyl.mat[indexstart:indexend,]
		#only plot if bigger than a single site
		
			formbeta<-c()
			for(j in 1:nrow(grabBeta)){
				tempbeta<-grabBeta[j,]
				formbeta<-c(formbeta,tempbeta)
			}
			formbeta<-as.numeric(formbeta)
			reppos<-rep(pos[indexstart:indexend],each=1)
			xmin<-min(reppos)-25
			xmax<-max(reppos)+25
			status<-rep( treatment,nrow(grabBeta))
			reppos.real.shift<-ifelse(status=="0",reppos+0.2,reppos)
			toplot<-data.frame((formbeta),reppos.real.shift,status)
			rownames(toplot)<-c()
			colnames(toplot)<-c("Beta","Position","Status")
			p<-ggplot(data=toplot,aes(x=(Position),y=(Beta),color=factor(Status)))+
				geom_point(size=1.5)+scale_colour_manual(values=c("black","red"))+
				theme(legend.title=element_blank(),panel.background=element_blank())+ 
				theme(legend.text = element_text(size = 20))+
				ylab("Percent Methylation(%)")+
				xlab(paste("Position on",DMR$chr[i],sep=" "))+
				scale_x_continuous(limits=c(xmin,xmax))+
				stat_smooth(method="loess",size=2,se=FALSE)+	#ggtitle(paste(paste(DMR$chr[i],DMR$start[i],DMR$end,sep="."),DMR_annotated$name[i],DMR_annotated$description[i],sep="\n"))+
				theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+ 
				theme(axis.title.x = element_text(size = rel(1.2), angle = 00)) + 
				theme(text = element_text(size=15))+
				theme(axis.text.y = element_text(size=17))+
				theme(axis.text.x = element_text(size=13))
			assign(paste("p",i,sep=""),p)
			plots[[i]] <- p	# add each plot into plot list
		}
	multiplot(plotlist = plots, cols = 3)
	dev.off()
}	
