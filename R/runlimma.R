#' DMS Analysis
#'
#' @param methylKit.obj methylBase object \code{methylKit.obj}
#' @param covariates covariate data frame \code{covariates}
#' @param variance cutoff for sites to use based on variance \code{variance}
#' 
#' @return single site analysis (using limma's lmFit) -- adjusts for covariates
#' and weights methylation values based on log10(coverage)
#' 
#' @keywords methylBase, methylKit, DMS, limma, variance
#'
#' @export
#' 
#' @examples
#' DMS<- runlimma(methylkit.obj,covariates,variance=0.75)

runlimma<-function(methylKit.obj,covariates,variance){
	require(limma)
	#get CpG counts (rowpercent methylation)
	mat=percMethylation(methylKit.obj)
	rownames(mat)<-paste(methylKit.obj$chr,methylKit.obj$start,sep=".")

	
	meth <- getData(methylKit.obj)
	coverage<-meth[,grep("coverage", colnames(methylKit.obj))]
	rownames(coverage)<-paste(methylKit.obj$chr,methylKit.obj$start,sep=".")
		
	# remove X,Y,M chromosome
	# chrY<-grep("chrY",rownames(mat))
	# chrX<-grep("chrX",rownames(mat))
	# chrM<-grep("chrM",rownames(mat))
	
	# chrs<-c(chrX,chrY,chrM)
	# mat_noXY<-mat[-chrs,]	
	# coverage_noXY <- coverage[-chrs,]
	# print(dim(coverage_noXY))
	# colnames(coverage_noXY) <- colnames(mat)
		
	# remove per-CpG outliers
	CpG=t(mat)
	for(i in 1:dim(CpG)[2]){
	CpG[,i][(abs(CpG[,i] - mean(CpG[,i], na.rm = T)) > 
		3 *sd(CpG[,i], na.rm = T))] <- NA
	}
	CpG=t(CpG)
	
	#variance across sites
	vars <- apply(CpG, 1, var, na.rm=TRUE)

	#only analyze variant sites
	CpG <- subset(CpG, vars > quantile(vars,1-variance)) 
	#mean.methylation.difference
	mean.meth.diff <- apply(CpG, 1, mean.meth,treatment=methylKit.obj@treatment)
	#get coverage
	coverage <- subset(coverage, rownames(coverage) %in% rownames(CpG))


# run limma's lmFit, weight by coverage; autosomal sites
# lmFit: Fit linear model for each gene given a series of arrays
# This function fits multiple linear models by weighted or generalized least
# squares. 
	DMRfit_log10cov <- lmFit(CpG,design=model.matrix(~methylKit.obj@treatment
		+as.matrix(covariates)), weights=log10(coverage))
	#get residuals for future plotting of out put data
	DMRnull_log10cov_null <-  lmFit(CpG,design=model.matrix(
		~as.matrix(covariates)), weights=log10(coverage))
	DMRfit_log10cov_fit_resid<-residuals(DMRnull_log10cov_null,CpG)
	#get mean methylation (residuals)
	mean.meth.resid<-apply(DMRfit_log10cov_fit_resid, 1, mean, na.rm=TRUE)
	fit <- eBayes(DMRfit_log10cov)
	Z<-fit$coefficients[,2] / fit$sigma
	mean.meth.diff.resid <- fit$coefficients[,2]
	meth_out<-cbind(mean.meth.resid,mean.meth.diff, mean.meth.diff.resid)
	DMR_fit<-topTable(fit,coef=2,number=nrow(CpG))
	#get in correct order
	meth_out_ordered <- meth_out[match(rownames(DMR_fit), rownames(meth_out)),]
	#combine data
	DMR_fit <-cbind(DMR_fit,meth_out_ordered)
	#get chr, start, end positions for graphing later on
	sites<-as.data.frame(matrix(unlist(strsplit(rownames(DMR_fit),"[.]")), 
		ncol=2, byrow=TRUE))
	colnames(sites)<-c("chr","start")
	sites$chr<-sites[,1]
	sites$start<-as.numeric(as.character(sites[,2]))
	sites$end<-as.numeric(as.character(sites[,2]))
	DMR_fit <- cbind(sites,DMR_fit)

	#sort output
	DMR_fit_sorted <- DMR_fit[order(DMR_fit$P.Value),]
	return(DMR_fit_sorted)
}


#' Get Hypermethylated Sites
#'
#' @param DMS output from DMS Analysis \code{DMS}
#' @param methdiff difference in methylation \code{methdiff}
#' @param pval pvalue cutoff \code{pval}
#' 
#' @return data frame with hypermethylated sites
#' 
#' @keywords methylBase,hyper, hypermethylation
#'
#' @export
#' 
#' @examples
#' DMS_hyper <- hyper(DMS,methdiff=0,pval=1)
hyper<-function(DMS,methdiff=0,pval=NULL){
	hyper=subset(DMS, DMS$mean.meth.diff.resid > methdiff & DMS$P.Value < pval)
	return(hyper)
}


#' Get Hypomethylated Sites
#'
#' @param DMS output from DMS Analysis \code{DMS}
#' @param methdiff difference in methylation \code{methdiff}
#' @param pval pvalue cutoff \code{pval}
#' 
#' @return data frame with hypomethylated sites
#' 
#' @keywords methylBase,hypo, hypomethylation
#'
#' @export
#' 
#' @examples
#'  DMS_hypo <- hypo(DMS,methdiff=0,pval=1)
hypo<-function(DMS,methdiff=0,pval=NULL){
	hypo=subset(DMS, DMS$mean.meth.diff.resid < -methdiff & DMS$P.Value < pval)
	return(hypo)
}
