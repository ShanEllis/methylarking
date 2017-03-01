#' DVS Analysis
#'
#' @param methylKit.obj methylBase object \code{methylKit.obj}
#' @param covariates covariate data frame \code{covariates}
#' @param variance cutoff for sites to use based on variance \code{variance}
#' 
#' return single site DVS analysis (using limma's lmFit) 
#' adjusts for covariates
#' and weights methylation values based on log10(coverage)
#' 
#' @keywords methylBase, methylKit, DVS, limma, variance
#'
#' @export
#' 
#' @examples
#' DVS<- runDVS(methylkit.obj,covariates=covariates_N63,variance=0.75)

runDVS<-function(methylKit.obj,covariates){
	require(limma)
	meth = getData(methylKit.obj)
	rownames(meth) <- paste(meth$chr,meth$start,sep=".")
	mat=percMethylation(methylKit.obj)
	rownames(mat)<-paste(methylKit.obj$chr,methylKit.obj$start,sep=".")

#get coverage
coverage<-meth[,grep("coverage", colnames(methylKit.obj))]
rownames(coverage)<-paste(methylKit.obj$chr,methylKit.obj$start,sep=".")


# remove per-CpG outliers
 CpG=mat
# for(i in 1:dim(CpG)[2]){
# CpG[,i][(abs(CpG[,i] - mean(CpG[,i], na.rm = T)) > 
	# 3 *sd(CpG[,i], na.rm = T))] <- NA
# }
# CpG=t(CpG)

#variance across sites
vars <- apply(CpG, 1, var, na.rm=TRUE)
#only analyze variant sites
#variance=variance
#CpG <- subset(CpG, vars > quantile(vars,1-variance))  #769,065

weights=log10(coverage)
#get same sites in CpG and coverage
coverage <- subset(coverage, rownames(coverage) %in% rownames(CpG))

#get measure of variability at each CpG
#A measure of variability is calculated for each CpG in each sample
#by subtracting out the group mean and taking the absolute or squared deviation.
means <- apply(CpG,1,mean,na.rm=TRUE)
varib <- abs(CpG-means)

mean.meth.diff <- apply(varib, 1, mean.meth,treatment=methylKit.obj@treatment)

### run the model
design=model.matrix(~methylKit.obj@treatment + as.matrix(covariates))

VMSfit_log10cov <- lmFit(varib,design=model.matrix(~methylKit.obj@treatment
		+as.matrix(covariates)), weights=log10(coverage))
	#get residuals for future plotting of out put data
	VMSnull_log10cov_null <-  lmFit(varib,design=model.matrix(
		~as.matrix(covariates)), weights=log10(coverage))
	VMSfit_log10cov_fit_resid<-residuals(VMSnull_log10cov_null,varib)
	#get mean methylation (residuals)
	mean.meth.resid<-apply(VMSfit_log10cov_fit_resid, 1, mean, na.rm=TRUE)
	fit <- eBayes(VMSfit_log10cov)
	mean.meth.diff.resid <- fit$coefficients[,2]
	meth_out<-cbind(mean.meth.resid,mean.meth.diff, mean.meth.diff.resid)
	VMS_fit<-topTable(fit,coef=2,number=nrow(varib))
	#get in correct order
	meth_out_ordered <- meth_out[match(rownames(VMS_fit), rownames(meth_out)),]
	#combine data
	VMS_fit <-cbind(VMS_fit,meth_out_ordered)
	#get chr, start, end positions for graphing later on
	sites<-as.data.frame(matrix(unlist(strsplit(rownames(VMS_fit),"[.]")), 
		ncol=2, byrow=TRUE))
	colnames(sites)<-c("chr","start")
	sites$chr<-sites[,1]
	sites$start<-as.numeric(as.character(sites[,2]))
	sites$end<-as.numeric(as.character(sites[,2]))
	VMS_fit <- cbind(sites,VMS_fit)

	#sort output
	VMS_fit_sorted <- VMS_fit[order(VMS_fit$P.Value),]
	return(VMS_fit_sorted)
}