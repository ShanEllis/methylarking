#' Boxplot of single sites from output from DMS Analysis
#'
#' @param methylKit.obj output methylkit object \code{methylKit.obj}
#' @param covariates covariate matrix \code{covariates}
#' @param variance what proportion of top variant sites to include \code{}
#' 
#' @return Rdata object for bootstrap
#' 
#' @keywords boostrap, limma, DMS 
#'
#' @export
#' 
#' @examples
#' bootstrap(meth_normalized_N63_count100_20L, covariates=covariates_N63, variance=0.75,JOBID=2)

bootstrap<-function(methylKit.obj,covariates,variance,JOBID){
	coef=2
	B=1
	require(limma)
	require(methylKit)
	treatment=methylKit.obj@treatment
	#get CpG counts (percent methylation)
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
	
	#get coverage
	coverage <- subset(coverage, rownames(coverage) %in% rownames(CpG))


# bootstrap	
# run limma's lmFit, weight by coverage; autosomal sites
# lmFit: Fit linear model for each gene given a series of arrays
# This function fits multiple linear models by weighted or generalized least
# squares. 

bootIndexes<-replicate(B, sample(1:ncol(mat),replace=TRUE),simplify=TRUE)

#fit the alternative model
message('fiting the alternative model')
	DMRfit_log10cov <- lmFit(CpG,design=model.matrix(~treatment +as.matrix(covariates)), weights=log10(coverage))
	#get residuals from full model
	resids_full <-  residuals(DMRfit_log10cov,CpG)
	
	full = eBayes(DMRfit_log10cov)
	full = topTable(full,coef=2,nrow(full))
	
	
#2 fit the null model 	
message('fiting the null model')

	DMRfit_null <- lmFit(CpG,design=model.matrix(~as.matrix(covariates)), weights=log10(coverage))
	estimates_null <-DMRfit_null$coefficients[,2:ncol(DMRfit_null$coefficients)]
	
#3 Resample residuals
message('resampling residuals from the alternative model')
	resids_resampled <- resids_full[,bootIndexes]

# generate pseudo-null data (combine fitts from null model with resampled residuals from akterntaive/full model
message('generating pseudonull data')

	pseudonull=matrix(NA,nrow=nrow(resids_resampled), ncol=ncol(resids_resampled))
	for (i in 1:nrow(resids_resampled)){
		pseudonull[i,] <- resids_resampled[i,]- estimates_null[i,2]*(covariates$age-mean(covariates$age)) - (estimates_null[i,3]*(covariates$sex))- (estimates_null[i,4]*(covariates$site))- estimates_null[i,4]*(covariates$SV1-mean(covariates$SV1))- estimates_null[i,5]*(covariates$SV2-mean(covariates$SV2))- estimates_null[i,6]*(covariates$SV3-mean(covariates$SV3))- estimates_null[i,6]*(covariates$SV4-mean(covariates$SV4))- estimates_null[i,7]*(covariates$SV5-mean(covariates$SV5))- estimates_null[i,8]*(covariates$SV6-mean(covariates$SV6))- estimates_null[i,9]*(covariates$SV7-mean(covariates$SV7))- estimates_null[i,10]*(covariates$SV8-mean(covariates$SV8))- estimates_null[i,11]*(covariates$SV9-mean(covariates$SV9))- estimates_null[i,12]*(covariates$SV10-mean(covariates$SV10))
	}
#4 refit full models with pseudonull as independent 
message('refit full model with pseudonull')
	DMRfit_pseudonull <- lmFit(pseudonull,design=model.matrix(~treatment +as.matrix(covariates)), weights=log10(coverage))
	
	pseudo = eBayes(DMRfit_pseudonull)
	#get methylation difference (From betas)
	differ=pseudo$coefficients[,"treatment"]
	#get out pvalues
	pseudo = topTable(pseudo,coef=coef,nrow(pseudonull))
	fit=pseudo
	pseudo = pseudo$P.Value
	
	names(pseudo) = rownames(CpG)
	names(differ) = rownames(CpG)
	
message('saving output')
	save(pseudo,bootIndexes,fit, differ ,file=paste("Bootstrap",JOBID,".rda",sep=""))
}