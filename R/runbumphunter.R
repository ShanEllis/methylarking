#' Run bumphunter to obtain DMRs (parameters stolen straight from bumphunter)
#'
#' @param methyl.mat matrix of samples by methylation \code{DMS}
#' @param covariates matrix (samples in rows, covar in column) \code{covariates}
#' @param pickCutoff Should bumphunter attempt to pick a cutoff using the 
#' permutation distribution?
#' @param pickCutoffQ The quantile used for picking the cutoff using the 
#' permutation distribution.
#' @param maxGap If cluster is not provided this maximum location gap will 
#' be used to define cluster via the clusterMaker function.
#' @param smooth A logical value. If TRUE the estimated profile will be 
#' smoothed with the smoother defined by smoothFunction
#' @param B An integer denoting the number of resamples to use when computing 
#' null distributions. This defaults to 0. If permutations is supplied that 
#' defines the number of permutations/bootstraps and B is ignored. \code{B}
#' @param nullMethod  Method used to generate null candidate regions, must be 
#' one of ‘bootstrap’ or ‘permutation’ (defaults to ‘permutation’). However, 
#' if covariates in addition to the outcome of interest are included in the 
#' design matrix (ncol(design)>2), the ‘permutation’ approach is not recommended. 
#' @param n.core number of cores to run in paralell 

#'  
#' 
#' @return annotated data frame; see bumphunter documentation for further 
#' explanation of output
#' 
#' @keywords bumphunter, DMRs 
#'
#' @export
#' 
#' @examples
#' DMS_annotated=annotateGenes(DMS)

runbumphunter <-function(methyl.mat,treatment, covariates,pickCutoff=TRUE,
	pickCutoffQ=0.975,maxGap=300,smooth=FALSE,B=250,nullMethod="bootstrap",
		n.core=20){
require(bumphunter)
require(doParallel)
registerDoParallel(cores=n.core)

covariates=covariates_N63
ff<-(~treatment+covariates$age+covariates$sex+covariates$site+covariates$SV1+
	covariates$SV2+covariates$SV3+covariates$SV4+covariates$SV5+covariates$SV6+
	covariates$SV7+covariates$SV8+covariates$SV9+covariates$SV10)
designmatrix<-model.matrix(ff,covariates)

#run bumphunter
chr<-grep("chr",unlist(strsplit(rownames(methyl.mat),"[.]")),val=T)
#position
pos<-as.numeric(grep("chr",unlist(strsplit(rownames(methyl.mat),"[.]")),val=T,
	invert=TRUE))
#runit
bump_out<-bumphunterEngine(methyl.mat,design=designmatrix,chr=chr,pos=pos,
	pickCutoff=pickCutoff,pickCutoffQ=pickCutoffQ,maxGap=maxGap,
	smooth=smooth,B=B,nullMethod=nullMethod)
return(bump_out$table)
}