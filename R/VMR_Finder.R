###############################
# VMR Finder
# http://biostatistics.oxfordjournals.org/content/early/2011/06/17/biostatistics.kxr013.full
# Andrew Jaffe
# Updated 10/24/2011
######################

# p = methylation matrix of Beta values: M probes down the rows by N samples across columns
# coi = covariate of interest. a vector of 1's if one-sample VMRs, or respective 
#		factor labels for 2 sample dVMRs'
# G = vector (length M) of median probe intensities. Or a matrix of intensities for each sample.
#		This is the green channel for CHARM
# pns = vector (length M) of probe group IDs (for smoothing and region finding).
#		generated from clusterMaker(chr,pos,maxGap). 
# chr = vector (length M) of chromosome information for each probe
# pos = vector (length M) of position information for each probe
# WBP = smoothing span parameter. For each probe group: 
#		span = WBP/median_probe_spacing/num_probes_in_group
# cutoffQuant = quantile of smoothed variance statistic used to define regions
# n_boot = number of bootstrap permutations 
# batches = NULL for no batch correction, or batch variable (like processing date)
#			If provided, strip plots with guide the number of batches to remove, 
#			with user input

vmrFinder = function(p, coi, G, pns, chr, pos,WBP=600,cutoffQuant=0.99,
	n_boot = 200, batches=NULL) {
	
	if(length(grep("chr", chr)) == 0) paste("chr",chr,sep="")
	
	# load libraries
	
	pckg = try(require(corpcor))
	if(!pckg) {
		cat("Installing 'corpcor' from CRAN\n")

		getPckg(corpcor)
		require(corpcor)
	}

	pckg = try(require(limma))
	if(!pckg) {
		cat("Installing 'limma' from Bioconductor\n")

		source("http://bioconductor.org/biocLite.R")
		biocLite("limma")
		require("limma")
	}
	
	pckg = try(require(Biobase))
	if(!pckg) {
		cat("Installing 'Biobase' from Bioconductor\n")
		source("http://bioconductor.org/biocLite.R")
		biocLite("Biobase")
		require("Biobase")
	}

	if(is.matrix(G)) G = rowMedians(G)
	
	# drop sex chr
	Index = which(!chr %in% c("chrX","chrY","chrM")) #added mitochondrial chr 
	p=p[Index,]; chr=chr[Index]; pos=pos[Index];pns=pns[Index]
	if(!is.null(G)) G = G[Index]

	
	type = length(unique(coi)) 
	
	stopifnot(type %in% c(1,2)) 
	
	vmrType = ifelse(type == 1, "vmr","dvmr")
	

	# one sample VMRs
	if(type == 1) {
		
		vmrName = unique(coi)
		if(is.character(vmrName)) vmrName = paste("group",vmrName,sep="")
		
		cleanp = p
		
		if(!is.null(batches)) {
			y = logit(p)
			d = factor(batches)
			mm=rowMedians(y)
			svd=fast.svd(y-mm, tol=0)
			
			len = ceiling(sqrt(ncol(y)))
			nr = min(c(len,4))
			mypar(nr, nr)
			
			for(k in 1:min(c(nr^2,ncol(y)))) {
				stripchart(svd$v[,k] ~ d,
					method = "jitter", vertical = TRUE,
					ylab = paste("pc",k))
			}
			
			cat("How many PCs to remove? Enter Number: ")
			nsv = scan(n=1)
				
			dd=svd$d;dd[c(1:as.numeric(nsv))]=0 # remove PCs 
			y=svd$u%*%diag(dd)%*%t(svd$v)
			cleanp = ilogit(y + mm)
		}
		
	
		# smooth out the effect of G
		s = rowMads(cleanp)
		fit = loessFit(log(s), log2(G), span = 0.05)
		newmad = fit$residuals

		sf = loessSmoothStat(newmad, pns,pos,WBP = WBP)
		#cutoff = quantile(sf, cutoffQuant) 		#original code
		cutoff = quantile(sf, cutoffQuant,na.rm=TRUE) #SE edit
		vmrs = regionFinderPos(sf, pns, chr,
					pos, cutoff = cutoff)
		
		stat = newmad
	} 

	if(type == 2) {
		cat("dVMR Finder\n")
		Indexes = splitit(coi)
		coiNames = names(Indexes)
		vmrName = paste(coiNames,collapse = "_minus_")

		cleanp = p
		pList = list()

		if(!is.null(batches)) {
		
			cat("Look at boxplots and remove batches at the prompts\n")

		
			for(i in seq(along=Indexes)) {
				y = logit(p[,Indexes[[i]]])
				# mod =  matrix(rep(1), nc = 1, nr = ncol(y))
				# nsv = num.sv(y,mod,method="be")$n
				d = factor(batches[Indexes[[i]]])
				mm=rowMedians(y)
				svd=fast.svd(y-mm, tol=0)
			
				len = ceiling(sqrt(ncol(y)))
				nr = min(c(len,4))
				mypar(nr, nr)
				
				for(k in 1:min(c(nr^2,ncol(y)))) {
					stripchart(svd$v[,k] ~ d,
					# stripchart(svd$v[,k],
						method = "jitter", vertical = TRUE,
						ylab = paste("pc",k))
				}
				
				# input = file("stdin")
				# nsv = print(readLines(n=1, ok=FALSE))
	
				# close(input)
				# nsv = readline(prompt="How many PCs to remove? Enter Number: ")
				cat("How many PCs to remove? Enter Number: ")
				nsv = scan(n=1)
				
				dd=svd$d;dd[c(1:as.numeric(nsv))]=0 # remove PCs 
				y=svd$u%*%diag(dd)%*%t(svd$v)
				pList[[i]] = ilogit(y + mm)
				cleanp[,Indexes[[i]]] = pList[[i]]
				
			}
		} else {
			pList[[1]] = cleanp[,Indexes[[1]]]
			pList[[2]] = cleanp[,Indexes[[2]]]
		}
			

		mm = rowMedians(cleanp)
		p1 = pList[[1]] - mm
		p2 = pList[[2]] - mm 
		
		
		# smooth out the effect of G
		s1 = rowMads(p1)
		fit = loessFit(log(s1), log2(G), span = 0.05)
		newmad1 = fit$residuals

		s2 = rowMads(p2)		
		fit = loessFit(log(s2), log2(G), span = 0.05)
		newmad2 = fit$residuals

		lev = newmad1 - newmad2
		levf = loessSmoothStat(lev, pns, pos, WBP=WBP)
		cat("\n")
		cutoff = quantile(abs(levf),cutoffQuant,na.rm=TRUE) #added na.rm SE
		cat("Finding VMRs\n")
		vmrs = regionFinder(levf, pns, chr,
			pos, cutoff = cutoff,verbose=FALSE)
		stat = lev
	}	
	
	# generic boostrapping
	# null=arimaBootSmoothStat(stat,n_boot,pns,pos,WBP) 
	# nullareas = bootArea(null,cutoff=cutoff, pns, type = vmrType)
	# vmrs$pval = edge.pvalue(vmrs$area, nullareas)
	# vmrs$fdr = edge.qvalue(vmrs$pval,lambda=0)$qval

	# sig = vmrs[vmrs$fdr < 0.05,]
	
	#return(sig)
	return(vmrs)
}

logit=function(x) log(x)-log(1-x)
ilogit=function(x) 1/(1+exp(-x))

# generates null areas
bootArea <- function(null, cutoff , pns, type="vmr") {

	area.dat <- list()
	for(j in 1:ncol(null)) {
		cat(j,",")
		if(type == "vmr") {
			hold.dat <- regionFinderPos(null[,j], pns,
					chr,pos, cutoff=cutoff,verbose=F)
		}
		if(type == "dmr") {
			hold.dat <- regionFinder(null[,j], pns,
					chr,pos, cutoff=cutoff,verbose=F)
		}
		if(type == "dvmr") {
			hold.dat <- regionFinder(null[,j], pns,
					chr,pos, cutoff=cutoff,verbose=F)
		}
		
		area.dat[[j]] <- hold.dat$area
	}
	area.dat <- unlist(area.dat)
	return(area.dat)
}

# via Rafael Irizarry
rowMads <- function(x,constant = 1.4826)constant*Biobase::rowMedians(abs(x-Biobase::rowMedians(x,na.rm=TRUE)),na.rm=TRUE) #allowed for NAs in computation

# via Rafael Irizarry
getSegments <- function(x,factor,cutoff=quantile(abs(x),0.99),verbose=TRUE){

  Indexes=split(seq(along=x),factor)
  regionID=vector("numeric",length(x))
  LAST = 0
  
  segmentation = vector("numeric", length(x))
  type = vector("numeric", length(x))

  for (i in seq(along = Indexes)) {
    if (verbose) if (i%%1000 == 0) cat(".")
    Index = Indexes[[i]]
    y = x[Index]
    z = sign(y) * as.numeric(abs(y) > cutoff)
    w = cumsum(c(1, diff(z) != 0)) + LAST
    segmentation[Index] = w
    type[Index] = z
    LAST = max(w)
  }
  if(verbose) cat("\n")
  ##add a vector of the pns
  res=list(upIndex=split(which(type>0),segmentation[type>0]),
    dnIndex=split(which(type<0),segmentation[type<0]),
    zeroIndex=split(which(type==0),segmentation[type==0]))
  names(res[[1]])<-NULL
  names(res[[2]])<-NULL
  names(res[[3]])<-NULL
  return(res)
}

# via Rafael Irizarry
clusterMaker <- function(chr,pos,order.it=TRUE,maxGap=300){
  nonaIndex=which(!is.na(chr) & !is.na(pos))
  Indexes=split(nonaIndex,chr[nonaIndex])
  clusterIDs=rep(NA,length(chr))
  LAST=0
  for(i in seq(along=Indexes)){
    Index=Indexes[[i]]
    x=pos[Index]

    if(order.it){ Index=Index[order(x)];x=pos[Index] }

    y=as.numeric(diff(x)>maxGap)
    z=cumsum(c(1,y))
    clusterIDs[Index]=z+LAST
    LAST=max(z)+LAST
  }
  clusterIDs
}  

# via Rafael Irizarry
##you can pass cutoff through the ...
regionFinder<-function(x,regionNames,chr,position,y=x,
                       summary=mean,ind=seq(along=x),order=TRUE,oneTable=TRUE,
                       ...){

  Indexes=getSegments(x[ind],regionNames[ind],...)
  
  res=vector("list",2)
  for(i in 1:2){
    res[[i]]=data.frame(chr=sapply(Indexes[[i]],function(Index) chr[ind[Index[1]]]),
         start=sapply(Indexes[[i]],function(Index) min(position[ind[Index]])),
         end=sapply(Indexes[[i]],function(Index) max(position[ind[Index]])),
         value=sapply(Indexes[[i]],function(Index) summary(y[ind[Index]])),
         area=sapply(Indexes[[i]],function(Index) abs(sum(y[ind[Index]]))),
         pns=sapply(Indexes[[i]],function(Index) regionNames[ind[Index]][1]),
         indexStart=sapply(Indexes[[i]],function(Index) min(ind[Index])),
         indexEnd=sapply(Indexes[[i]],function(Index) max(ind[Index])))
    res[[i]]$L=res[[i]]$indexEnd-res[[i]]$indexStart+1
  }
  names(res)=c("up","dn")
  if(order & !oneTable){
    if(nrow(res$up)>0) res$up=res$up[order(-res$up$area),]
    if(nrow(res$dn)>0) res$dn=res$dn[order(-res$dn$area),]
  }
  if(oneTable){
    res=rbind(res$up,res$dn)
    if(order & nrow(res)>0) res=res[order(-res$area),]
  }
  return(res)
}



# loess smoothing of statistic
loessSmoothStat <- function(stat,pns,pos,WBP=600) {
	require(limma)
	stat.smooth <- stat
	Indexes=split(seq(along=pns),pns)
	for(i in seq(along=Indexes)){
		if(i%%10000 == 0) cat(".")
		Index=Indexes[[i]]
		lens=median(diff(pos[Index]))
		spans=WBP/lens/length(Index) ##this is actually span for loess
		if(length(Index) > 3 & lens > 0) {
		    stat.smooth[Index] <- loessFit(stat[Index],pos[Index], span = spans)$fitted
		} else stat.smooth[Index] = median(stat[Index])
	}
	return(stat.smooth)
}

##you can pass cutoff through the ...
# doesn't take abs(area)
regionFinderPos <-function(x,regionNames,chr,position,y=x,
                       summary=mean,ind=seq(along=x),order=TRUE,
                       ...){

  tmp = regionFinder(x=x,regionNames=regionNames,
		chr=chr,position=position,y=x,
                       summary=mean,ind=seq(along=x),order=TRUE,oneTable=FALSE,
                       ...)
  
	out = tmp[["up"]]
	return(out)

}

# via Jeffrey Leek
# pvalue calculator
edge.pvalue <- function(stat, stat0, pool=TRUE) {
  err.func <- "edge.pvalue"
  m <- length(stat)
  if(pool==TRUE) {
    if(is.matrix(stat0)) {stat0 <- as.vector(stat0)}
    m0 <- length(stat0) 
    v <- c(rep(T, m), rep(F, m0))
    v <- v[order(c(stat,stat0), decreasing = TRUE)]
    u <- 1:length(v)
    w <- 1:m
    p <- (u[v==TRUE]-w)/m0
    p <- p[rank(-stat)]
    p <- pmax(p,1/m0)
  } else {
    if(is.vector(stat0)) {
      err.msg(err.func,"stat0 must be a matrix.")
      return(invisible(1))
    }
    if(ncol(stat0)==m) {stat0 <- t(stat0)}
    if(nrow(stat0)!=m){
      err.msg(err.func,"Number of rows of stat0 must equal length of stat.")
      return(invisible(1))
    }
    stat0 <- (stat0 - matrix(rep(stat,ncol(stat0)),byrow=FALSE,nrow=m)) >= 0
    p <- apply(stat0,1,mean)
    p <- pmax(p,1/ncol(stat0))
  }
  return(p)
}

# via Jeffrey Leek
# qvalue calculator that doesn't use tcltk
edge.qvalue <- function(p,lambda = seq(0, 0.9, 0.05), pi0.method = "smoother",
                         fdr.level = NULL, robust = FALSE,smooth.df = 3, smooth.log.pi0 = FALSE, ...) {

  err.func <- "edge.qvalue"
  if (min(p) < 0 || max(p) > 1) {
    err.msg(err.func,"P-values not in valid range.")
    return(invisible(1))
  }
  if (length(lambda) > 1 && length(lambda) < 4) {
    err.msg(err.func,"If length of lambda greater than 1, you need at least 4 values.")
    return(invisible(1))
  }
  if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >= 1)) {
    err.msg(err.func,"Lambda must be in [0,1).")
    return(invisible(1))
  }
  m <- length(p)
  if (length(lambda) == 1) {
    if (lambda < 0 || lambda >= 1) {
      err.msg(err.func,"Lambda must be in [0,1).")
      return(invisible(1))
    }
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0 <- min(pi0, 1)
  } else {
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
    }
    if (pi0.method == "smoother") {
      if (smooth.log.pi0){ 
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0 <- predict(spi0, x = max(lambda))$y
      }
      if (smooth.log.pi0) {
        pi0 <- exp(pi0)
      }
      pi0 <- min(pi0, 1)
    }
    else if (pi0.method == "bootstrap") {
      minpi0 <- min(pi0)
      mse <- rep(0, length(lambda))
      pi0.boot <- rep(0, length(lambda))
      for (i in 1:100) {
        p.boot <- sample(p, size = m, replace = TRUE)
        for (i in 1:length(lambda)) {
          pi0.boot[i] <- mean(p.boot > lambda[i])/(1 - lambda[i])
        }
        mse <- mse + (pi0.boot - minpi0)^2
      }
      pi0 <- min(pi0[mse == min(mse)])
      pi0 <- min(pi0, 1)
    }
    else {
      err.msg(err.func,"'pi0.method' must be one of 'smoother' or 'bootstrap'")
      return(invisible(1))
    }
  }
  if (pi0 <= 0) {
    err.msg(err.func,"The estimated pi0 <= 0. Check that you have valid\np-values or use another lambda method.")
    return(invisible(1))
  }
  if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
    err.msg(err.func,"'fdr.level' must be within (0,1].")
    return(invisible(1))
  }
  u <- order(p)
  qvalue.rank <- function(x) {
    idx <- sort.list(x)
    fc <- factor(x)
    nl <- length(levels(fc))
    bin <- as.integer(fc)
    tbl <- tabulate(bin)
    cs <- cumsum(tbl)
    tbl <- rep(cs, tbl)
    tbl[idx] <- tbl
    return(tbl)
  }
  v <- qvalue.rank(p)
  qvalue <- pi0 * m * p/v
  if (robust) {
    qvalue <- pi0 * m * p/(v * (1 - (1 - p)^m))
  }
  qvalue[u[m]] <- min(qvalue[u[m]], 1)
  for (i in (m - 1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
  }
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, pvalues = p, fdr.level = fdr.level, significant = (qvalue <= fdr.level), lambda = lambda)
  }
  else {
    retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, pvalues = p, lambda = lambda)
  }
  class(retval) <- "qvalue"
  return(retval)
}


# arima bootstrapping/sampling, loess smoothing
arimaBootSmoothStat <- function(stat,n_boot,pns,pos,WBP=300) {
	null.dat <- matrix(nr = length(stat), nc = n_boot)

	# lstat = log(stat)
	Indexes=split(seq(along=pns),pns)

	
	# AR process
	phi = matrix(NA, ncol = 4, nrow = length(Indexes))
	cat("calculating AR coefficient.","\n")
	for(i in seq(along = Indexes)) {
		if(i%%10000 == 0) cat(".")
		Index = Indexes[[i]]
		if(length(Index) > 100) {
			phi[i,] = ar(stat[Index],order.max=4, aic = F)$ar
		}
	}
	
	# ar model parameters
	ar1 = mean(phi[,1], na.rm = T)
	ar.sd = sd(stat, na.rm=T) #added na.rm=T
	
	cat("\nAR Bootstrapping.","\n")
	for(j in 1:n_boot) {
		cat(".")
		null.dat[,j] = arima.sim(list(ar = ar1), n = length(stat), 
			sd = ar.sd)
	}
	

	cat("\nSmoothing of null data.","\n")
	
	require(limma)
	for(i in seq(along=Indexes)){
		if(i%%10000 == 0) cat(".")
		Index=Indexes[[i]]
		lens=median(diff(pos[Index]))
		spans=WBP/lens/length(Index)
		for(j in 1:n_boot) {
			if(length(Index) > 3 & lens > 0) {
				null.dat[Index,j] <- loessFit(null.dat[Index,j],pos[Index], 
					span = spans)$fitted
			}
		}
	}
	return(null.dat)
}

getPckg = function(pckg) install.packages(pckg, repos = "http://cran.r-project.org")
splitit <- function(x) split(seq(along=x),x)
	
# via Rafael Irizarry
library(RColorBrewer)
mypar <- function(a=1,b=1,brewer.n=8,brewer.name="Dark2",...){
 par(mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
 par(mfrow=c(a,b),...)
 palette(brewer.pal(brewer.n,brewer.name))
}