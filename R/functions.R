

##setMethod(".calls_ICE", c("integer", "HmmOptions"),
.splitAnnotation <- function(object, ...){  ##, P.CHOM.Normal, P.CHOM.LOH, SAMPLE=1){
	pkgs <- strsplit(annotation(object), ",")[[1]]	
	fn <- featureNames(object)
	if(length(pkgs) > 1){
		require(RSQLite) || stop("RSQLite package not available")
		sql <- "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
		object1 <- object2 <- object
		annotation(object1) <- pkgs[1]
		annotation(object2) <- pkgs[2]
		tmp1 <- dbGetQuery(db(object1), sql)
		tmp2 <- dbGetQuery(db(object2), sql)

		idx.pkg1 <- match(tmp1[["man_fsetid"]], fn)
		idx.pkg1 <- idx.pkg1[!is.na(idx.pkg1)]

		idx.pkg2 <- match(tmp2[["man_fsetid"]], fn)
		idx.pkg2 <- idx.pkg2[!is.na(idx.pkg2)]

		featureData(object)$platform <- rep(NA, nrow(object))
		featureData(object)$platform[idx.pkg1] <- pkgs[1]
		featureData(object)$platform[idx.pkg2] <- pkgs[2]
	}
	object
}
	
	
findBreaks <- function(x, states, position, chromosome, sample,
		       lik1, lik2){
	if(length(chromosome)==1) chromosome <- rep(chromosome, length(position))
	splitby <- factor(cumsum(c(1, diff(x) != 0)))
	indices <- split(1:length(x), splitby)
	len <- sapply(indices, length)
	S <- states[sapply(split(x, splitby), unique)]
	pos <- t(sapply(split(position, splitby), range))
	size <- apply(t(sapply(split(position, splitby), range)), 1, diff)
	chr <- sapply(split(chromosome, splitby), unique)
	breaks <- data.frame(matrix(NA, length(chr), 7))
	colnames(breaks) <- c("sample", "chr", "start", "end", "nbases", "nprobes", "state")
	breaks$sample <- rep(sample, length(chr))
	breaks$chr <- chr
	breaks$start <- pos[, 1]
	breaks$end <- pos[, 2]
	breaks$nbases <- size
	breaks$nprobes <- len
	breaks$state <- S
	if(!missing(lik1) & !missing(lik2)){
		likdiff <- function(index, lik1, lik2, state){
			state <- unique(state[index])
			i <- range(index)
			if(min(i) > 1) i[1] <- i[1]-1
			if(max(x) < nrow(lik1)) i[2] <- i[2]+1
			##the more positive the better
			d1 <- diff(lik1[i, state])
			d2 <- diff(lik2[i, "N"])
			LR <- d1-d2
			return(LR)
		}
		LR <- as.numeric(sapply(indices, likdiff, lik1=lik1, lik2=lik2, state=x))
	}
	breaks <- breaks[sapply(chr, length) == 1, ]
	breaks$chr <- unlist(breaks$chr)
	return(breaks)
}





getChromosomeArm <- function(snpset){
	data(chromosomeAnnotation, package="SNPchip", envir=environment())
	chrAnn <- as.matrix(chromosomeAnnotation)
	tmp <- position(snpset) <= chromosomeAnnotation[chromosome(snpset), "centromereStart"]
	chromosomeArm <- as.character(tmp)
	chromosomeArm[chromosomeArm == "TRUE"] <- "p"
	chromosomeArm[chromosomeArm == "FALSE"] <- "q"
	chromosomeArm
}


addFeatureData <- function(snpset){
	featureNames <- featureNames(snpset)

	if(!("arm" %in% fvarLabels(snpset))){
		arm <- getChromosomeArm(snpset)
	}
	fD <- data.frame(cbind(chromosome(snpset), arm=chromosomeArm), row.names=featureNames,
			 stringsAsFactors=FALSE)
	fD$position <- position(snpset) 
	colnames(fD) <- c("chromosome", "arm", "position")
	varMetadata <- data.frame(labelDescription=c("chromosome", "chromosomal arm", "physical position"), row.names=colnames(fD))
	featureData <- new("AnnotatedDataFrame", data=fD, varMetadata=varMetadata)
	featureData
}
	
##.calculateYlim <- function(object, op){
##  if("copyNumber" %in% ls(assayData(object))){
##    ##only print this if there is more than one sample or more than 1 chromosome to plot
##    if(length(unique(chromosome(object))) > 1 || ncol(object) > 1){
##      print("one.ylim is FALSE. Calculating ylim based on the percentiles of the copy number distribution")
##    }
##    if(op$log == "y"){
##      ylim <- range(copyNumber(object), na.rm=TRUE)
##    } else{
##      ##use 1 and 99 quantiles to avoid outliers
##      ylim <- c(quantile(copyNumber(object), prob=0.001, na.rm=TRUE),
##                quantile(copyNumber(object), prob=0.999, na.rm=TRUE))
##    }
##  } else{
##    ##ylimits for genotypes??
##    y <- .getY(object)
##    ##jitter the genotype calls
##    y <- jitter(y, amount=0.05)
##    ylim <- range(y)
##  }
##  ylim
##}

##For each of the breaks in xBreak, find the cytoband(s)
.locateCytoband <- function(breaks, cytoband){
  cytoband <- cytoband[cytoband[, "chrom"] == breaks["chromosome"], ]
  ##include cytobands that begin before the break and end in the break or after the break
  bands <- cytoband[cytoband[, "chromEnd"] > as.numeric(breaks["start"]) & cytoband[, "chromStart"] < as.numeric(breaks["last"]), "name"]
  bands <- paste(bands, collapse=", ")
}

##setMethod("calculateBreakpoints", "HmmParameter",
##          function(object, x, position, chromosome, states, digits=2, sampleNames, ...){
.calculatebreaks <- function(object, digits=2){
	if(class(object) != "HmmPredict") stop("arg object in .calcualtebreaks should be an object of class HmmPredict")
	position <- position(object)
	chromosome <- unique(chromosome(object))
	if(length(chromosome) > 1) stop("only 1 chromosome at a time")
	states <- states(object)
	sampleNames <- sampleNames(object)
	if(length(sampleNames) > 1) stop("only 1 sample at a time")
	x <- predictions(object)
	if(any(is.na(x))){
		Nmissing <- sum(is.na(x))
	} else Nmissing <- 0
	N <- length(x)
	colnames <- c("id", "chr", "state", "size", 
		      "N", "start", "last", "prev", "next")
	
	if(length(unique(x)) == 1){
		breaks <- data.frame(sampleNames,
				     chromosome,
				     states[unique(x)],
				     (max(position) - min(position))/1e6,
				     length(position),
				     min(position)/1e6,
				     max(position)/1e6,
				     NA,
				     NA)
		colnames(breaks) <- colnames
		return(breaks)
	}
	
	d <- diff(x)
	if(sum(abs(d), na.rm=TRUE) < 1){
		breaks <- data.frame(sampleNames,
				     chromosome,
				     states[unique(x)],
				     (max(position) - min(position))/1e6,
				     length(position),
				     min(position)/1e6,
				     max(position)/1e6,
				     NA, NA)
		colnames(breaks) <- colnames
		return(breaks)		    
	}
	index <- c(1, (2:N)[d != 0 & !is.na(d)])
	index <- cbind(index[-length(index)], index[2:length(index)])
	index[, 2] <- index[, 2] - 1
	lastBreak <- index[nrow(index), 2]
	index <- rbind(index, c(lastBreak+1, N))  

	##index is a matrix of indices in object where breaks occured.
	##Replace by physical position.
	physical.positions <- matrix(position[as.vector(index)], nrow(index), ncol(index))
	physical.positions <- physical.positions/1e6
	colnames(physical.positions) <- c("start", "last")
	size <- (physical.positions[, "last"] - physical.positions[, "start"])
	N <- index[, 2] - index[, 1] + 1
	physical.positions <- cbind(physical.positions, size, N)
	colnames(physical.positions)[3:4] <- c("size", "N")
	predictedState <- function(i, x){
		if(any(is.na(i))){
			stop("missing values")
		}
		pred <- x[(i[1]:i[2])]
		pred <- pred[!is.na(pred)]
		pred <- unique(pred)
		if(length(pred) > 1) {
			stop("predictions not unique")
		}
		pred
	}
	##The vector of predicted states (integer)
	ps <- t(apply(index, 1, predictedState, x))
	states <- states[as.numeric(ps)]
	breaks <- data.frame(physical.positions)
	breaks$state <- states
	breaks$chrom <- chromosome
	##breaks$id <- rep(sampleNames(object), nrow(breaks))
	
	##---------------------------------------------------------------------------
	##Find positions of adjacent SNPs
	##---------------------------------------------------------------------------	
	adjacentSnps <- function(x, position, chromosome, digits){
		start <- as.numeric(x["start"])*1e6
		last <- as.numeric(x["last"])*1e6
		chrom <- x["chrom"]
		if(any(position < start & chromosome == chrom)){
			prevSnp <- max(position[position < start & chromosome == chrom])
		} else prevSnp <- NA
		if(any(position > last & chromosome == chrom)){                 
			nextSnp <- min(position[position > last & chromosome == chrom])
		} else nextSnp <- NA
		if(!is.na(prevSnp)) prevSnp <- prevSnp/1e6
		if(!is.na(nextSnp)) nextSnp <- nextSnp/1e6
		adj <- c(prevSnp, nextSnp)
		adj
	}
	adjacent <- t(apply(breaks, 1, adjacentSnps, position=position,
			    chromosome=chromosome, digits=digits))
	colnames(adjacent) <- c("prev", "next")
	breaks <- cbind(breaks, adjacent)
	breaks$id <- rep(sampleNames, nrow(breaks))
	colnames(breaks) <- c("start", "last", "size", "N", "state", "chr", "prev", "next", "id")
	##Add columns for the name of the first SNP and the name of the last SNP
	breaks <- breaks[, colnames]
	breaks
}
	

calculateCnSE <- function(object,
                          referenceSet,
                          epsilon=0.1){
	if(min(copyNumber(object), na.rm=TRUE) > 0){
		print("Robust estimates of the standard error are on the log2 scale.")
		print("Transforming copy number in object to log2 scale")
		copyNumber(object) <- log2(copyNumber(object))
	} else{
		warning("Negative values in the copy number.  Assume that data in the copyNumber element of assayData has been suitably transformed and is approximately Gaussian")
	}
	if(!missing(referenceSet)){
		if(min(copyNumber(referenceSet), na.rm=TRUE) > 0){
			print("Transforming copy number in referenceSet to log2 scale")
			copyNumber(referenceSet) <- log2(copyNumber(referenceSet))
		} else{

			warning("Negative values in the copy number.  Calculations assume that data in the copyNumber element of assayData has been suitably transformed and is approximately Gaussian")
			
		}
		
	}
	if(missing(referenceSet)){
		print("using object to compute across sample standard deviations...assumes reasonable sample size")
		referenceSet <- object
	}
	require(genefilter) || stop("genefilter not available")
	is.autosome <- chromosome(object) %in% as.character(1:22) 
	object <- object[is.autosome, ]  
	referenceSet <- referenceSet[match(featureNames(object), featureNames(referenceSet)), ]
	robustSD <- function(X) diff(quantile(X, probs=c(0.16, (1-0.16)), na.rm=TRUE))/2 
	within.sd <- diag(apply(copyNumber(object), 2, robustSD))
	across.sd <- apply(copyNumber(referenceSet), 1, robustSD)
	##  across.sd <- rowSds(copyNumber(referenceSet), na.rm=TRUE)
	across.sd <- matrix(across.sd, nrow=nrow(object), ncol=ncol(object), byrow=FALSE)
	##scale across.sd by the median sd of the sample
	median.across.sd <- median(across.sd, na.rm=TRUE)
	std.across.sd <- across.sd/median.across.sd
	SE <- std.across.sd %*% within.sd
	SE[SE == 0] <- epsilon
	rownames(SE) <- featureNames(object)
	SE
}

scaleTransitionProbability <- function(states, SCALE=2, normalLabel="N"){##HmmOptions
	matrix(1, length(states), length(states))
}

viterbi <- function(initialStateProbs,
		    emission,
		    tau,
		    arm,
		    tau.scale,
		    verbose=TRUE,
		    returnLikelihood=FALSE){
	if(class(emission) == "data.frame")
		emission <- as.matrix(emission)
	S <- ncol(emission)
	T <- nrow(emission)
	if(length(initialStateProbs) != S){
		stop("initialStateProbs (the initial state probabilities, should be a numeric vector of length S, where S is the number of hidden states")
	}		
	if(any(is.na(emission))){
		##if(verbose) message("Converting missing values in the emission matrix to 0")
		emission[is.na(emission)] <- 0
	}
	if(any(is.nan(emission))){
		if(verbose) message("some of the log emission probabilities are NaN.  Replacing with 0") 		
		emission[is.nan(emission)] <- 0
	}
	if(any(is.infinite(emission))){
		if(verbose) message("some of the log emission probabilities are infinite.  Replacing with 0's") 		
		emission[is.infinite(emission)] <- 0
	}
	if(missing(arm)) arm <- rep(as.integer(0), T)
	if(length(arm) != T) {
		if(verbose) message("arm not the right length.  assuming all values on same chromosomal arm")
		arm <- rep(as.integer(0), T)
	}
	if(missing(tau.scale)){
		##if(verbose) message("tau.scale not specified.  not scaling the genomic distance")
		tau.scale <- matrix(1, S, S)
	}
	if(any(!(dim(tau.scale) == S))){
		if(verbose) message("tau.scale not the right dimension.  not scaling the genomic distance")
		tau.scale <- matrix(1, S, S)
	}
	result <- vector("integer", T)
	delta <- matrix(as.double(0), nrow=T, ncol=S)
	tmp <- list(as.matrix(as.double(as.matrix(emission))),
		    as.double(as.matrix(initialStateProbs)),
		    as.matrix(as.double(tau)),
		    as.character(arm),
		    as.matrix(as.double(tau.scale)),
		    as.integer(S),
		    as.integer(T),
		    result,
		    as.matrix(as.double(delta)))
	fit <- .C("viterbi",
		  tmp[[1]],
		  tmp[[2]],
		  tmp[[3]],
		  tmp[[4]],
		  tmp[[5]],
		  tmp[[6]],
		  tmp[[7]],
		  tmp[[8]],
		  tmp[[9]])
	if(!returnLikelihood){
		return(fit[[8]])
	} else{
		tmp <- list(states=fit[[8]], likelihood=fit[[9]])
		return(tmp)
	}
}

	
##Code this in C
##viterbi <- function(beta, pi, tau, arm, tau.scale, normalIndex, distance=NULL, constant=1e8){
##	##beta: log emission probabilities
##	##pi: log initial state probabilities (vector of length S)
##	##tau: transition probabilities (original scale)
##	##tau.scale: matrix on original scale
##	S <- ncol(beta)
##	T <- nrow(beta)  ##matrix of T rows and S columns
##	delta <- psi <- matrix(NA, T, S)
##	delta[1, ] <- pi + beta[1, ]
##	psi[1, ] <- rep(0, S)
##	i <- which(rowSums(is.na(beta)) > 0)
##	beta[i, ] <- 0
##	if(!is.null(distance)){
##		tau <- exp(-2*distance/1e8)
##		tauNormal <- exp(-2*distance/constant)		
##	}
##	for(t in 2:T){
##		if(t %% 10000 == 0) cat(".")
##		if(arm[t] != arm[t-1]){
##			delta[t, ] <- pi + beta[t, ]
##			psi[t, ] <- rep(0, S)
##			next()
##		}
##		AA <- matrix(tau[t-1], nr=S, nc=S)  ##
##		epsilon <- (1-tau[t-1])/(S-1)
##		eNormal <- (1-tauNormal[t-1])/(S-1)		
##		AA[upper.tri(AA)] <- AA[lower.tri(AA)] <- epsilon
##		AA[normalIndex, ] <- rep(eNormal, S)
##		AA[normalIndex, normalIndex] <- tauNormal[t-1]
##			
##		AA <- log(AA*tau.scale)  
##		for(j in 1:S){
##			tmp <- delta[t-1, ] + AA[, j]
##			delta[t, j] <- max(tmp) + beta[t, j]
##			psi[t, j] <- order(tmp, decreasing=TRUE)[1]
##		}
##	}
##	Pstar <- max(delta[nrow(delta), ])
##	qhat <- rep(NA, nrow(delta))
##	qhat[T] <- order(delta[T, ], decreasing=TRUE)[1]
##	for(t in (T-1):1){
##		if(arm[t] != arm[t+1]){
##			qhat[t] <- order(delta[t, ], decreasing=TRUE)[1]
##		} else {
##			qhat[t] <- psi[t+1, qhat[t+1]]
##		}
##	}
##	return(qhat)
##}


##How often is a SNP altered across samples
##makeTable <- function(object, state, by, unit=1000, digits=3){
##	stop("function needs work")
##	if(missing(state)) stop("must specify state")
##	if(missing(by)) by <- "chr"
##	breaks <- list()
##	for(j in 1:ncol(object)){
##		breaks[[j]] <- findBreaks(predictions(object)[, j],
##					  chromosome=chromosome(object),
##					  position=position(object),
##					  states=states(object),
##					  sample=sampleNames(object)[j])
##		size <- breaks[[j]][, "end"]-breaks[[j]][, "start"]
##		breaks[[j]]$size <- 
##	}
##	breaks <- do.call("rbind", breaks)
##	X <- breaks[breaks[, "state"] == state, ]
##	##split by chromosome
##	X <- split(X, X[, by])
##	##Frequency of alterations of a specific type per chromosome
##	if(length(X) == 1){
##		freq <- nrow(X[[1]])
##	} else{
##		freq <- sapply(X, nrow)
##	}
##	##Length of alterations
##	med.L <- sapply(X, function(x) median(as.numeric(x[, "size"])*unit, na.rm=TRUE))
##	avg.L <- sapply(X, function(x) mean(as.numeric(x[, "size"])*unit, na.rm=TRUE))
##	sd.L <- sapply(X, function(x) sd(as.numeric(x[, "size"])*unit, na.rm=TRUE))
##	##Number of SNPs per alteration
##	med.N <- sapply(X, function(x) median(as.numeric(x[, "N"]), na.rm=TRUE))
##	avg.N <- sapply(X, function(x) mean(as.numeric(x[, "N"]), na.rm=TRUE))
##	sd.N <- sapply(X, function(x) sd(as.numeric(x[, "N"]), na.rm=TRUE))	
##	stats <- as.data.frame(cbind(freq, med.L, avg.L, sd.L, med.N, avg.N, sd.N))
##	colnames(stats) <- c("Freq", "med(length)", "avg(length)", 
##			     "sd(length)", "med(n.snp)", "avg(n.snp)", 
##			     "sd(n.snp)")
##	stats[, c(3, 4, 6, 7)] <- round(stats[, c(3, 4, 6, 7)], digits)
##	if("X" %in% rownames(stats)){
##		rownames(stats)["X"] <- 23
##	}
##	stats <- stats[order(as.numeric(rownames(stats))), ]
##	
##	##Marginal stats
##	X <- breaks[breaks[, "state"] == state, ]
##	if(nrow(stats) <= 1){
##		return(stats)
##	}
##	freq <- c(median(stats[, "Freq"]),
##		  mean(stats[, "Freq"]),
##		  sd(stats[, "Freq"]))
##	med.L <- c(median(stats[, "med(length)"]),
##		   mean(stats[, "med(length)"]),
##		   sd(stats[, "med(length)"]))
##	avg.L <- c(median(stats[, "avg(length)"]),
##		   mean(stats[, "avg(length)"]),
##		   sd(stats[, "avg(length)"]))
##	avg.L <- round(avg.L, digits)
##	sd.L <- c(median(stats[, "sd(length)"]),
##		   mean(stats[, "sd(length)"]),
##		   sd(stats[, "sd(length)"]))
##	sd.L <- round(sd.L, digits)
##	med.N <- c(median(stats[, "med(n.snp)"]),
##		   mean(stats[, "med(n.snp)"]),
##		   sd(stats[, "med(n.snp)"]))		
##	avg.N <- c(median(stats[, "avg(n.snp)"]),
##		   mean(stats[, "avg(n.snp)"]),
##		   sd(stats[, "avg(n.snp)"]))
##	avg.N <- round(avg.N, digits)
##	sd.N <- c(median(stats[, "sd(n.snp)"]),
##		   mean(stats[, "sd(n.snp)"]),
##		   sd(stats[, "sd(n.snp)"]))
##	sd.N <- round(sd.N, digits)
##	overall <- cbind(freq, med.L, avg.L, sd.L, med.N, avg.N, sd.N)
##	rownames(overall) <- c("median", "avg", "sd")
##	stats <- rbind(stats, overall)
##	stats
##}

scaleTransitionToState <- function(cols, AA, SCALE){
	epsilon <- AA[, cols]
	epsilon.tilde <- epsilon/SCALE
	AA[, cols] <- epsilon.tilde
	d <- rowSums(epsilon-epsilon.tilde)
	diag(AA) <- diag(AA)+d
	AA
}


.copynumberEmission <- function(copynumber,
				states,
				mu,
				uncertainty,
				takeLog,
				verbose=TRUE){
	cne <- copynumber
	location <- mu
	if(missing(takeLog)) stop("must specify whether to take the log2 of the copy number matrix")
	if(missing(uncertainty)) stop("must supply uncertainty estimates")	
	if(takeLog){
		if(verbose) message("Calculating emission probabilities of log(copy number)")
		i <- which(!is.na(as.vector(cne)))
		##Take log of hidden state
		if(any(location < 0)){
			if(verbose) message("negative values in copy number will be NA's")
		}
		tmp <- rep(NA, length(as.vector(cne)))
		tmp[i] <- log2(cne[i])
		cne <- matrix(tmp, nrow=nrow(snpset), ncol=ncol(snpset))		
		location <- log2(location)
	} else{
		if(verbose) message("no transformation of copy number")
	}
	fn <- rownames(cne)
	S <- length(states)
	cne <- array(cne, dim=c(nrow(cne), ncol(cne), S))
	dimnames(cne) <- list(rownames(cne), colnames(cne), states)
	##assume true copy number mean is the same for all samples
	##Easiest to keep everything the same dimension, even if not needed at this point
	location <- aperm(array(location, dim=c(S, ncol(cne), nrow(cne))))
	
	##Use SNP-specific standard errors
	i <- which(!(is.na(as.vector(1/uncertainty))) & !is.na(as.vector(cne)))
	if(all(1/uncertainty > 0, na.rm=TRUE)){
		if(verbose) message("Using uncertainty as standard deviation for the copy number")
		scale <- array(uncertainty, dim=dim(cne))
	} else{
		stop("confidence scores in slot cnConfidence must be positive")
	}
	k <- which(!is.na(as.vector(cne)))
	if(!identical(dim(cne), dim(location))) stop("dimensions must be the same")
	emission.cn <- rep(NA, length(as.vector(location)))
	emission.cn[k] <- dnorm(as.vector(cne)[k], as.vector(location)[k], as.vector(scale)[k])
	emission.cn <- array(emission.cn, dim=dim(cne))
	dimnames(emission.cn) <- list(rownames(cne), colnames(cne), states)
	if(verbose) message("returning the emission probability on the log scale")
	return(log(emission.cn))
}

.genotype.emission <- function(calls, states, ICE=FALSE, probHomCall, pkgs){
##	if(!("calls" %in% assayDataElementNames(object@snpset))){
##		warning("calls not an element of the assayData for object@snpset.")
##		return()
##	}
##	if(!validObject(object)){
##		stop("HmmOptions object is not valid")
##	}
##	snpset <- object@snpset
##	states <- object@states
	S <- length(states)	
##	gte <- array(calls(snpset), c(nrow(snpset), ncol(snpset), S))
	gte <- array(calls, c(nrow(calls), ncol(calls), S))
##	dimnames(gte) <- list(featureNames(snpset), sampleNames(snpset), states)
	dimnames(gte) <- list(rownames(calls), colnames(calls), states)
##	if(!object@calls.ICE) {
	if(!ICE){
		##for algorithms that produce bi-allelic calls, only 3
		##calls are possible.  The emission probability for
		##the call is the emission probability for a
		##homozygous call versus the emission probability for
		##a heterozygous call
		i <- which(as.vector(gte) == 1 | as.vector(gte) == 3)
		j <- which(as.vector(gte) == 2)
		k <- which(as.vector(gte) == 4 | is.na(as.vector(gte)))
		##prob. for each state.
		##cn.location <- aperm(array(object@copyNumber.location, dim=c(S, ncol(snpset), nrow(snpset))))
		##probHom <-  aperm(array(object@probHomCall, dim=c(S, ncol(snpset), nrow(snpset))))
		probHom <- aperm(array(probHomCall, dim=c(S, ncol(calls), nrow(calls))))
		emission.gt <- vector("numeric", length=length(as.vector(gte)))
		emission.gt[i] <- as.vector(probHom)[i]
		emission.gt[j] <- (1-as.vector(probHom))[j]
		emission.gt[k] <- NA
		emission.gt[c(i, j)] <- log(emission.gt[c(i, j)])
		emission.gt <- array(emission.gt, dim=dim(gte))
	}  else {
		##pkgs <- strsplit(annotation(object@snpset), ",")[[1]]
##		if(length(pkgs) > 1){
##			object@snpset <- .splitAnnotation(snpset)
##		} else{
##			featureData(object@snpset)$platform <- NA
##			featureData(object@snpset)$platform <- annotation(snpset)
##		}
		##emission.gt <- array(NA, dim=c(nrow(object@snpset), ncol(object@snpset), 2))
		emission.gt <- array(NA, dim=c(nrow(calls), ncol(calls), 2))
		dimnames(emission.gt) <- list(rownames(calls), colnames(calls), c("LOH", "N"))
		for(i in 1:length(pkgs)){
			j <- which(pkgs[[i]])
			##j <- which(fData(object@snpset)$platform == pkgs[i])
			##annotation(object@snpset) <- pkgs[1]
			stop("need to fix this function")
			emission.gt[j, , ] <- .getCallEmission(object[j, ])
		}
		emission.gt <- log(emission.gt)
	}
	return(emission.gt)##log scale
}










