breaks <- function(x, states, position, chromosome, chromosomeAnnotation=NULL, verbose=FALSE){
	x <- as.matrix(x)
	if(missing(states)) stop("must specify labels for the hidden states")
	if(missing(position)) stop("must provide physical position of loci")
	if(missing(chromosome)) stop("must indicate chromosome")
	results <- vector("list", ncol(x))
	sns <- colnames(x)
	for(j in 1:ncol(x)){
		results[[j]] <- findBreaks(x[, j], states=states,
					   position=position,
					   chromosome=chromosome,
					   sample=colnames(x)[j],
					   chromosomeAnnotation=chromosomeAnnotation,
					   verbose=verbose)
		results[[j]]$sample <- sns[j]
	}
	if(length(results) == 1){
		results <- results[[1]]
	} else  {
		results <- do.call("rbind", results)
	}
	return(results)
}


findBreaks <- function(x, states, position, chromosome, sample, chromosomeAnnotation=NULL, verbose=FALSE){
	if(is.null(chromosomeAnnotation)){
		data(chromosomeAnnotation, package="SNPchip", envir=environment())
		chrAnn <- as.matrix(chromosomeAnnotation)
	}			
	chrAnn <- chromosomeAnnotation
	if(is.matrix(x)) if(ncol(x) > 1) stop("x should be a vector or matrix with 1 column")	
	if(!is.integer(chromosome)) {
		chromosome <- chromosome2integer(chromosome)
	}
	if(!all(chromosome %in% 1:24)){
			warning("Chromosome annotation is currently available for chromosomes 1-22, X and Y")
			message("Please add/modify data(chromosomeAnnotation, package='SNPchip') to accomodate special chromosomes")
			stop()
	}
	if(!is.integer(position)) {
		if(verbose) message("Coerced position to an integer.")
		position <- as.integer(position)
	}	
	##ensure that the reported breaks do not span the centromere
	chromosome <- integer2chromosome(chromosome)
	uchrom <- unique(chromosome)
	positionList <- split(position, chromosome)
	positionList <- positionList[match(uchrom, names(positionList))]
	arm <- list()
	for(i in seq(along=uchrom)){
		arm[[i]] <- as.integer(ifelse(positionList[[i]] <= chrAnn[uchrom[i], "centromereStart"], 0, 1))	
	}		
	arm <- unlist(arm)
	if(length(chromosome)==1) chromosome <- rep(chromosome, length(position))
	splitby <- factor(cumsum(c(1, diff(x) != 0 | diff(arm) != 0)))
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
	breaks$nbases <- size+1
	breaks$nprobes <- len
	breaks$state <- S
##	if(!missing(lik1) & !missing(lik2)){
##		likdiff <- function(index, lik1, lik2, state){
##			state <- unique(state[index])
##			i <- range(index)
##			if(min(i) > 1) i[1] <- i[1]-1
##			if(max(x) < nrow(lik1)) i[2] <- i[2]+1
##			##the more positive the better
##			d1 <- diff(lik1[i, state])
##			d2 <- diff(lik2[i, "N"])
##			LR <- d1-d2
##			return(LR)
##		}
##		LR <- as.numeric(sapply(indices, likdiff, lik1=lik1, lik2=lik2, state=x))
##	}
	breaks <- breaks[sapply(chr, length) == 1, ]
	breaks$chr <- unlist(breaks$chr)
	return(breaks)
}


genotypeEmission <- function(genotypes, conf, states, probHomCall,
			     probMissing, verbose=TRUE){
	##function to replace .getCallEmission
	##message("Use genotypeEmissionCrlmm if using crlmm-processed data on affy250k or affy 6.0")
	if(!is.numeric(genotypes)) stop("genotypes must be integers (1=AA, 2=AB, 3=BB, 4=missing")
	emission <- array(genotypes, dim=c(nrow(genotypes), ncol(genotypes), length(states)))
	for(s in seq(along=states)){
		tmp <- genotypes
		tmp[tmp == 1 | tmp == 3] <- probHomCall[s]
		tmp[tmp == 2] <- 1-probHomCall[s]
		if(!missing(probMissing)) tmp[tmp == 4 | is.na(tmp)] <- probMissing[s]
		emission[, , s] <- tmp
	}
	dimnames(emission) <- list(rownames(genotypes),
				   colnames(genotypes),
				   states)
	logemit <- log(emission)
	return(logemit)
}

genotypeEmissionCrlmm <- function(genotypes, conf,
				  pHetCalledHom=0.001,
				  pHetCalledHet=0.995,
				  pHomInNormal=0.8,
				  pHomInRoh=0.999,  ##Pr(AA or BB | region of homozygosity)
				  annotation){
	if(missing(annotation)) stop("must specify annotation")
	##data(list=paste(annotation, "Conf", sep=""), package="VanillaICE", envir=environment())
	if(length(pHomInNormal) == nrow(conf)){  ##convert to vector
		pHomInNormal <- as.numeric(matrix(pHomInNormal, nrow(conf), ncol(conf), byrow=FALSE))
	} else pHomInNormal <-  rep(pHomInNormal, length(as.integer(genotypes)))
	loader(paste(annotation, "Conf.rda", sep=""), .vanillaIcePkgEnv, "VanillaICE")
	hapmapP <- getVarInEnv("reference")
	hapmapP[, 2] <- 1-exp(-hapmapP[, 2]/1000)
	##p = 1-exp(-X/1000)
	##1000*log(1-p)=X
	confidence <- 1-exp(-conf/1000)	
	i11 <- hapmapP[, 1] == 3  ##called homozygous truth homozygous
	i12 <- hapmapP[, 1] == 4  ##called homozygous truth heterozygous
	i21 <- hapmapP[, 1] == 1  ##called het truth hom
	i22 <- hapmapP[, 1] == 2  ##called het truth het
	f11 <- density(hapmapP[i11, 2], from=0, to=1, n=1e3)
	f12 <- density(hapmapP[i12, 2], from=0, to=1, n=1e3)
	f21 <- density(hapmapP[i21, 2], from=0, to=1, n=1e3)
	f22 <- density(hapmapP[i22, 2], from=0, to=1, n=1e3)

	##-------------------------------------------------------------------------
	##distribution of observed call probabilities when the call is homozygous
	##-------------------------------------------------------------------------

	##-------------------------------------------------------------------------
	##P(phat | LOH, gte)  
	##-------------------------------------------------------------------------	
	GT <- as.integer(genotypes)
	confidence <- as.numeric(confidence)
	pTruthIsNormal <- pTruthIsRoh <- rep(NA, length(GT))	
	confidence[confidence==0] <- 0.01 ##Otherwise, NA's result
	hom <- which(GT == 1 | GT == 3)
	observedPcalledHom <- cut(confidence[hom], breaks=f11$x, labels=FALSE)
	pTruthIsRoh[hom] <- f11$y[observedPcalledHom]
	het <- which(GT == 2)
	observedPcalledHet <- cut(confidence[het], breaks=f11$x, labels=FALSE)
	pTruthIsRoh[het] <- f21$y[observedPcalledHet]

	##-------------------------------------------------------------------------
	##Calculate P(phat | Normal, HOM)
	##-------------------------------------------------------------------------					
	chet1 <- f22$y[cut(confidence[het], breaks=f22$x, labels=FALSE)] 
	chet2 <- f21$y[cut(confidence[het], breaks=f21$x, labels=FALSE)]
	##term5[1]=P(true genotype is HET | genotype call is AB, state is normal)
	pTruthIsNormal[het] <- chet1*pHetCalledHet + chet2*(1-pHetCalledHet)
	##chom1=called homozygous truth heterozygous
	chom1 <- f12$y[cut(confidence[hom], breaks=f12$x, labels=FALSE)]
	##chom2=called homozygous truth homozygous
	chom2 <- f11$y[cut(confidence[hom], breaks=f11$x, labels=FALSE)]
	##chom4 <- 0.9999    ##P(HOM|CHOM)
	##probability that the true state is HOM when genotype call is homozygous
	##pHetCalledHom = P(true genotype is HET | calls is AA or BB, state is normal)
	pTruthIsNormal[hom] <- chom1*pHetCalledHom + chom2*(1-pHetCalledHom)
	fNormal <- fLoh <- rep(NA, length(GT))
	fNormal[hom] <- pHomInNormal[hom] * pTruthIsNormal[hom]
	fNormal[het] <- (1-pHomInNormal[het]) * pTruthIsNormal[het]
	fLoh[hom] <- pHomInRoh * pTruthIsRoh[hom]
	fLoh[het] <- (1-pHomInRoh) * pTruthIsRoh[het]

	f <- array(NA, dim=c(nrow(genotypes), ncol(genotypes), 2))
	dimnames(f) <- list(rownames(genotypes),
			    colnames(genotypes),
			    c("norm", "ROH"))	
	f[, , "norm"] <- matrix(fNormal, nrow(genotypes), ncol(genotypes))
	f[, , "ROH"] <- matrix(fLoh, nrow(genotypes), ncol(genotypes))
	f[f  == 0] <- min(f[f > 0], na.rm=TRUE)
	f <- log(f)

	return(f)	
}

copynumberEmission <- function(copynumber,
			       states,
			       mu,
			       sds,
			       takeLog,
			       verbose=TRUE, na.rm=TRUE){
	##if(!is.matrix(mu)) mu <- matrix(mu, nrow(copynumber), length(states), byrow=TRUE)
	if(missing(takeLog)) stop("must specify whether to take the log2 of the copy number matrix")
	fn <- rownames(copynumber)
	S <- length(states)
	if(takeLog){
		if(any(copynumber == 0, na.rm=TRUE)){
			message("0 copy number estimates are replaced with 0.05 to avoid nan's when taking the log")
			copynumber[copynumber == 0] <- 0.05
		}
	}
	if(any(rowSums(is.na(copynumber)) == ncol(copynumber))){
		stop("Some rows have all missing values. Exclude these before continuing.")
	}
	if(any(colSums(is.na(copynumber)) == nrow(copynumber))){
		stop("Some samples have all missing values. Exclude these samples before continuing.")
	}
	if(!missing(sds)){
		if(!is.matrix(sds)) sds <- matrix(sds, nrow(copynumber), ncol(copynumber), byrow=TRUE)
	} else {
		sds <- robustSds(copynumber, takeLog=takeLog, na.rm=na.rm)
	}
	tmp <- rowSums(!is.finite(sds))
	if(any(sds == 0, na.rm=TRUE)){
		warning("some sds were zero.  Replacing with typical value")
		sds[sds == 0] <- median(sds, na.rm=TRUE)
	}
	if(takeLog){
		copynumber <- log2(copynumber)
		mu <- log2(mu)
	}
	emission.cn <- array(NA, dim=c(nrow(copynumber), ncol(copynumber), S))
	if(!is.matrix(mu))
		mu <- matrix(mu, nrow(copynumber), length(mu), byrow=TRUE)	
	for(j in 1:ncol(copynumber)){
		cn <- matrix(copynumber[, j], nrow(copynumber), ncol(mu))
		sd <- matrix(sds[, j], nrow(copynumber), ncol(mu))
		k <- which(!is.na(as.numeric(cn)))
		##emission.cn <- rep(NA, length(as.vector(mu)))
		tmp <- rep(NA, length(as.numeric(mu)))
		tmp[k] <- dnorm(as.numeric(cn)[k], as.numeric(mu)[k], as.numeric(sd)[k])
		emission.cn[, j, ] <- tmp
		##emission.cn[k] <- dnorm(as.vector(copynumber)[k], as.vector(mu)[k], as.vector(scale)[k])
		##emission.cn <- array(emission.cn, dim=dim(copynumber))
		##dimnames(emission.cn) <- list(rownames(copynumber), colnames(copynumber), states)
		##if(verbose) message("returning the emission probability on the log scale")
	}
	dimnames(emission.cn) <- list(rownames(copynumber),
				      colnames(copynumber),
				      states)
	return(log(emission.cn))
}

rowMAD <- function(x, y, ...){
	notna <- !is.na(x)
	sds <- 1.4826*rowMedians(abs(x-rowMedians(x, ...)), ...)
	return(sds)
}

	


robustSds <- function(x, takeLog=FALSE, ...){
	if(!is.matrix(x)) stop("the copy number estimates should be provided as a matrix")
	if(takeLog) x <- log2(x)
	if(ncol(x) > 3){
		sds1 <- rowMAD(x, ...)
		sds1 <- matrix(sds1, nrow(x), ncol(x))
		sds2 <- apply(x, 2, "mad", ...)
		sds2 <- sds2/median(sds2, na.rm=TRUE)
		sds <- t(t(sds1) * sds2)
	} else {
		sds <- apply(x, 2, "mad", ...)
		sds <- matrix(sds, nrow(x), ncol(x), byrow=TRUE)
	}
	return(sds)
}


					

viterbi <- function(emission,
		    tau,
		    arm,
		    initialStateProbs,		    
		    verbose=FALSE,
		    chromosome,
		    position,
		    sampleNames,
		    locusNames,
		    normalIndex,
		    returnLikelihood=FALSE,
		    normal2altered=1,
		    altered2normal=1,
		    altered2altered=1){
	if(class(emission) != "array") stop("emission probabilities must be an array: snps, samples, states. ")
	if(missing(normalIndex)) stop("Must specify integer for normalIndex")
	if(!is.numeric(normalIndex)) stop("normalIndex should be numeric")
	if(normal2altered <= 0) stop("normal2altered must be > 0")
	if(altered2normal <= 0) stop("altered2normal must be > 0")
	if(altered2altered <= 0) stop("altered2altered must be > 0")	
	if(returnLikelihood){
		if(missing(normalIndex)) stop("Must specify index for 'normal' state")
	}
	if(!is.array(emission)){
		stop("the emission probabilities should be an array")
	}
	results <- matrix(NA, nrow(emission), ncol(emission))
	S <- dim(emission)[3]
	T <- nrow(emission)
	if(missing(initialStateProbs)){
		initialStateProbs <- log(rep(1/S, S))
	}
	if(length(initialStateProbs) != S){
		stop("initialStateProbs (the initial state probabilities, should be a numeric vector of length S, where S is the number of hidden states")
	}
	if(!all(initialStateProbs <= 0)){
		if(all(initialStateProbs >= 0 & initialStateProbs <= 1)){
			initialStateProbs <- log(initialStateProbs)
		} else stop("initial state probabilities should be a probability or a log probability")
	}
	if(any(is.na(emission))){
		if(verbose) message("Converting missing values in the emission matrix to 0")
		emission[is.na(emission)] <- 0
	}
	if(any(is.nan(emission))){
		message("some of the log emission probabilities are NaN.  Replacing with 0") 		
		emission[is.nan(emission)] <- 0
	}
	if(any(emission < -50)){
		message("some of the log emission probabilities are very small -- probable outliers.  Replacing with a small value (-10)") 		
		emission[emission < -50] <- -50
	}
	if(missing(arm)){
		message("chromosome arm not specified...HMM is not fit separately to each chromosomal arm")
		arm <- rep(as.integer(0), T)
	}
	if(length(arm) != T) {
		message("arm not the right length.  assuming all values on same chromosomal arm")
		arm <- rep(as.integer(0), T)
	}
	if(missing(tau)){
		stop("transition probabilities not specified")
	}
	if(length(tau) != T) stop("tau must have length T")
	## The last index is arbitrary and, by default, is NA. Must replace this by a number-- C can not handle NA's
	tau[is.na(tau)] <- 0
	delta <- matrix(as.double(0), nrow=T, ncol=S)
	lik <- array(NA, dim=c(T, ncol(results), S))
	for(j in 1:ncol(results)){
		result <- rep(as.integer(0), T)		
		tmp <- list(as.matrix(as.double(as.matrix(emission[, j, ]))),##beta
			    as.double(as.matrix(initialStateProbs)),##initialP
			    as.matrix(as.double(tau)),##tau
			    as.integer(arm),##arm
			    as.integer(S),##number of states
			    as.integer(T),##number of loci
			    result,##placeholder for results
			    as.matrix(as.double(delta)),##delta
			    normal2altered,##c1
			    altered2normal,##c2
			    altered2altered,##c3
			    as.integer(normalIndex),
			    as.double(rep(0, S^2)))##normalState
		##inputs <- tmp
		##.GlobalEnv[["inputs"]] <- inputs
		tmp2 <- .C("viterbi",		
			   tmp[[1]],
			   tmp[[2]],
			   tmp[[3]],
			   tmp[[4]],
			   tmp[[5]],
			   tmp[[6]],
			   tmp[[7]],
			   tmp[[8]],
			   tmp[[9]],
			   tmp[[10]],
			   tmp[[11]],
			   tmp[[12]],
			   tmp[[13]])  ##can verify that the transition prob. matrix is correct (for last snp)
		##check transition probs.
		M <- matrix(tmp2[[13]], S, S)
		if(!all(is.finite(M))) stop("Infinite values in transition prob. matrix")
		if(!all.equal(rowSums(exp(M)), rep(1, S))){
			warning("Rows of the transition probability matrix do not sum to 1")
		}
		results[, j] <- tmp2[[7]]
		lik[, j, ] <- matrix(tmp2[[8]], nrow(results))
	}
	if(returnLikelihood){
		lrdiff <- vector("list", ncol(results))
		for(j in 1:ncol(results)){
			indices <- c(0, cumsum(diff(results[, j]) != 0 | diff(arm) != 0))
			startLoc <- which(c(0, diff(arm)) != 0)
			stopLoc <- c(startLoc -1, nrow(results))
			lr <- rep(NA, length(table(indices)))
			for(k in seq(along=table(indices))){
				index <- range(which(indices == (k-1)))
				thisState <- unique(results[index, j])
				##if(thisState == normalIndex) next()
				tmp <- lik[index, j, thisState] - lik[index, j, normalIndex]
				lr[k] <- sum(tmp)
			}
			lrdiff[[j]] <- lr
		}
		results <- list(results, lrdiff)
	}
	if(!is.null(rownames(emission)))
		rownames(results) <- rownames(emission)
	if(!is.null(colnames(emission)))	
		colnames(results) <- colnames(emission)
	return(results)		
}

transitionProbability <- function(chromosome, position, TAUP=1e8, chromosomeAnnotation, verbose=FALSE){
	if(!is.integer(chromosome)) {
		chromosome <- chromosome2integer(chromosome)
	}
	if(!all(chromosome %in% 1:24)){
			warning("Chromosome annotation is currently available for chromosomes 1-22, X and Y")
			message("Please add/modify data(chromosomeAnnotation, package='SNPchip') to accomodate special chromosomes")
			stop()
	}
	if(!is.integer(position)) {
		if(verbose) message("Coerced position to an integer.")
		position <- as.integer(position)
	}
	if(length(chromosome) != length(position)) stop("chromosome and position arguments must be the same length")
	if(missing(chromosomeAnnotation)){
		if(verbose) message("chromosomeAnnotation not specified... using centromere locations from SNPchip")
		data(chromosomeAnnotation, package="SNPchip", envir=environment())
		chromosomeAnnotation <- as.matrix(chromosomeAnnotation)		
	}
	chrAnn <- chromosomeAnnotation
	uchrom <- unique(integer2chromosome(chromosome))
	chromosomeArm <- tau <- vector("list", length(uchrom))
	positionList <- split(position, chromosome)
	positionList <- positionList[match(unique(chromosome), names(positionList))]
	for(i in seq(along=uchrom)){
		chromosomeArm[[i]] <- as.integer(ifelse(positionList[[i]] <= chrAnn[uchrom[i], "centromereEnd"], 0, 1))
		##probability SNP is informative.  Note, the last index is assigned zero -- this is arbitrary.  
		tau[[i]] <- c(exp(-2 * diff(positionList[[i]])/TAUP), NA)
		if(any(tau[[i]] > 1, na.rm=TRUE)){
			##values greater than one occur when the diff is negative.
			stop("check that the physical position is ordered from smallest to largest")
		}
		if(any(tau[[i]] < 0, na.rm=TRUE)) stop("some of the computed transition probabilities less than zero.")		
	}
	chromosomeArm <- unlist(chromosomeArm)
	chromosomeArm <- cumsum(c(0, diff(chromosomeArm) != 0 | diff(chromosome) != 0))
	tau <- unlist(tau)
	annotation <- cbind(chromosome, position, chromosomeArm, tau)
	colnames(annotation) <- c("chromosome", "position", "arm", "transitionPr")
	range.tau <- range(annotation[, "transitionPr"], na.rm=TRUE)
	##check range
	if(range.tau[1] < 0 | range.tau[2] > 1) stop("Transition probabilities are not valid.  Check that the snpset object has been ordered by chromosome and physical position.")
	##assign a minimum pr. that the snp is informative
	tau[tau < 0.5] <- 0.5
	annotation[, "transitionPr"] <- tau
	return(annotation)
}

	
##Code this in C
viterbiR <- function(emission, initialP, tau, arm){
##	emission: log emission probabilities
##	initialP: log initial state probabilities (vector of length S)
##	tau: transition probabilities (original scale)
##	tau.scale: matrix on original scale
	S <- ncol(emission)
	T <- nrow(emission)  ##matrix of T rows and S columns
	delta <- psi <- matrix(NA, T, S)
	delta[1, ] <- initialP + emission[1, ]
	psi[1, ] <- rep(0, S)
	i <- which(rowSums(is.na(emission)) > 0)
	tau.scale <- 1
	#emission[i, ] <- 0
	for(t in 2:T){
		if(t %% 10000 == 0) cat(".")
		if(arm[t] != arm[t-1]){
			delta[t, ] <- initialP + emission[t, ]
			psi[t, ] <- 0
			next()
		}
##		AA <- matrix(NA, nr=S, nc=S)
##		for(i in 1:S){
		AA <- matrix(tau[t-1], nr=S, nc=S)  
		epsilon <- (1-tau[t-1])/(S-1)
		##eNormal <- (1-tauNormal[t-1])/(S-1)		
		AA[upper.tri(AA)] <- AA[lower.tri(AA)] <- epsilon
		##AA[normalIndex, ] <- rep(eNormal, S)
		##AA[normalIndex, normalIndex] <- tauNormal[t-1]
		AA <- log(AA*tau.scale)  
		for(j in 1:S){
			tmp <- delta[t-1, ] + AA[, j]
			delta[t, j] <- max(tmp) + emission[t, j]
			psi[t, j] <- order(tmp, decreasing=TRUE)[1]
		}
	}
	Pstar <- max(delta[nrow(delta), ])
	qhat <- rep(NA, nrow(delta))
	qhat[T] <- order(delta[T, ], decreasing=TRUE)[1]
	for(t in (T-1):1){
		if(arm[t] != arm[t+1]){
			qhat[t] <- order(delta[t, ], decreasing=TRUE)[1]
		} else {
			qhat[t] <- psi[t+1, qhat[t+1]]
		}
	}
	return(qhat)
}


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

hmm <- function(object,
		states,
		mu=NULL,
		probs=NULL,
		takeLog=FALSE,
		initialP,
		returnSegments=TRUE,
		TAUP=1e8,
		verbose=FALSE,
		ice=FALSE,
		envir,
		normalIndex){
	if(missing(envir)) envir <- new.env()	
	if(!all(c("position", "chromosome") %in% fvarLabels(object))){
		stop("'position' and 'chromosome' must be in fvarLabels(object), or transitionPr must be provided")
	}
	object <- object[order(chromosome(object), position(object)), ]
	envir[["locusset"]] <- object
	if(missing(states)) stop("must specify states")
	envir[["states"]] <- states
	if(missing(initialP)) initialP <- rep(1, length(states))/length(states)
	envir[["initialP"]] <- initialP
	envir[["mu"]] <- mu
	if(is.null(probs)){
		if(ice){
			probs <- c(0.05, 0.99, 0.7, 0.999)
		}
	}
	envir[["probs"]] <- probs
	envir[["takeLog"]] <- takeLog
	envir[["returnSegments"]] <- returnSegments
	envir[["TAUP"]] <- TAUP
	envir[["verbose"]] <- verbose
	envir[["ice"]] <- ice
	tau <- transitionProbability(chromosome=chromosome(object),
				     position=position(object),
				     TAUP=TAUP,
				     verbose=verbose)
	arm <- tau[, "arm"]
	transitionPr <- tau[, "transitionPr"]
	envir[["transitionPr"]] <- transitionPr
	envir[["arm"]] <- arm
	if(takeLog){
		copyNumber(object) <- log2(copyNumber(object))
		mu <- log2(mu)
	}		
	if(verbose) message("Calculating emission probabilities")
	calculateEmission(object=object,
			  mu=mu,
			  probs=probs,
			  envir=envir,
			  states=states,
			  verbose=verbose,
			  ice=ice)
	if(!is.null(envir[["emission.gt"]])){
		emission <- envir[["emission.cn"]]+envir[["emission.gt"]]
	} else {
		emission <- envir[["emission.cn"]]
	}
	envir[["emission"]] <- emission
	##emission <- .GlobalEnv[["emission.cn"]]+.GlobalEnv[["emission.gt"]]
	viterbiResults <- viterbi(initialStateProbs=log(initialP),
				  emission=emission,
				  tau=transitionPr,
				  arm=arm,
				  verbose=verbose,
				  normalIndex=normalIndex)
	dimnames(viterbiResults) <- list(featureNames(object), sampleNames(object))	
	if(returnSegments){
		viterbiResults <- breaks(x=viterbiResults,
					 states=states,
					 position=position(object),
					 chromosome=chromosome(object),	
					 verbose=verbose)
	}
	return(viterbiResults)
}

setMethod("update", "environment", function(object, ...){
	if(length(ls(object)) == 0) stop("nothing to update")
	hmm(object=object[["locusset"]],
	    states=object[["states"]],
	    mu=object[["mu"]],
	    probs=object[["probs"]],
	    takeLog=object[["takeLog"]],
	    initialP=object[["initialP"]],
	    returnSegments=object[["returnSegments"]],
	    TAUP=object[["TAUP"]],
	    verbose=object[["verbose"]],
	    ice=object[["ice"]],
	    envir=object)	
})
	












