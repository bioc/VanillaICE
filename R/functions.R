rowMAD <- function(x, y, ...){
	notna <- !is.na(x)
	sds <- 1.4826*rowMedians(abs(x-rowMedians(x, ...)), ...)
	return(sds)
}

robustSds <- function(x, takeLog=FALSE, ...){
	if(!is.matrix(x)) stop("x is not a matrix")
	if(takeLog) x <- log2(x)
	sds <- apply(x, 2, "mad", ...)
	sds <- matrix(sds, nrow(x), ncol(x), byrow=TRUE)
	dimnames(sds) <- dimnames(x)
	return(sds)
}

stackRangedData <- function(object){
	if(is(object, "list")){
		if(length(object)==1){
			object <- object[[1]]
		} else {
			object <- RangedDataList(object)
			object <- stack(object)
			ix <- match("sample", colnames(object))
			if(length(ix) > 0) object <- object[, -ix]
		}
		if(is(object, "RangedDataHMM")) return(object)
	}
	rangedData <- RangedDataHMM(ranges=ranges(object),
				    chromosome=object$chrom,
				    sampleId=object$sampleId,
				    state=object$state,
				    coverage=object$coverage,
				    LLR=object$LLR)
	return(rangedData)
}

hmm.setup <- function(...) .Deprecated("hmm.setup function is deprecated.  hmm will run directly on a class of oligoSnpSet, SnpSet, CopyNumberSet, or BeadStudioSet. See ?hmm for details.")

tnorm <- function(x, mean, sd, lower=0, upper=1){
       phi <- function(x, mu, sigma) dnorm(x, mu, sigma)
       Phi <- function(x, mu, sigma) pnorm(x, mu, sigma)
       res <- phi(x, mean, sd)/(Phi(upper, mean, sd)-Phi(lower, mean, sd))
       res
}

trnorm <- function(x){
	function(mu, sigma, upper=1, lower=0){
		dnorm(x, mu, sigma)/(pnorm(upper, mu, sigma) - pnorm(lower, mu, sigma))
	}
}

getColClasses <- function(filename, lrr.colname, baf.colname){
	tmp <- read.table(filename,
			  row.names=NULL,
			  header=TRUE,
			  stringsAsFactors=FALSE,
			  sep="\t", nrows=50)
	j <- grep(lrr.colname, colnames(tmp))
	k <- grep(baf.colname, colnames(tmp))
	if(length(j) == 0 || length(k)==0)
		stop("lrr.colname or baf.colname not in header")
	colClasses <- as.character(sapply(tmp[1, ], class))
	index <- setdiff(seq_along(colClasses), c(1, j, k))
	colClasses[index] <- rep("NULL", length(index)) ## don't read in the other columns
	colClasses
}

## copied from MinimumDistance
read.bsfiles <- function(path="", filenames, ext="", row.names=1,
			 sep="\t",
			 lrr.colname="Log.R.Ratio",
			 baf.colname="B.Allele",
			 drop=FALSE,
			 colClasses,
			 nrows=1.8e6,
			 ...){
	if(path != ""){
		fnames <- file.path(path, paste(filenames, ext, sep=""))
	} else fnames <- paste(filenames, ext, sep="")
	stopifnot(all(file.exists(fnames)))
	if(missing(colClasses)){
		colClasses <- getColClasses(fnames[1], lrr.colname, baf.colname)
	}
	for(i in seq_along(fnames)){
		##cat(".")
		tmp <- read.table(fnames[i],
				  row.names=row.names,
				  sep=sep,
				  nrows=nrows,
				  header=TRUE,
				  stringsAsFactors=FALSE,
				  colClasses=colClasses,
				  check.names=FALSE,
				  comment.char="", ...)
		tmp <- as.matrix(tmp)
		if(i==1){
			dat <- array(NA, dim=c(nrow(tmp), 2, length(fnames)))
			j <- grep(lrr.colname, colnames(tmp))
			k <- grep(baf.colname, colnames(tmp))
			stopifnot(length(j)==1)
			stopifnot(length(k)==1)
			if(!drop){
				dimnames(dat) <- list(rownames(tmp),
						      c("lrr", "baf"),
						      basename(filenames))
			}
			##lrr.data <- matrix(NA, nrow(tmp), length(filenames))
			##baf.data <- matrix(NA, nrow(tmp), length(filenames))
		}
		dat[, 1, i] <- tmp[, j]
		dat[, 2, i] <- tmp[, k]
	}
	##cat("\n")
	return(dat)
}

getProbB <- function(cdfname, featurenames){
	cdfpath <- system.file("extdata", package=cdfname)
	if(file.exists(file.path(cdfpath, "pb_gw6.rda"))){
		load(file.path(cdfpath, "pb_gw6.rda"))
		probB <- rep(NA, length(featurenames))
		pb <- pb[names(pb) %in% featurenames]
		index <- match(names(pb), featurenames)
		stopifnot(!any(is.na(index)))
		probB[index] <- pb
	} else probB <- rep(0.5, length(featurenames))
	return(probB)
}

keyOffFirstFile <- function(filename, cdfname, genome, lrr.colname, baf.colname, ...){
	## read in one file
	## return feature matrix in chromosome, position order
	dat <- read.bsfiles(filenames=filename, lrr.colname=lrr.colname, baf.colname=baf.colname, ...)
	cdfpath <- system.file("extdata", package=cdfname)
	fnames <- list.files(cdfpath)
	mult.build <- length(grep("_hg1[89].rda$", fnames))>=1
	if(mult.build){
		load(file.path(cdfpath, paste("snpProbes_", genome, ".rda", sep="")))
		load(file.path(cdfpath, paste("cnProbes_", genome, ".rda", sep="")))
	} else {
		if(genome=="hg18") stop("Currently only genome build hg19 is available in this annotation package")
		load(file.path(cdfpath, "snpProbes.rda"))
		load(file.path(cdfpath, "cnProbes.rda"))
	}
	snpProbes <- get("snpProbes")
	cnProbes <- get("cnProbes")
	features <- rbind(snpProbes, cnProbes)
	keep.index <- which(rownames(features) %in% rownames(dat))
	features <- features[keep.index, ]

	index.order <- order(features[, "chrom"], features[, "position"])
	features <- features[index.order, ]

	issnp <- as.logical(rownames(features) %in% rownames(snpProbes))
	probB <- as.integer(getProbB(cdfname, rownames(features))*100)
	arm <- oligoClasses:::.getArm(features[, "chrom"], features[, "position"], genome)
	index <- match(rownames(features), rownames(dat))

	identical(rownames(dat)[index], rownames(features))
	##features2 <- cbind(features, issnp, probB, arm, index)
	features2 <- data.frame(chrom=features[, "chrom"],
				position=features[, "position"],
				probB=probB,
				isSnp=issnp,
				arm=arm,
				index=index)
	##colnames(features2) <- c(colnames(features), "isSnp", "probB", "arm", "index")
	return(features2)
}

## ad-hoc.  Do we want to put priors on the means?
constrainMu <- function(mu, is.log){
	if(is.log){
		is.ratio <- all.equal(mu[3], 0, tolerance=0.2) == TRUE
		if(is.ratio){
			mu[3] <- ifelse(mu[3] < -0.2, -0.2, mu[3])
			mu[3] <- ifelse(mu[3] > 0.2, 0.2, mu[3])
			mu[1] <- ifelse(mu[1] > -1, -1, mu[1])
			mu[2] <- ifelse(mu[2] > -0.25, -0.25, mu[2])
			mu[4] <- ifelse(mu[4] < 0.25, 0.25, mu[4])
			mu[4] <- ifelse(mu[4] > 0.6, 0.6, mu[4])
			mu[5] <- ifelse(mu[5] < 0.65, 0.65, mu[5])
		} else {
			mu[3] <- ifelse(mu[3] < 0.75, 0.75, mu[3])
			mu[3] <- ifelse(mu[3] > 1.25, 1.25, mu[3])
			mu[1] <- ifelse(mu[1] > -0.5, -0.5, mu[1])
			mu[2] <- ifelse(mu[2] > 0.7, 0.7, mu[2])
			mu[4] <- ifelse(mu[4] < 1.25, 1.25, mu[4])
			mu[4] <- ifelse(mu[4] > 2, 2, mu[4])
			mu[5] <- ifelse(mu[5] < 2, 2, mu[5])
		}
	} else {
		mu[1] <- ifelse(mu[1] > 0.8, 0.8, mu[1])
		mu[2] <- ifelse(mu[2] > 1.7, 1.7, mu[2])
		mu[2] <- ifelse(mu[2] < 0.9, 0.9, mu[2])
		mu[3] <- ifelse(mu[3] > 2.2, 2.2, mu[3])
		mu[3] <- ifelse(mu[3] < 1.8, 1.8, mu[3])
		mu[4] <- ifelse(mu[4] < 2.3, 2.3, mu[4])
		mu[4] <- ifelse(mu[4] > 3, 3, mu[4])
		mu[5] <- ifelse(mu[5] < 3, 3, mu[5])
	}
	return(mu)
}

constrainSd <- function(sigma){
	sigma[1] <- min(3*sigma[3], max(sigma[1], sigma[3]))
	sigma[6] <- min(3*sigma[3], max(sigma[6], sigma[3]))
	sigma[5] <- min(3*sigma[3], max(sigma[5], sigma[3]))
	sigma[2:4] <- sigma[3]
	return(sigma)
}


computeTransitionProb <- function(x, TAUP, S, tauMAX=1-5e-6){
	p <- exp(-2*diff(x)/TAUP)
	minimum <- 1-1/((S-1)) + 0.01
	p <- pmax(p, minimum)
	p <- pmin(p, tauMAX)
	return(as.matrix(p))
}

copyNumberLimits <- function(is.log) if(is.log) return(c(-2.5, 3)) else return(c(0,10))

thresholdCopyNumber <- function(object, limits){
	object <- pmin(object, limits[2])
	object <- pmax(object, limits[1])
	object
}

centerCopyNumber <- function(object, is.snp){
	snp.index <- which(is.snp)
	mu.snp <- apply(object[snp.index, , drop=FALSE], 2, median, na.rm=TRUE)
	object[snp.index, ] <- sweep(object[snp.index, , drop=FALSE],  2, mu.snp)
	if(any(!is.snp)){
		np.index <- which(!is.snp)
		mu.np <- apply(object[np.index, , drop=FALSE], 2, median, na.rm=TRUE)
		object[np.index, ] <- sweep(object[np.index, , drop=FALSE], 2, mu.np)
	}
	object
}

guessCnScale <- function(x){
	is.log <- if(all(x >= 0, na.rm=TRUE)) FALSE else TRUE
	is.log
}

copyNumberStates <- function(normalCn) {
	is.ratio <- all.equal(normalCn, 1, tolerance=0.2) == TRUE
	if(is.ratio){
		return(c(0, 1/2, 1, 1, 3/2, 4/2))
	} else{
		is.absolute <- all.equal(normalCn, 2, tolerance=0.2)==TRUE
		if(is.absolute){
			return(c(0, 1, 2, 2, 3, 4))
		} else{
			stop("median copy number is not near 1 or 2")
		}
	  }
}

rescale <- function(x, l, u){
	b <- 1/(u-l)
	a <- l*b
	(x+a)/b
}

makeNonDecreasing <- function(x){
	d <- diff(x)
	if(all(d >= 0)) return(x)
	index <- which(d < 0)
	if(index[1] == 1){
		x[1] <- x[2]
		index <- index[-1]
		if(length(index) ==0)
			return(x)
	}
	l <- length(index)
	if(l == length(x)-1){
		i <- index[l]
		x[length(x)] <- x[length(x)-1]
		index <- index[-l]
		if(length(index) == 0) return(x)
	}
	x[index+1] <- x[index]
	return(x)
}

stackGRangesList <- function(fit, build){
	L <- length(fit[[1]])
	gr <- list()
	for(i in seq_len(L)){
		rdl <- lapply(fit, "[[", i)
		tmp <- stack(GRangesList(rdl))
		## above adds column called sample
		##ids <- unlist(sapply(fit, sampleNames))
		gr[[i]] <- tmp[, -match("sample", colnames(values(tmp)))]
		j <- match("sample.1", colnames(values(gr[[i]])))
		if(length(j) > 0){
			colnames(values(gr[[i]]))[j] <- "sample"
		}
	}
	sl <- getSequenceLengths(build)
	rd <- GRangesList(gr)#, seqlengths=sl)
	nms <- unique(as.character(runValue(unlist(seqnames(rd)))))
	seqlengths(rd) <- sl[match(nms, names(sl))]
	names(rd) <- names(fit[[1]])
	metadata(rd) <- list(genome=build)
	rd
}

generatorFun <- function(r, b, gt, is.snp, cnStates,
			 normalIndex, TAUP, limits, center,
			 prOutlierBAF, p.hom, position, is.log, computeLLR, chrom, verbose=FALSE){
	S <- length(cnStates)
	nc <- ncol(r)
	b <- b[is.snp, , drop=FALSE]
	nr <- nrow(r)
	nb <- nrow(b)
	np.index <- which(!is.snp)
	names(prOutlierBAF)[1] <- "prOutlier"
	CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")
	##
	## params for copy number
	sds <- apply(r, 2, mad, na.rm=TRUE)
	if(center) r <- centerCopyNumber(r, is.snp) + cnStates[normalIndex]
	initialCnParams <- function(j){
		mus <- cnStates
		sigmas <- rep(sds[j], S)
		## only works if the chromosome is not a 2-copy gain
		## iHets <- which(b[, j] > 0.45 & b[, j] < 0.55)
		## pHets <- length(iHets)/nr
		## if(pHets > 0.1) mus[c(3,4)] <- median(r[iHets], na.rm=TRUE)
		##mus[5] <- mus[3] + 2*sds[j]
		##mus[6] <- mus[3] + 3*sds[j]
		##p <- matrix(c(0.99,  0.01), S, 2, byrow=TRUE)
		p <- 0.01
		paramsCN <- list(mu=mus, sigmas=sigmas, p=p)
	}
	## params for BAF
	allele.prob <- getPrB()
	initialBafParams <- function(){
		musBAF <- c(0, 0.1, 1/3, 0.5, 2/3, 0.9, 1)
		sdsBAF <- c(0.02, rep(0.05, 5), 0.02)
		names(musBAF) <- names(sdsBAF) <- c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B")
		paramsBAF <- list(mus=musBAF, sigmas=sdsBAF, prOutlier=prOutlierBAF)
	}
	g2 <- g1 <- matrix(NA, nb, S)
	d1 <- matrix(NA, nb, S)
	isHom <- b < 0.02 | b > 0.98
	indexHom <- list()
	for(j in seq_len(nc)) indexHom[[j]] <- which(isHom[, j])
	emitr <- matrix(NA, nr, S)
	emitb <- matrix(NA, nb, S)
	cn.dunif <- dunif(r, limits[1], limits[2])
	Gtotal <- matrix(NA, nb, 4)
	updateBafEmission <- function(bf, params, j){
		mus <- params[["mus"]]
		sds <- params[["sigmas"]]
		##prOutlier <- params$prOutlier$initial
		prOutlier <- params$prOutlier$prOutlier
		computeEmitPr <- function(dlist, prOutlier, allele.prob){
			sumd <- 0
			q <- 1-prOutlier
			for(i in seq_along(dlist)){
				sumd <- sumd + dlist[[i]] * q*allele.prob[i]
			}
			##p <- matrix(prOutlier, nc, nr, byrow=TRUE)
			t(sumd+prOutlier)
		}
		trNormal <- trnorm(bf)
		##
		##---------------------------------------------------------------------------
		## copy number 1
		dA <- trNormal(mus["A"], sds["A"])
		dB <- trNormal(mus["B"], sds["B"])
		emitb[, 2] <- computeEmitPr(list(dA,dB), prOutlier, allele.prob[, c("A", "B")])
		emitb[, 4] <- emitb[, 2]
		##
		##
		##---------------------------------------------------------------------------
		## copy number 2
		dAB <- trNormal(mus["AB"], sds["AB"])
		emitb[, 3] <- computeEmitPr(list(dA, dAB, dB), prOutlier=prOutlier, allele.prob=allele.prob[, c("AA", "AB", "BB")])
		##
		##---------------------------------------------------------------------------
		## copy number 3
		dAAB <- trNormal(mus["AAB"], sds["AAB"])
		dABB <- trNormal(mus["ABB"], sds["ABB"])
		emitb[, 5] <- computeEmitPr(list(dA, dAAB, dABB, dB), prOutlier, allele.prob[, c("AAA", "AAB", "ABB", "BBB")])
		##
		##---------------------------------------------------------------------------
		## copy number 4
		dAAAB <- trNormal(mus["AAAB"], sds["AAAB"])
		dABBB <- trNormal(mus["ABBB"], sds["ABBB"])
		emitb[, 6] <- computeEmitPr(list(dA, dAAAB, dAB, dABBB, dB), prOutlier, allele.prob[, c("AAAA", "AAAB", "AABB", "ABBB", "BBBB")])
		##
		##
		## small p.hom makes homozygous genotypes less informative
		##
		if(p.hom < 1){
			i <- indexHom[[j]]
			if(length(i) > 0) emitb[i, c(3, 5, 6)] <- (1-p.hom)*emitb[i, 2] + p.hom*emitb[i, c(3, 5, 6)]
		}
		##if(length(np.index) > 0) emitb[np.index, ] <- 1
		nas <- anyMissing(emitb)
		if(nas) emitb[is.na(emitb)] <- 1
		##index <- which(rowSums(is.infinite(emitb)) > 0)
		return(emitb)
	}

	## might try dt here instead
	mydnorm <- function(x){
		function(mean, sd){
			mus <- matrix(mean, nr, S, byrow=TRUE)
			sds <- matrix(sd, nr, S, byrow=TRUE)
			dnorm(x, mus, sds)
		}
	}
	updateCnEmission <- function(cn, params, j){
		mus <- params[["mu"]]
		sds <- params[["sigmas"]]
		p <- params[["p"]]
		q <- 1-p
		normalDens <- mydnorm(cn)
		d <- normalDens(mus, sds)
		emitr <- q*d + p*cn.dunif[, j]
		emitr.na <- is.na(emitr)
		hasna <- any(emitr.na)
		if(hasna) emitr[emitr.na] <- 1
		return(emitr)
	}

	theoreticalMeans <- initialBafParams()$mus
	names(theoreticalMeans) <- c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B")
	##strictHets <- apply(b > 0.4 & b < 0.6, 2, mean, na.rm=TRUE) < 0.05
	##fewHets <- apply(b > 0.05 & b < 0.45 | b > 0.55 & b < 0.95, 2, mean, na.rm=TRUE) < 0.05
	updateBafParams <- function(bf, params, h, j) {
		weightedMoments <- function(w) {
			if(sum(w > 0.5, na.rm=TRUE) >= 10){
				w[w < 0.01 | is.na(w)] <- 0
				sum.w <- sum(w)
				mu <- sum(w*bf, na.rm=TRUE)/sum.w
				sigma <- sqrt((sum(w*(bf-mu)^2, na.rm=TRUE))/sum.w)
				return(c(mu, sigma))
			} else c(NA, NA)
		}
		mus <- params[["mus"]]
		if(verbose) print(mus)
		sigmas <- params[["sigmas"]]
		pOut <- params[["prOutlier"]]
		pOutlierMax <- pOut[["max"]]
		pOutlierMaxROH <- pOut[["maxROH"]]
		p <- pOut[["prOutlier"]]
		p <- pmin(p, pOutlierMax)
		p <- 1-p ## probability not an outlier
		m <- rowSums(h, na.rm=TRUE)
		h <- h/m
		calculateGamma <- function(d, h, allele.prob){
			sumd <- 0
			## conditional on not an outlier
			G <- length(d)
			## integrate over G possible genotypes
			## (marginal probability for CN given not an outlier)
			for(j in seq_len(G)){
				d[[j]] <- allele.prob[j]*d[[j]]
				sumd <- sumd+d[[j]]
			}
			## integrate out the outlier indicator to get the marginal probability for the copy number
			## (marginal probability for CN)
			## denom <- p*sumd + (1-p)
			lapply(d, function(d) d*h/sumd)
		}
		updateBMoments <- function(gammas, nms){
			moments <- t(sapply(gammas, weightedMoments))
			rownames(moments) <- nms
			moments
		}
		nms <- names(mus)
		trNormal <- trnorm(bf)
		##
		## For cn = 0, BAF~Unif(0,1)
		gamma0 <- dunif(bf, 0,1) * h[, 1]
		gammaTotal <- function(g) {
			sumg <- 0
			for(i in seq_len(length(g))) sumg <- sumg+g[[i]]
			sumg
		}
		## BAF means for cn 1
		dA <- trNormal(mus["A"], sigmas["A"])
		dB <- trNormal(mus["B"], sigmas["B"])
		gamma1 <- calculateGamma(list(dA, dB), h[, 2], allele.prob[, c("A", "B")])
		## these estimates will be unreliable if there is weak
		## evidence of hemizygous deletion as the moment estimatore includes
		## the forward/backward probabilities
		m1 <- updateBMoments(gamma1, nms=c("A", "B"))
		##
		## BAF means for cn 2
		dAB <- trNormal(mus["AB"], sigmas["AB"])
		gamma2 <- calculateGamma(list(dA, dAB, dB), h[, 3], allele.prob[, c("AA", "AB", "BB")])
		m2 <- updateBMoments(gamma2, nms=c("A", "AB", "B"))##["AB", , drop=FALSE]
		##
		## BAF means for cn 3
		dAAB <- trNormal(mus["AAB"], sigmas["AAB"])
		dABB <- trNormal(mus["ABB"], sigmas["ABB"])
		gamma3 <- calculateGamma(list(dA, dAAB, dABB, dB), h[, 5], allele.prob[, c("AAA", "AAB", "ABB", "BBB")])
		m3 <- updateBMoments(gamma3, nms=c("A", "AAB", "ABB", "B"))##[c("AAB", "ABB"), ]
		##
		## Update mixture components.
		dAAAB <- trNormal(mus["AAAB"], sigmas["AAAB"])
		dABBB <- trNormal(mus["ABBB"], sigmas["ABBB"])
		gamma4 <- calculateGamma(list(dA, dAAAB, dAB, dABBB, dB), h[, 6], allele.prob[,c("AAAA", "AAAB", "AABB", "ABBB", "BBBB")])
		m4 <- updateBMoments(gamma4, nms=c("A", "AAAB", "AB", "ABBB", "B"))##[c("AAAB", "ABBB"), ]
		##
		## use fv/bv to find which state is likely to give the best estimate of the means/variances for A, and B
		totalh <- apply(h, 2, sum, na.rm=TRUE)
		imax <- which.max(totalh)
		mA_B <- switch(paste("state", imax, sep=""),
			       state1=m1,#?? homozygous null
			       state2=m1,
			       state3=m2,
			       state4=m1,## use means of A and B for copy-neutral ROH
			       state5=m3,
			       state6=m4)
		Gtotal[, 1] <- gammaTotal(gamma1) ## A or B
		Gtotal[, 2] <- gammaTotal(gamma2) ## AB
		Gtotal[, 3] <- gammaTotal(gamma3) ## AAB, ABB
		Gtotal[, 4] <- gammaTotal(gamma4) ## AAAB, ABBB
		##mus.new <- c(mA_B[c("A", "B"),1], m2["AB", 1], m3[c("AAB", "ABB"), 1], m4[c("AAAB", "ABBB"), 1])
		mus.new <- c(0, 1, m2["AB", 1], m3[c("AAB", "ABB"), 1], m4[c("AAAB", "ABBB"), 1])
		names(mus.new)[1:3] <- c("A", "B", "AB")
		if(any(is.na(mus.new))){
			mus.new <- mus.new[match(names(theoreticalMeans),names(mus.new))]
			mus.new[is.na(mus.new)] <- theoreticalMeans[is.na(mus.new)]
		}
		mus.new <- makeMusBafNondecreasing(mus.new)
		if(mus.new["AB"] < 0.45 | mus.new["AB"] > 0.55) mus.new["AB"] <- params$mus[["AB"]]
		##mus.new <- constrainMuBaf(mus.new)
		## sd for hets. For hom null, hemizygous null, and regions of homozygosity, use the prior
		if(imax %in% c(3,5,6)){
			sAB <- switch(paste("state", imax, sep=""),
				       state3=m2["AB", 2],
				       state5=mean(c(m3["AAB", 2], m3["ABB", 2])),
				       state6=mean(c(m4["AAAB", 2], m4["AB", 2], m4["ABBB", 2])))
		} else sAB <- params$sigmas[["AB"]]
		## if there are few hets, the sd estimate for AB BAFs will be unreliable
		sigmas.new <- rep(sAB, 7)
		names(sigmas.new) <- names(params$sigmas)
		##sigmas.new <- c(mA_B[c("A", "B"),2], m2["AB", 2], m3[c("AAB", "ABB"), 2], m4[c("AAAB", "ABBB"), 2])
		sigmas.new[c("A", "B")] <- mA_B[c("A", "B"),2]
		if(any(is.na(sigmas.new))){
			sigmas.new <- sigmas.new[match(names(params$sigmas), names(sigmas.new))]
			sigmas.new[is.na(sigmas.new)] <- params$sigmas[is.na(sigmas.new)]
		}
		sigmas.new[c("A", "B")] <- pmax(sigmas.new[c("A", "B")], 0.001)
		##names(sigmas.new)[3] <- "AB"
		## do not allow any of the sd estimates to be less than the sd for the homozygous genotypes
		##
		## proportion outliers (how to account for homozygous deletions?)
		##    -- exclude consecutive outliers
		has.nas <- anyMissing(Gtotal)
		if(has.nas) {
			nas <- is.na(Gtotal[,1])
			Gtotal <- Gtotal[-which(nas), ]
		}
		isout <- rowMax(Gtotal) < 0.9
		p <- mean(diff(isout) != 0)
		params[["mus"]] <- mus.new
		params[["sigmas"]] <- sigmas.new
		pOut[["prOutlier"]] <- min(p, 0.001)
		params[["prOutlier"]] <- pOut
		return(params)
	}
	cn.sd.new <- rep(NA, S)
	updateCnParams <- function(cn, params, h, j) {
		mu <- params[["mu"]]
		if(verbose) print(mu)
		sigma <- params[["sigmas"]]
		p <- params[["p"]]; q <- 1-p
		min.sd <- sds[j]/2
		##h <- matrix(fv*bv, nr, S)
		##h <- matrix(h, nr, S)
		m <- rowSums(h, na.rm=TRUE)
		h <- h/m
		d.unif <- cn.dunif[, j]
		d <- q*dnorm(cn, mean=matrix(mu, nr, S, byrow=TRUE), sd=matrix(sigma, nr, S, byrow=TRUE)) + p*d.unif
		d.total <- d+d.unif
		d1 <- d/(d+d.unif)
		d2 <- 1-d1
		##pout <- mean(rowMax(d1) < 0.5) ## uniform distribution has higher probability than any of the cn states
		g1 <- h * d1  ## emission from observed states
		g1[g1 < 0.01] <- 0
		##g2 <- 1-g1
		##g2 <- h * (1-d1) ## emission from outlier state
		estimable <- apply(g1 > 0.5, 2, sum, na.rm=TRUE) >= 10
		if(!any(estimable)) {
			mu[c(3,4)] <- median(cn, na.rm=TRUE)
			mu[5] <- mu[3] + sds[j]
			mu[6] <- mu[3] + 2*sds[j]
			mu[2] <- mu[3] - 1*sds[j]
			mu[1] <- params$mu[1]
		}
		totalh <- apply(g1, 2, sum, na.rm=TRUE)
		##round(total.g1,3)
		## if total.g1 is small, then there's little data to estimate the means.
		mu.new <- colSums(g1*cn, na.rm=TRUE)/totalh
		if(any(!estimable))
			mu.new[!estimable] <- params[["mu"]][!estimable]
		imax <- which.max(totalh)
		if(imax == 4){
			mu.new[3] <- mu.new[4]
		} else mu.new[4] <- mu.new[3]
		##
		##
		## without this constraint, the means for the gain
		## states will creep towards normal copy number
		##
		##
		mu.new <- makeNonDecreasing(mu.new)
		one.sd <- params$sigmas[3]
		nsd <- 1
		diploid <- max(totalh[3], totalh[4])
		if(diploid > totalh[5]){
			mu.new[5] <- max(mu.new[3]+nsd*one.sd, mu.new[5])
			mu.new[6] <- max(mu.new[5]+nsd*one.sd, mu.new[6])
		} else{
			mu.new[3] <- min(mu.new[5]-nsd*one.sd, mu.new[3])
			mu.new[6] <- max(mu.new[5]+nsd*one.sd, mu.new[6])
		}
		mu.new[1] <- min(mu.new[2] - one.sd, mu.new[1])
		mu.new[4] <- mu.new[3]
		##if(fewHets[j]) mu.new[c(5,6)] <- params[["mu"]][c(5,6)]
		## For loop
		for(s in seq_len(S)) {
			cn.sd.new[s] <- sqrt(sum(g1[, s]*(cn-mu.new[s])^2, na.rm=TRUE)/totalh[s])
		}
		## assume sds are the same for non-null copy number states.
		## use the state that has the highest overall probability for estimating the sd.
		cn.sd.new[2:6] <- cn.sd.new[imax]
		cn.sd.new <- pmax(cn.sd.new, min.sd)
		if(any(!estimable))
			cn.sd.new[!estimable] <- params$sigmas[!estimable]
		params[["mu"]] <- mu.new
		params[["sigmas"]] <- cn.sd.new
		##total.g2 <- apply(g2, 2, sum, na.rm=TRUE)
		##denom.eq52 <- apply(g1 + g2, 2, sum, na.rm=TRUE)
		p <- mean(rowSums(g1, na.rm=TRUE) < 0.1)
		##p <- min(total.g2/denom.eq52)  ## altered states will have higher values of total.g2/denom.eq52.
		params[["p"]] <- min(p, 0.1)
		return(params)
	}

	transitionPr <- function(TAUP, tauMAX=1-5e-8){
		p <- exp(-2*diff(position)/TAUP)
		minimum <- 1-1/((S-1)) + 0.01
		p <- pmax(p, minimum)
		p <- pmin(p, tauMAX)
		return(as.matrix(p))
	}
	tau <- transitionPr(TAUP)

	arm <- integer(nr)
	statePath <- integer(nr)
	scaleFactor <- rep(0, nr)
	bv <- fv <- matrix(0.0, nr*S, 1L)
	initialProb <- rep(1/S, S)

	fitViterbi <- function(emit){
		emit <- as.matrix(as.numeric(emit))
		res <- .C("viterbi2", emit=emit, pi=initialProb, tau=tau,
			  arm=arm, S=S, nr=nr, statePath=statePath,
			  fv=fv, bv=bv, 1, 1, 1, 3L, scaleFactor=scaleFactor)[c("statePath", "fv", "bv", "scaleFactor")]
		statePath <- res[["statePath"]]
		sf <- res[["scaleFactor"]]
		fv <- res[["fv"]]
		bv <- res[["bv"]]
		res <- list(statePath=statePath, fv=fv, bv=bv, sf=sf)
	}

	## computeLLR
	CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")
	toGRanges <- function(statePath, j){
		id <- colnames(r)[j]
		if(is.null(id)) id <- paste("sample", length(j), sep="")
		rl <- Rle(statePath)
		starts <- position[start(rl)]
		ends <- position[end(rl)]
		states <- statePath[start(rl)]
		if(!is.null(id)){
			gr <- GRanges(CHR, IRanges(starts, ends), numberProbes=width(rl), state=states, sample=id)
		} else {
			gr <- GRanges(CHR, IRanges(starts, ends), numberProbes=width(rl), state=states)
		}
	}

	if(computeLLR){
		log.initial <- log(initialProb)
		tauC <- 1-tau
		lP.N2N <- log(1-(tauC*(S-1))) ##probability normal -> normal
		lP.N2A <- log(tauC) ##probability normal -> altered
		P.A2A <- sapply(1-(tauC*(1+(S-2))), function(x) max(x, 0.01))
		lP.A2A <- log(P.A2A) ## probability altered to same altered state
		lP.A2N <- lP.N2A ##probability altered -> normal
		lP.A2Astar <- lP.N2A ## probability altered -> different altered state
		## featureRanges
		fr <- GRanges(rep(CHR, length(position)), IRanges(position, width=1))
	}
	computeLogLikRatio <- function(gr, emit){
		gr2 <- gr
		alt.index <- which(state(gr) != normalIndex & numberProbes(gr) > 1)
		if(length(alt.index)==0) return(rep(0.0, length(gr)))
		gr <- gr[alt.index, ]
		log.emission <- log(emit)
		L <- length(gr)
		LLR <- rep(NA,  L)
		olaps <- findOverlaps(gr, fr)
		index <- subjectHits(olaps)
		indexList <- split(index, queryHits(olaps))
		starts <- sapply(indexList, min)
		ends <- sapply(indexList, max)
		statePath <- as.integer(state(gr))
		T <- length(statePath)
		rangeLogLik <- function(from, to, thisState){
			index <- seq(from, to)
			index2 <- index[-1]## t=2, ...., t*
			if(from == 1){
				if(to < T){
					logLik.vit <- log.initial[thisState] + log.emission[from, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) + lP.A2N[to+1] + log.emission[to+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) + lP.N2N[to+1] + log.emission[to+1, normalIndex]
				} else {
					logLik.vit <- log.initial[thisState] + log.emission[from, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) ##+ lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) ##+ lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
				}
			} else { ## first index > 1
				index2 <- index[-1]
				if(to < T){
					logLik.vit <- lP.N2A[from] +
						sum(lP.A2A[index2]) +
							lP.A2N[to+1] +
								log.emission[from, thisState] +
									sum(log.emission[index2, thisState]) +
										log.emission[to+1, normalIndex]
					logLik.null <-
						sum(lP.N2N[index]) +
							lP.N2N[to+1]  +
								sum(log.emission[index, normalIndex]) +
									log.emission[to+1, normalIndex]
				} else {
					logLik.vit <- lP.N2A[from] + log.emission[from, thisState] + sum(lP.A2A[index2] + log.emission[index2, thisState])
					logLik.null <- lP.N2N[from] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex])
				}
			}
			LLR <- logLik.vit - logLik.null
			return(LLR)
		}
		states <- state(gr)
		for(k in seq_len(L)) LLR[k] <- rangeLogLik(from=starts[k], to=ends[k], thisState=states[k])
		res <- rep(0.0, length(gr2))
		res[alt.index] <- LLR
		return(res)
	}

	list(initialCnParams=initialCnParams,
	     initialBafParams=initialBafParams,
	     updateBafEmission=updateBafEmission,
	     updateCnEmission=updateCnEmission,
	     updateBafParams=updateBafParams,
	     updateCnParams=updateCnParams,
	     trPr=transitionPr,
	     fitViterbi=fitViterbi,
	     toGRanges=toGRanges,
	     computeLogLikRatio=computeLogLikRatio)
}

getPrB <- function(pB=NULL){
	##n <- length(x)
	if(is.null(pB)){
		pB <- 0.5
	}
	pA <- 1-pB
	pAA <- pA^2
	pAB <- 2*pA*pB
	pBB <- 1-pAA-pAB
	pAAA <- pA^3
	pAAB <- 3*pA^2*pB
	pABB <- 3*pA*pB^2
	pBBB <- 1-pAAA-pAAB-pABB
	pAAAA <- pA^4
	pAAAB <- 4*pA^3*pB
	pAABB <- 6*pA^2*pB^2
	pABBB <- 4*pA*pB^3
	pBBBB <- pB^4
	p <- matrix(c(pB,
		      pA,
		      pAA,
		      pAB,
		      pBB,
		      pAAA,
		      pAAB,
		      pABB,
		      pBBB,
		      pAAAA,
		      pAAAB,
		      pAABB,
		      pABBB,
		      pBBBB),
		    byrow=FALSE, nrow=length(pB), ncol=14)
	colnames(p) <- c("B", "A", "AA", "AB", "BB",
			 "AAA", "AAB", "ABB", "BBB",
			 "AAAA", "AAAB", "AABB", "ABBB", "BBBB")
	return(p)
}

generatorFunG <- function(r, gt, is.snp, cnStates,
			  normalIndex, TAUP, limits, center,
			  position, is.log, computeLLR, chrom, verbose=FALSE){
	S <- length(cnStates)
	nc <- ncol(r)
	nr <- nrow(r)
	np.index <- which(!is.snp)
	CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")

	## params for copy number
	sds <- apply(r, 2, mad, na.rm=TRUE)
	if(center) r <- centerCopyNumber(r, is.snp) + cnStates[normalIndex]
	r <- thresholdCopyNumber(r, limits=limits)

	initialCnParams <- function(j){
		mus <- cnStates
		sigmas <- rep(sds[j], S)
		##p <- matrix(c(0.99,  0.01), S, 2, byrow=TRUE)
		p <- 0.01
		paramsCN <- list(mu=mus, sigmas=sigmas, p=p)
	}

	g2 <- g1 <- matrix(NA, nr, S)
	d1 <- matrix(NA, nr, S)
	emitr <- matrix(NA, nr, S)
	cn.dunif <- dunif(r, limits[1], limits[2])
	Gtotal <- matrix(NA, nr, 4)

	## might try dt here instead
	mydnorm <- function(x){
		function(mean, sd){
			mus <- matrix(mean, nr, S, byrow=TRUE)
			sds <- matrix(sd, nr, S, byrow=TRUE)
			dnorm(x, mus, sds)
		}
	}

	updateCnEmission <- function(cn, params, j){
		mus <- params[["mu"]]
		sds <- params[["sigmas"]]
		p <- params[["p"]]
		q <- 1-p
		normalDens <- mydnorm(cn)
		d <- normalDens(mus, sds)
		emitr <- q*d + p*cn.dunif[, j]
		emitr.na <- is.na(emitr)
		hasna <- any(emitr.na)
		if(hasna) emitr[emitr.na] <- 1
		return(emitr)
	}
	cn.sd.new <- rep(NA, S)
	updateCnParams <- function(cn, params, fv, bv, j) {
		mu <- params[["mu"]]
		sigma <- params[["sigmas"]]
		p <- params[["p"]]; q <- 1-p
		min.sd <- sds[j]/2

		h <- matrix(fv*bv, nr, S)
		m <- rowSums(h, na.rm=TRUE)
		h <- h/m

		d.unif <- cn.dunif[, j]
		d <- q*dnorm(cn, mean=matrix(mu, nr, S, byrow=TRUE), sd=matrix(sigma, nr, S, byrow=TRUE)) + p*d.unif
		d1 <- d/(d+d.unif)

		g1 <- h * d1
		g2 <- h * (1-d1)

		totalh <- apply(g1, 2, sum, na.rm=TRUE)
		mu.new <- colSums(g1*cn, na.rm=TRUE)/totalh
		mu.new[-4] <- constrainMu(mu.new[-4], is.log)
		mu.new[4] <- mu.new[normalIndex]
		mu.new <- makeNonDecreasing(mu.new)

		## For loop
		for(s in seq_len(S))  cn.sd.new[s] <- sqrt(sum(g1[, s]*(cn-mu.new[s])^2, na.rm=TRUE)/totalh[s])
		cn.sd.new <- pmax(cn.sd.new, min.sd)
		cn.sd.new <- constrainSd(cn.sd.new)
		params[["mu"]] <- mu.new
		params[["sigmas"]] <- cn.sd.new

		## update p as the proportion of values not fit by any of the states
		total.g2 <- apply(g2, 2, sum, na.rm=TRUE)
		denom.eq52 <- apply(g1 + g2, 2, sum, na.rm=TRUE)
		p <- min(total.g2/denom.eq52)  ## altered
		params[["p"]] <- p
		params
	}

	transitionPr <- function(TAUP, tauMAX=1-5e-8){
		p <- exp(-2*diff(position)/TAUP)
		minimum <- 1-1/((S-1)) + 0.01
		p <- pmax(p, minimum)
		p <- pmin(p, tauMAX)
		return(as.matrix(p))
	}
	tau <- transitionPr(TAUP)

	arm <- integer(nr)
	statePath <- integer(nr)
	scaleFactor <- rep(0, nr)
	bv <- fv <- matrix(0.0, nr*S, 1L) ## forward and backward variables
	initialProb <- rep(1/S, S)

	fitViterbi <- function(emit){
		emit <- as.matrix(as.numeric(emit))
		res <- .C("viterbi2", emit=emit, pi=initialProb, tau=tau,
			  arm=arm, S=S, nr=nr, statePath=statePath,
			  fv=fv, bv=bv, 1, 1, 1, 3L,
			  scaleFactor=scaleFactor)[c("statePath", "fv", "bv", "scaleFactor")]
		statePath <- res[["statePath"]]
		sf <- res[["scaleFactor"]]
		fv <- res[["fv"]]
		bv <- res[["bv"]]
		res <- list(statePath=statePath, fv=fv, bv=bv, sf=sf)
	}

	## computeLLR
	CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")
	toGRanges <- function(statePath, j){
		id <- colnames(r)[j]
		rl <- Rle(statePath)
		starts <- position[start(rl)]
		ends <- position[end(rl)]
		states <- statePath[start(rl)]
		GRanges(CHR, IRanges(starts, ends), numberProbes=width(rl), state=states, sample=id)
	}

	if(computeLLR){
		log.initial <- log(initialProb)
		tauC <- 1-tau
		lP.N2N <- log(1-(tauC*(S-1))) ##probability normal -> normal
		lP.N2A <- log(tauC) ##probability normal -> altered
		P.A2A <- sapply(1-(tauC*(1+(S-2))), function(x) max(x, 0.01))
		lP.A2A <- log(P.A2A) ## probability altered to same altered state
		lP.A2N <- lP.N2A ##probability altered -> normal
		lP.A2Astar <- lP.N2A ## probability altered -> different altered state
		## featureRanges
		fr <- GRanges(rep(CHR, length(position)), IRanges(position, width=1))
	}
	computeLogLikRatio <- function(gr, emit){
		gr2 <- gr
		alt.index <- which(state(gr) != normalIndex & numberProbes(gr) > 1)
		if(length(alt.index)==0) return(rep(0.0, length(gr)))
		gr <- gr[alt.index, ]
		log.emission <- log(emit)
		L <- length(gr)
		LLR <- rep(NA,  L)
		olaps <- findOverlaps(gr, fr)
		index <- subjectHits(olaps)
		indexList <- split(index, queryHits(olaps))
		starts <- sapply(indexList, min)
		ends <- sapply(indexList, max)
		statePath <- as.integer(state(gr))
		T <- length(statePath)
		rangeLogLik <- function(from, to, thisState){
			index <- seq(from, to)
			index2 <- index[-1]## t=2, ...., t*
			if(from == 1){
				if(to < T){
					logLik.vit <- log.initial[thisState] + log.emission[from, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) + lP.A2N[to+1] + log.emission[to+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) + lP.N2N[to+1] + log.emission[to+1, normalIndex]
				} else {
					logLik.vit <- log.initial[thisState] + log.emission[from, thisState]  + sum(lP.A2A[index2] + log.emission[index2, thisState]) ##+ lP.A2N[last.index+1] + log.emission[last.index+1, normalIndex]
					logLik.null <- log.initial[normalIndex] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex]) ##+ lP.N2N[last.index+1] + log.emission[last.index+1, normalIndex]
				}
			} else { ## first index > 1
				index2 <- index[-1]
				if(to < T){
					logLik.vit <- lP.N2A[from] +
						sum(lP.A2A[index2]) +
							lP.A2N[to+1] +
								log.emission[from, thisState] +
									sum(log.emission[index2, thisState]) +
										log.emission[to+1, normalIndex]
					logLik.null <-
						sum(lP.N2N[index]) +
							lP.N2N[to+1]  +
								sum(log.emission[index, normalIndex]) +
									log.emission[to+1, normalIndex]
				} else {
					logLik.vit <- lP.N2A[from] + log.emission[from, thisState] + sum(lP.A2A[index2] + log.emission[index2, thisState])
					logLik.null <- lP.N2N[from] + log.emission[from, normalIndex] + sum(lP.N2N[index2] + log.emission[index2, normalIndex])
				}
			}
			LLR <- logLik.vit - logLik.null
			return(LLR)
		}
		states <- state(gr)
		for(k in seq_len(L)) LLR[k] <- rangeLogLik(from=starts[k], to=ends[k], thisState=states[k])
		res <- rep(0.0, length(gr2))
		res[alt.index] <- LLR
		return(res)
	}
	list(initialCnParams=initialCnParams,
	     updateCnEmission=updateCnEmission,
	     updateCnParams=updateCnParams,
	     trPr=transitionPr,
	     fitViterbi=fitViterbi,
	     toGRanges=toGRanges,
	     computeLogLikRatio=computeLogLikRatio)
}

gcPercent <- function(feature.gr, width){
	##feature.gr <- feature.gr[chromosome(feature.gr) == chrom, ]
	bsgenomepkg <- paste("BSgenome.Hsapiens.UCSC.", metadata(feature.gr)$genome, sep="")
	library(bsgenomepkg, character.only=TRUE)
	Hsapiens <- get("Hsapiens")
	chrom <- chromosome(feature.gr)[1]
	chr <- Hsapiens[[chrom]]
	v <- suppressWarnings(Views(chr, start=start(feature.gr), end=end(feature.gr)))
	gcFunction <- function(x){
		alf <- alphabetFrequency(x, as.prob=TRUE)
		rowSums(alf[, c("G", "C"), drop=FALSE])
	}
	probeGC <- suppressWarnings(gcFunction(v))
	probeGC <- as.integer(probeGC*100)
	return(probeGC)
}

gcCorrect <- function(object, width){
	feature.gr <- makeFeatureGRanges(object)
	feature.gr <- sort(feature.gr)
	feature.gr <- resize(feature.gr, width=width, fix="center")
	metadata(feature.gr) <- list(genome=genomeBuild(object))

	arm <- getArm(object)
	arm <- factor(arm, levels=unique(arm))
	feature.grl <- split(feature.gr, arm)

	gcP <- lapply(feature.grl, gcPercent, width=width)
	gcP <- gcP[elementLengths(gcP) > 100]
	## not sure how many cuts we want to make...
	binGCcontent <- function(x) bins <- cut(x, breaks=seq(0,100, by=1), labels=FALSE, include.lowest=TRUE)
	##binGCcontent <- function(x) bins <- cut(x, breaks=quantile(x, probs=seq(0, 1, 0.1)), labels=FALSE, include.lowest=TRUE)
	binList <- lapply(gcP, binGCcontent)

	## Calculate mean lrr for each bin.
	## Subtract the mean lrr for the gc bins.  Add back the grand median so that the
	## center for the chromosome is unchanged
	avgRBin <- function(r, bin, grand.median) {
		bins <- as.integer(factor(bin))
		ravg <- sapply(split(r, bins), mean, na.rm=TRUE)
		## what about QN by middle bins..
		ravg <- ravg[bins]
		r.adj <- r - ravg + grand.median
		##boxplot(split(r/100, bins), pch=".")
		##boxplot(split(r.adj/100, bins), pch=".")
		return(r.adj)
	}
	for(j in seq_len(ncol(object))){
		rList <- split(lrr(object)[, j], arm)  ## can do this on the integer scale
		rList <- rList[elementLengths(rList) > 100]
		medianR <- sapply(rList, median, na.rm=TRUE)
		adjR <- mapply(avgRBin, rList, binList, medianR)
		rvec <- unlist(adjR)
		if(length(rvec) < nrow(object)){
			index <- which(arm %in% names(rList))
			lrr(object)[index, j] <- rvec
		} else lrr(object)[, j] <- rvec
	}
	return(object)
}


validChromosomeIndex <- function(object){
	index <- which(chromosome(object) <= 24 & !is.na(chromosome(object)) & !is.na(position(object)))
	if(length(index) < 1){
		stop("Chromosome must be 1-24 and physical position \n
                      can not be missing.  See chromosome() and position().\n
                      See integer2chromosome(1:24) for integer codes used by\n
                      VanillaICE.")
	}
	return(index)
}

makeMusBafNondecreasing <- function(mus){
	nms <- names(mus)
	index <- match(c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B"), names(mus))
	mus2 <- mus[index]
	mus2[2] <- if(mus2[2] > mus2[3]) mus2[3] else mus2[2]
	mus2[3] <- if(mus2[3] > mus2[4]) mus2[4] else mus2[3]
	mus2[5] <- if(mus2[5] < mus2[4]) mus2[4] else mus2[5]
	mus2[6] <- if(mus2[6] < mus2[5]) mus2[5] else mus2[6]
	mus2
}

##invalidCnConfidence <- function(x){
##	is.na(x) | x <= 0 | is.nan(x) | is.infinite(x)
##}






