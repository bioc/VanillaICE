generatorViterbiBeadStudioSet <- function(R, B, chr, center,
					  snp.index, anyNP,
					  cnStates){
	chr <- chr
	R <- R
	B <- B
	toMatrix <- function(j){
		r <- R[, j, drop=FALSE]/100
		b <- B[, j, drop=FALSE]/1000
		if(anyNP) b <- b[snp.index, , drop=FALSE]
		rownames(r) <- rownames(b) <- NULL
		if(center){
			if(all(as.integer(chr) <= 22)){
				meds <- colMedians(r, na.rm=TRUE)
				r2 <- sweep(r, 2, meds)
				r <- r2+cnStates[3]
			} else{
				autosome.index <- which(chr <= 22)
				rauto <- r[autosome.index, ]
				meds <- colMedians(rauto, na.rm=TRUE)
				rauto <- sweep(rauto, 2, meds)
				rauto <- rauto+cnStates[3]
				r[autosome.index, ] <- rauto
			}
		}
		close(R)
		close(B)
		list(r=r, b=b)
	}
	return(toMatrix)
}

generatorViterbiCnGt <- function(R, G, chr, center,
				 snp.index, anyNP,
				 cnStates){
	chr <- chr
	R <- R
	G <- G
	toMatrix <- function(j){
		r <- R[, j, drop=FALSE]/100
		g <- G[, j, drop=FALSE]
		if(anyNP) g <- g[snp.index, , drop=FALSE]
		rownames(r) <- rownames(g) <- NULL
		if(center){
			if(all(as.integer(chr) <= 22)){
				meds <- colMedians(r, na.rm=TRUE)
				r2 <- sweep(r, 2, meds)
				r <- r2+cnStates[3]
			} else{
				autosome.index <- which(chr <= 22)
				rauto <- r[autosome.index, ]
				meds <- colMedians(rauto, na.rm=TRUE)
				rauto <- sweep(rauto, 2, meds)
				rauto <- rauto+cnStates[3]
				r[autosome.index, ] <- rauto
			}
		}
		close(R)
		close(G)
		list(r=r, g=g)
	}
	return(toMatrix)
}

generatorViterbiSnpSet <- function(G, chr, isff){
	chr <- chr
	G <- G
	toMatrix <- function(j){
		g <- G[, j, drop=FALSE]
		if(isff) close(G)
		return(g)
	}
	return(toMatrix)
}

generatorViterbiSnpSetIce <- function(G, GP, chr, isff){
	chr <- chr
	G <- G
	GP <- GP
	toMatrix <- function(j){
		g <- G[, j, drop=FALSE]
		gp <- GP[, j, drop=FALSE]
		gp <- i2p(gp)
		if(isff) {
			close(G)
			close(GP)
		}
		return(list(g=g, gt.conf=gp))
	}
	return(toMatrix)
}

generatorGRanges <- function(chrom, position, build, ids, TAUP, tauMAX){
	S <- 6
	CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")
	chrarm <- oligoClasses:::.getArm(as.integer(chrom), position, build)
	chrarm <- factor(chrarm, unique(chrarm))
	sl <- getSequenceLengths(build)
	sl <- sl[unique(CHR)]

	toGRanges <- function(statePath, id){##, j){
		sp <- split(statePath, chrarm)
		rlist <- lapply(sp, Rle)
		pos <- split(position, chrarm)
		x <- x2 <- x1 <- rl <- states <- NULL
		starts <- foreach(rl=rlist, x=pos) %do% x[start(rl)]
		ends <- foreach(rl=rlist, x=pos) %do% x[end(rl)]
		statelist <- foreach(states=sp, rl=rlist) %do% states[start(rl)]
		chrlist <- split(CHR, chrarm)
		chrl <- foreach(chrom=chrlist, rl=rlist) %do% chrom[start(rl)]
		gr <- foreach(states=statelist,
			      x1=starts,
			      x2=ends,
			      rl=rlist,
			      chrom=chrl) %do% {
				      GRanges(chrom, IRanges(x1, x2),
					      numberProbes=width(rl),
					      sample=id,
					      state=states,
					      seqlengths=sl)
			      }
		gr <- unlist(GRangesList(gr))
	}
	if(missing(tauMAX)) tauMAX <- 1-5e-8
	d <- diff(position)
	p <- exp(-2*d/TAUP)
	minimum <- 1-1/((S-1)) + 0.01
	p <- pmax(p, minimum)
	p <- pmin(p, tauMAX)
	##
	## 1-(1-tau)*(S-1)*c1, c1=1.  can't be less than 0.8
	##
	if(any(d < 0)) p[d < 0] <- 0.8 ## different chromosomes
	tau <- p; rm(p)

	initialProb <- rep(1/S, S)
	##log.initial <- log(initialProb)
	tauC <- 1-tau ## probability that state t+1 != state t
	lP.N2N <- log(tau)
	lP.N2A <- log(tauC) ##probability normal -> altered
	lP.A2A <- lP.N2N
	lP.A2N <- lP.N2A ##probability altered -> normal
	lP.A2Astar <- lP.N2A ## probability altered -> different altered state
	if(length(CHR)==1){
		fr <- GRanges(rep(CHR, length(position)), IRanges(position, width=1))
	} else fr <- GRanges(CHR, IRanges(position, width=1))
	transitionProbs <- list(lP.N2N=lP.N2N,
				lP.N2A=lP.N2A,
				lP.A2A=lP.A2A,
				lP.A2N=lP.A2N,
				lP.A2Astar=lP.A2Astar,
				fr=fr)
	funs <- list(toGRanges=toGRanges,
				 tau=tau,
				 transitionProbs=transitionProbs)
	return(funs)
}

generatorViterbi <- function(rlist, blist, chr, center, snp.index, anyNP){
	chr <- chr
	rlist <- rlist
	blist <- blist
	toMatrix <- function(j){
		rl <- lapply(rlist, function(x, j) x[, j, drop=FALSE], j=j)
		bl <- lapply(blist, function(x, j) x[, j, drop=FALSE], j=j)
		r <- do.call("rbind", rl)/100
		b <- do.call("rbind", bl)/1000
		if(anyNP) b <- b[snp.index, , drop=FALSE]
		rownames(r) <- rownames(b) <- NULL
		if(center){
			if(all(as.integer(chr) <= 22)){
				meds <- colMedians(r, na.rm=TRUE)
				r <- sweep(r, 2, meds)
			} else{
				autosome.index <- which(chr <= 22)
				rauto <- r[autosome.index, , drop=FALSE]
				meds <- colMedians(rauto, na.rm=TRUE)
				r[autosome.index, ] <- sweep(rauto, 2, meds)
			}
		}
		sapply(rlist, close)
		sapply(blist, close)
		list(r=r, b=b)
	}
	return(toMatrix)
}


generatorFun <- function(r, b, gt, snp.index, cnStates,
			 normalIndex, tau, limits, center,
			 prOutlierBAF, p.hom, is.log,
			 computeLLR, verbose=FALSE,
			 transitionProbs){
	S <- length(cnStates)
	nc <- ncol(r)
	## remove nonpolymorphic markers
	## b <- b[snp.index, , drop=FALSE]
	nr <- nrow(r)
	nb <- nrow(b)
	##np.index <- which(!is.snp)
	names(prOutlierBAF)[1] <- "prOutlier"
	##CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")
	##
	## params for copy number
	sds <- rep(NA, ncol(r))
	for(j in seq_along(sds)){
		i <- r[, j] != 0
		sds[j] <- sd(r[i,j], na.rm=TRUE)
	}
	##sds <- apply(r, 2, mad, na.rm=TRUE)
	##if(center) r <- centerCopyNumber(r, is.snp) + cnStates[normalIndex]
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
	bindex <- b[, 1] > 0 & b[, 1] < 1
	initialBafParams <- function(){
		musBAF <- c(0, 0.1, 1/3, 0.5, 2/3, 0.9, 1)
		## we have already removed the nonpolymorphic markers
		index <- b[, 1] < 0.25
		sdA <- sd(b[index & bindex, 1], na.rm=TRUE)
		sdsBAF <- c(sdA, rep(sdA*1.5, 5), sdA)
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
		bindex <- bf > 0 & bf < 1
		weightedMoments <- function(w, homozygousA=FALSE, homozygousB=FALSE) {
			if(sum(w > 0.5, na.rm=TRUE) >= 10){
				##w[w < 0.01 | is.na(w)] <- 0
				w[is.na(w)] <- 0
				w <- w[bindex]
				bf <- bf[bindex]
				sum.w <- sum(w)
				if(homozygousA) mu <- 0
				if(homozygousB) mu <- 1
				if(!homozygousA & !homozygousB)
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
			moments <- matrix(NA, length(nms), 2)
			for(i in seq_len(length(nms))){
				homozygousA <- nms[i] == "A"
				homozygousB <- nms[i] == "B"
				moments[i, ] <- weightedMoments(gammas[[i]],
								homozygousA=homozygousA,
								homozygousB=homozygousB)
			}
##			if(nms==c("A", "B")){
##				moments <- matrix(NA, 2, 2)
##				moments[1, ] <- weightedMoments(gammas, homozygousA=TRUE)
##				moments[2, ] <- weightedMoments(gammas, homozygousB=TRUE)
##			} else {
##				moments <- t(sapply(gammas, weightedMoments))
##			}
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
		rmax <- tryCatch(rowMax(Gtotal), error=function(e) NULL)
		if(!is.null(rmax)){
			isout <- rmax < 0.9
			p <- mean(diff(isout) != 0)
		} else p <- pOut
		params[["mus"]] <- mus.new
		params[["sigmas"]] <- sigmas.new
		pOut[["prOutlier"]] <- min(p, 0.001)
		params[["prOutlier"]] <- pOut
		return(params)
	}
	cn.sd.new <- rep(NA, S)
	updateCnParams <- function(cn, params, h, j) {
		lindex <- cn != 0
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
		if(diploid > totalh[2]){
			mu.new[2] <- min(mu.new[3]-one.sd, mu.new[2])
		} else{
			mu.new[3] <- max(mu.new[2]+one.sd, mu.new[3])
		}
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
		if(any(is.na(cn.sd.new))) {
			naindex <- which(is.na(cn.sd.new))
			cn.sd.new[naindex] <- params$sigmas[naindex]
		}
		cn.sd.new <- pmax(cn.sd.new, min.sd)
		##if(any(!estimable))
		##	cn.sd.new[!estimable] <- params$sigmas[!estimable]
		params[["mu"]] <- mu.new
		params[["sigmas"]] <- cn.sd.new
		##total.g2 <- apply(g2, 2, sum, na.rm=TRUE)
		##denom.eq52 <- apply(g1 + g2, 2, sum, na.rm=TRUE)
		p <- mean(rowSums(g1, na.rm=TRUE) < 0.1)
		##p <- min(total.g2/denom.eq52)  ## altered states will have higher values of total.g2/denom.eq52.
		params[["p"]] <- min(p, 0.1)
		return(params)
	}
	##tau <- transitionPr(TAUP)
	arm <- integer(nr)
	statePath <- integer(nr)
	scaleFactor <- rep(0, nr)
	bv <- fv <- matrix(0.0, nr*S, 1L)
	initialProb <- rep(1/S, S)
	fitViterbi <- function(emit){
		emit <- as.matrix(as.numeric(emit))
		res <- .C("viterbi2",
			  emit=emit,
			  pi=initialProb,
			  tau=tau,
			  arm=arm,
			  S=S,
			  nr=nr,
			  statePath=statePath,
			  fv=fv,
			  bv=bv,
			  1,
			  1,
			  1,
			  3L,
			  scaleFactor=scaleFactor)[c("statePath", "fv", "bv", "scaleFactor")]
		statePath <- res[["statePath"]]
		sf <- res[["scaleFactor"]]
		fv <- res[["fv"]]
		bv <- res[["bv"]]
		res <- list(statePath=statePath, fv=fv, bv=bv, sf=sf)
	}

	## this should not be within the for loop.
	## computeLLR
##	CHR <- paste("chr", oligoClasses::integer2chromosome(chrom), sep="")
##	chrarm <- oligoClasses:::.getArm(as.integer(chrom), position, build)
##	chrarm <- factor(chrarm, unique(chrarm))
##	toGRanges <- function(statePath, j){
##		id <- colnames(r)[j]
##		if(is.null(id)) id <- paste("sample", length(j), sep="")
##		sp <- split(statePath, chrarm)
##		rlist <- lapply(sp, Rle)
##		pos <- split(position, chrarm)
##		##rl <- Rle(statePath)
##		starts <- foreach(rl=rlist, x=pos) %do% x[start(rl)]
##		ends <- foreach(rl=rlist, x=pos) %do% x[end(rl)]
##		##starts <- position[start(rl)]
##		##ends <- position[end(rl)]
##		statelist <- foreach(states=sp, rl=rlist) %do% states[start(rl)]
##		chrlist <- split(CHR, chrarm)
##		chrl <- foreach(chrom=chrlist, rl=rlist) %do% chrom[start(rl)]
##		gr <- foreach(states=statelist,
##			      x1=starts,
##			      x2=ends,
##			      rl=rlist,
##			      chrom=chrl) %do% {
##				      GRanges(chrom, IRanges(x1, x2),
##					      numberProbes=width(rl), sample=id,
##					      state=states)
##			      }
##		gr <- unlist(GRangesList(gr))
##		##states <- statePath[start(rl)]
####		if(!is.null(id)){
####			gr <- GRanges(CHR, IRanges(starts, ends), numberProbes=width(rl), state=states, sample=id)
####		} else {
####			gr <- GRanges(CHR, IRanges(starts, ends), numberProbes=width(rl), state=states)
####		}
##	}
	log.initial <- log(rep(1/S, S))

	computeLogLikRatio <- function(gr, emit){
		lP.N2N <- transitionProbs$lP.N2N
		lP.N2A <- transitionProbs$lP.N2A
		lP.A2A <- transitionProbs$lP.A2A
		lP.A2N <- transitionProbs$lP.A2N
		lP.A2Astar <- transitionProbs$lP.A2Astar
		fr <- transitionProbs$fr
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
			result <- logLik.vit - logLik.null
			return(result)
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
	     ##trPr=transitionPr,
	     fitViterbi=fitViterbi,
	     ##toGRanges=toGRanges,
	     computeLogLikRatio=computeLogLikRatio)
}

generatorFunG2 <- function(r, cnStates,
			   normalIndex, tau, limits,
			   p.hom, is.log,
			   computeLLR, verbose=FALSE,
			   transitionProbs){
	## why are we passing g?
	S <- length(cnStates)
	nc <- ncol(r)
	nr <- nrow(r)

	## params for copy number
	sds <- rep(NA, ncol(r))
	for(j in seq_along(sds)){
		i <- r[, j] != 0
		sds[j] <- sd(r[i,j], na.rm=TRUE)
	}

	initialCnParams <- function(j){
		mus <- cnStates
		sigmas <- rep(sds[j], S)
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
	#
	arm <- integer(nr)
	statePath <- integer(nr)
	scaleFactor <- rep(0, nr)
	bv <- fv <- matrix(0.0, nr*S, 1L)
	initialProb <- rep(1/S, S)

	fitViterbi <- function(emit){
		emit <- as.matrix(as.numeric(emit))
		res <- .C("viterbi2", emit=emit, pi=initialProb,
			  tau=tau,
			  arm=arm, S=S, nr=nr, statePath=statePath,
			  fv=fv, bv=bv, 1, 1, 1, 3L,
			  scaleFactor=scaleFactor)[c("statePath", "fv", "bv", "scaleFactor")]
		statePath <- res[["statePath"]]
		sf <- res[["scaleFactor"]]
		fv <- res[["fv"]]
		bv <- res[["bv"]]
		res <- list(statePath=statePath, fv=fv, bv=bv, sf=sf)
	}
	log.initial <- log(rep(1/S, S))
	computeLogLikRatio <- function(gr, emit){
		lP.N2N <- transitionProbs$lP.N2N
		lP.N2A <- transitionProbs$lP.N2A
		lP.A2A <- transitionProbs$lP.A2A
		lP.A2N <- transitionProbs$lP.A2N
		lP.A2Astar <- transitionProbs$lP.A2Astar
		fr <- transitionProbs$fr
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
	     fitViterbi=fitViterbi,
	     computeLogLikRatio=computeLogLikRatio)
}

generatorFunSnpSet <- function(g, normalIndex, tau,
			       computeLLR, verbose=FALSE,
			       transitionProbs,
			       S){
	## why are we passing g?
	nc <- ncol(g)
	nr <- nrow(g)

	#
	arm <- integer(nr)
	statePath <- integer(nr)
	scaleFactor <- rep(0, nr)
	bv <- fv <- matrix(0.0, nr*S, 1L)
	initialProb <- rep(1/S, S)

	fitViterbi <- function(emit){
		emit <- as.matrix(as.numeric(emit))
		res <- .C("viterbi2", emit=emit, pi=initialProb,
			  tau=tau,
			  arm=arm, S=S, nr=nr, statePath=statePath,
			  fv=fv, bv=bv, 1, 1, 1, 3L,
			  scaleFactor=scaleFactor)[c("statePath", "fv", "bv", "scaleFactor")]
		statePath <- res[["statePath"]]
		sf <- res[["scaleFactor"]]
		fv <- res[["fv"]]
		bv <- res[["bv"]]
		res <- list(statePath=statePath, fv=fv, bv=bv, sf=sf)
	}
	log.initial <- log(rep(1/S, S))
	computeLogLikRatio <- function(gr, emit){
		lP.N2N <- transitionProbs$lP.N2N
		lP.N2A <- transitionProbs$lP.N2A
		lP.A2A <- transitionProbs$lP.A2A
		lP.A2N <- transitionProbs$lP.A2N
		lP.A2Astar <- transitionProbs$lP.A2Astar
		fr <- transitionProbs$fr
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
	list(fitViterbi=fitViterbi,
	     computeLogLikRatio=computeLogLikRatio)
}
