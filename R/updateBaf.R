

updateBaf <- function(x,
		      cnStates,
		      mus,
		      sigmas,
		      normalIndex,
		      rohIndex=normalIndex+1,
		      nUpdates=3,
		      fv,
		      bv,
		      prOutlier,
		      prOutlierMax,
		      prB){ ## probability outlier
	if(!missing(prB)){
		prB <- prB/100
		prB[is.na(prB)] <- 0.5
	} else prB <- rep(0.5, length(x))
	if(missing(mus)){
		mus <- c(0, 1/4, 1/3, 1/2, 2/3, 3/4, 1)
	}
	nms <- c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B")
	names(mus) <- nms
	if(missing(sigmas)){
		sigmas <- c(0.02, 0.05, 0.05, 0.05, 0.05 , 0.05, 0.02)
	}
	names(sigmas) <- names(mus)
	if(nUpdates==0) return(mus)
	s <- sigmas
	S <- length(cnStates)
	e <- cbind(1-prOutlier, prOutlier)
	g2 <- g1 <- matrix(NA, length(x), S)
	missing.fv <- missing(fv)
	if(!missing.fv){
		nr <- if(is(x, "matrix")) nrow(x) else length(x)
		fv <- matrix(fv, nr, S)
		bv <- matrix(bv, nr, S)
		h <- fv*bv
		m <- matrix(rowSums(h, na.rm=TRUE), length(x), S, byrow=FALSE)
		h <- h/m
		##rm(m, fv, bv);
	}
	d1 <- matrix(NA, length(x), S)
	na.index <- which(is.na(x))
	anyNA <- length(na.index) > 0
	pB <- prB
	pA <- 1-prB
	pAA <- pA^2
	pAB <- 2*pA*pB
	pBB <- 1-pAA-pAB
	pAAA <- pA^3
	pAAB <- 3*pA^2*pB
	pABB <- 3*pA*pB^2
	pBBB <- pA^3-pAAA-pAAB-pABB
	pAAAA <- pA^4
	pAAAB <- 4*pA^3*pB
	pAABB <- 6*pA^2*pB^2
	pABBB <- 4*pA*pB^3
	pBBBB <- pB^4
	for(iter in seq_len(nUpdates)){
		if(iter > nUpdates) break()
		d.unif <- dunif(x, 0,1)
		d0 <- d.unif
		d <- matrix(NA, length(x), length(mus), dimnames=list(NULL, names(mus)))
		d[, "A"] <- tnorm(x, mus["A"], sigmas["A"])
		d[, "B"] <- tnorm(x, mus["B"], sigmas["B"])
		d[, "AB"] <- tnorm(x, mus["AB"], sigmas["AB"])
		d[, "AAB"] <- tnorm(x, mus["AAB"], sigmas["AAB"])
		d[, "ABB"] <- tnorm(x, mus["ABB"], sigmas["ABB"])
		d[, "AAAB"] <- tnorm(x, mus["AAAB"], sigmas["AAAB"])
		d[, "ABBB"] <- tnorm(x, mus["ABBB"], sigmas["ABBB"])
		##emission probabilities
		p1 <- e[2,1]*pA
		p2 <- e[2,1]*pB
		p3 <- e[2,2]
		mA <- p1*d[, "A"]/(p1*d[, "A"] + p2*d[, "B"] + p3)
		mB <- p2*d[, "B"]/(p1*d[, "A"] + p2*d[, "B"] + p3)
		m1out <- 1-mA-mB
		p1 <- e[3,1]*pAA
		p2 <- e[3,1]*pAB
		p3 <- e[3,1]*pBB
		p4 <- e[3,2]
		mAA <- p1*d[, "A"]/(p1*d[, "A"] + p2*d[,"AB"]+p3*d[, "B"] +p4)
		mAB <- p2*d[, "AB"]/(p1*d[, "A"] + p2*d[,"AB"]+p3*d[, "B"] +p4)
		mBB <- p3*d[, "B"]/(p1*d[, "A"] + p2*d[,"AB"]+p3*d[, "B"] +p4)
		m2out <- 1-mAA-mAB-mBB
		p1 <- e[5,1]*pAAA
		p2 <- e[5,1]*pAAB
		p3 <- e[5,1]*pABB
		p4 <- e[5,1]*pBBB
		p5 <- e[5,2]
		denom <- (p1*d[, "A"] + p2*d[, "AAB"] + p3*d[, "ABB"] + p4*d[, "B"] + p5)
		mAAA <- p1*d[, "A"]/denom
		mAAB <- p2*d[, "AAB"]/denom
		mABB <- p3*d[, "ABB"]/denom
		mBBB <- p4*d[, "B"]/denom
		m3out <- 1-mAAA-mAAB-mABB-mBBB
		p1 <- e[6,1]*pAAAA
		p2 <- e[6,1]*pAAAB
		p3 <- e[6,1]*pAABB
		p4 <- e[6,1]*pABBB
		p5 <- e[6,1]*pBBBB
		p6 <- e[6,2]
		denom <- (p1*d[, "A"] + p2*d[, "AAAB"] + p3*d[, "AB"] + p4*d[, "ABBB"] + p5*d[, "B"] + p6)
		mAAAA <- p1*d[,"A"]/denom
		mAAAB <- p2*d[,"AAAB"]/denom
		mAABB <- p3*d[,"AB"]/denom
		mABBB <- p4*d[,"ABBB"]/denom
		mBBBB <- p5*d[,"B"]/denom
		m4out <- 1-mAAAA-mAAAB-mAABB-mABBB-mBBBB
		##.index=subjectHits(findOverlaps(ir[1,], featureData(trioSet)))
		##h[.index, ]
		g0 <- h[, 1]*d0
		g1A <- h[, 2]*mA
		g1B <- h[, 2]*mB
		g1out <- h[, 2]*m1out
		g2AA <- h[, 3]*mAA
		g2AB <- h[, 3]*mAB
		g2BB <- h[, 3]*mBB
		g2out <- h[,3]*m2out
		g2rohAA <- h[, 4]*mA
		g2rohBB <- h[, 4]*mB
		g3AAA <- h[,5]*mAAA
		g3AAB <- pmax(h[,5]*mAAB,1e-9)
		g3ABB <- pmax(h[,5]*mABB,1e-9)
		g3BBB <- h[,5]*mBBB
		g3out <- h[,5]*m3out
		g4AAAA <- h[,6]*mAAAA
		g4AAAB <- pmax(h[,6]*mAABB, 1e-9)
		g4AABB <- pmax(h[,6]*mAABB, 1e-9)
		g4ABBB <- pmax(h[,6]*mABBB, 1e-9)
		g4BBBB <- h[,6]*mBBBB
		g4out <- h[, 6]*m4out
		sigma.new <- mu.new <- rep(NA, length(mus))
		names(sigma.new) <- names(mu.new) <- names(mus)
		mu.new["A"] <- sum(g1A*x, na.rm=TRUE)/sum(g1A, na.rm=TRUE)
		mu.new["B"] <- sum(g1B*x, na.rm=TRUE)/sum(g1B, na.rm=TRUE)
		mu.new["AB"] <- sum(g2AB*x,na.rm=TRUE)/sum(g2AB, na.rm=TRUE)
		mu.new["AAB"] <- sum(g3AAB*x,na.rm=TRUE)/sum(g3AAB, na.rm=TRUE)
		mu.new["ABB"] <- sum(g3ABB*x, na.rm=TRUE)/sum(g3ABB, na.rm=TRUE)
		mu.new["AAAB"] <- sum(g4AAAB*x,na.rm=TRUE)/sum(g4AAAB, na.rm=TRUE)
		mu.new["ABBB"] <- sum(g4ABBB*x,na.rm=TRUE)/sum(g4ABBB, na.rm=TRUE)
		sigma.new["A"] <- sqrt(sum(g1A*(x-mu.new["A"])^2, na.rm=TRUE)/sum(g1A, na.rm=TRUE))
		sigma.new["B"] <- sqrt(sum(g1B*(x-mu.new["B"])^2, na.rm=TRUE)/sum(g1B, na.rm=TRUE))
		sigma.new["AB"] <- sqrt(sum(g2AB*(x-mu.new["AB"])^2, na.rm=TRUE)/sum(g2AB, na.rm=TRUE))
		sigma.new["AAB"] <- sqrt(sum(g3AAB*(x-mu.new["AAB"])^2, na.rm=TRUE)/sum(g3AAB, na.rm=TRUE))
		sigma.new["ABB"] <- sqrt(sum(g3ABB*(x-mu.new["ABB"])^2, na.rm=TRUE)/sum(g3ABB, na.rm=TRUE))
		sigma.new["AAAB"] <- sqrt(sum(g4AAAB*(x-mu.new["AAAB"])^2, na.rm=TRUE)/sum(g4AAAB, na.rm=TRUE))
		sigma.new["ABBB"] <- sqrt(sum(g4ABBB*(x-mu.new["ABBB"])^2, na.rm=TRUE)/sum(g4ABBB, na.rm=TRUE))
		mus <- mu.new
		mus <- constrainMuBaf(mus)
		mus <- makeMusBafNondecreasing(mus)
		sigmas <- sigma.new
		## ensures that .05, 0.95 are at least the .975 or 0.025 quantile
		sigmas <- pmax(sigmas, 0.025)
		p1bar <- sum(g1out,na.rm=TRUE)/sum(g1A+g1B+g1out, na.rm=TRUE)
		e[4, 2] <- e[2, 2] <- p1bar
		p2bar <- sum(g2out,na.rm=TRUE)/sum(g2AA+g2AB+g2BB+g2out, na.rm=TRUE)
		e[3, 2] <- p2bar
		p3bar <- sum(g3out,na.rm=TRUE)/sum(g3AAA+g3AAB+g3ABB+g3BBB+g3out, na.rm=TRUE)
		e[5, 2] <- p3bar
		p4bar <- sum(g4out,na.rm=TRUE)/sum(g4AAAA+g4AAAB+g4AABB+g4ABBB+g4BBBB+g4out, na.rm=TRUE)
		e[6, 2] <- p4bar
		pr.min <- max(e[, 2], 1e-5)
		e[, 2] <- pmax(e[, 2], pr.min)
		e[, 2] <- pmin(e[, 2], prOutlierMax)
		e[, 1] <- 1-e[, 2]
	}
	names(mus) <- names(sigmas) <- nms
	return(list(mu=mus, sigma=sigmas, p=e[, 2]))
}

constrainMuBaf <- function(mus){
	mus["A"] <- ifelse(mus["A"] > 0.005, 0.005, mus["A"])
	mus["B"] <- ifelse(mus["B"] < 0.995, 0.995, mus["B"])
	mus["AAAB"] <- ifelse(mus["AAAB"] > 0.35, 0.35, mus["AAAB"])
	mus["AAAB"] <- ifelse(mus["AAAB"] < 0.1, 0.1, mus["AAAB"])
	mus["AAB"] <- ifelse(mus["AAB"] < 0.2, 0.2, mus["AAB"])
	mus["AAB"] <- ifelse(mus["AAB"] > 0.45, 0.45, mus["AAB"])
	mus["AB"] <- ifelse(mus["AB"] < 0.45, 0.44, mus["AB"])
	mus["AB"] <- ifelse(mus["AB"] > 0.55, 0.55, mus["AB"])
	mus["ABBB"] <- ifelse(mus["ABBB"] > 0.9, 0.9, mus["ABBB"])
	mus["ABBB"] <- ifelse(mus["ABBB"] < 0.65, 0.65, mus["ABBB"])
	mus["ABB"] <- ifelse(mus["ABB"] > 0.8, 0.8, mus["ABB"])
	mus["ABB"] <- ifelse(mus["ABB"] < 0.55, 0.55, mus["ABB"])
	## symmetry around AB
	delta1 <- mus["AB"] - mus["AAB"]
	delta2 <- mus["ABB"] - mus["AB"]
	avg <- (delta1+delta2)/2
	mus["AAB"] <- mus["AB"]-avg
	mus["ABB"] <- mus["AB"]+avg
	delta1 <- mus["AB"] - mus["AAAB"]
	delta2 <- mus["ABBB"] - mus["AB"]
	avg <- (delta1+delta2)/2
	mus["AAAB"] <- mus["AB"] - avg
	mus["ABBB"] <- mus["AB"] + avg
	mus
}

makeMusBafNondecreasing <- function(mus){
	nms <- names(mus)
	index <- match(c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B"), names(mus))
	mus2 <- mus[index]
	mus2[2] <- if(mus2[2] > mus2[3]) mus2[3] else mus2[2]
	mus2[3] <- if(mus2[3] > mus2[4]) mus2[4] else mus2[3]
	mus2[5] <- if(mus2[5] < mus2[4]) mus2[4] else mus2[5]
	mus2[6] <- if(mus2[6] < mus2[5]) mus2[5] else mus2[6]
	##mus2 <- makeNonDecreasing(mus[index])
	mus <- mus2[match(nms, names(mus2))]
	mus
}
isNonDecreasing <- function(mus){
	nms <- names(mus)
	index <- match(c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B"), names(mus))
	all(diff(mus[index]) > 0)
}

