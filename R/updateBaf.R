## creator function



##updateBaf <- function(x,
##		      cnStates,
##		      mus,
##		      sigmas,
##		      normalIndex,
##		      rohIndex=normalIndex+1,
##		      nUpdates=3,
##		      fv,
##		      bv,
##		      prOutlierBAF,
##		      prB,
##		      constraint,
##		      allele.prob){
##	constrainMuBaf <- constraint[["mus.baf"]]
##	constrainSdBaf <- constraint[["sds.baf"]]
##	if(missing(mus)){
##		mus <- c(0, 1/4, 1/3, 1/2, 2/3, 3/4, 1)
##	}
##	nms <- names(mus) <- c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B")
##	if(missing(sigmas)){
##		sigmas <- c(0.02, 0.05, 0.05, 0.05, 0.05 , 0.05, 0.02)
##	}
##	names(sigmas) <- names(mus)
##	if(nUpdates==0) return(mus)
##	s <- sigmas
##	S <- length(cnStates)
##	## outlier probabilities
##	e <- cbind(1-prOutlierBAF[["initial"]], prOutlierBAF[["initial"]])
##	## ensure initial values satisfy constraints
##	pOutlierMax <- prOutlierBAF[["max"]]
##	pOutlierMaxROH <- prOutlierBAF[["maxROH"]]
##	e[, 2] <- pmin(e[, 2], pOutlierMax)
##	e[2,2] <- min(e[2,2], pOutlierMaxROH)
##	g2 <- g1 <- matrix(NA, length(x), S)
##	missing.fv <- missing(fv)
##	if(!missing.fv){
##		nr <- if(is(x, "matrix")) nrow(x) else length(x)
##		h <- matrix(fv*bv, nr, S)
##		m <- rowSums(h, na.rm=TRUE)
##		##m <- matrix(rowSums(h, na.rm=TRUE), length(x), S, byrow=FALSE)
##		h <- h/m
##		##rm(m, fv, bv);
##	}
##	d1 <- matrix(NA, length(x), S)
##	na.index <- which(is.na(x))
##	anyNA <- length(na.index) > 0
##	ishet <- which(x > 0.2 & x < 0.8)
##	emp.sd <- sd(x[ishet], na.rm=TRUE)
##	d <- matrix(NA, length(x), length(mus), dimnames=list(NULL, names(mus)))
##	outBounds.index <- which(x < 0 | x > 1)
##	is.out <- length(outBounds.index) > 0
##	trNormal <- trnorm(x)
##	for(iter in seq_len(nUpdates)){
##		d.unif <- dunif(x, 0,1)
##		d0 <- d.unif
##		d[, "A"] <- trNormal(mus["A"], sigmas["A"])
##		d[, "B"] <- trNormal(mus["B"], sigmas["B"])
##		d[, "AB"] <- trNormal(mus["AB"], sigmas["AB"])
##		d[, "AAB"] <- trNormal(mus["AAB"], sigmas["AAB"])
##		d[, "ABB"] <- trNormal(mus["ABB"], sigmas["ABB"])
##		d[, "AAAB"] <- trNormal(mus["AAAB"], sigmas["AAAB"])
##		d[, "ABBB"] <- trNormal(mus["ABBB"], sigmas["ABBB"])
##		if(is.out) d[outBounds.index, ] <- 0
##		##emission probabilities
##		p1 <- e[2,1]*allele.prob[, "A"]
##		p2 <- e[2,1]*allele.prob[, "B"]
##		p3 <- e[2,2]
##		denom <- (p1*d[, "A"] + p2*d[, "B"] + p3)
##		mA <- p1*d[, "A"]/denom
##		mB <- p2*d[, "B"]/denom
##		m1out <- 1-mA-mB
##		p1 <- e[3,1]*allele.prob[, "AA"]
##		p2 <- e[3,1]*allele.prob[, "AB"]
##		p3 <- e[3,1]*allele.prob[, "BB"]
##		p4 <- e[3,2]
##		denom <- (p1*d[, "A"] + p2*d[,"AB"]+p3*d[, "B"] +p4)
##		mAA <- p1*d[, "A"]/denom
##		mAB <- p2*d[, "AB"]/denom
##		mBB <- p3*d[, "B"]/denom
##		m2out <- 1-mAA-mAB-mBB
##		p1 <- e[5,1]*allele.prob[, "AAA"]
##		p2 <- e[5,1]*allele.prob[, "AAB"]
##		p3 <- e[5,1]*allele.prob[, "ABB"]
##		p4 <- e[5,1]*allele.prob[, "BBB"]
##		p5 <- e[5,2]
##		denom <- (p1*d[, "A"] + p2*d[, "AAB"] + p3*d[, "ABB"] + p4*d[, "B"] + p5)
##		mAAA <- p1*d[, "A"]/denom
##		mAAB <- p2*d[, "AAB"]/denom
##		mABB <- p3*d[, "ABB"]/denom
##		mBBB <- p4*d[, "B"]/denom
##		m3out <- 1-mAAA-mAAB-mABB-mBBB
##		p1 <- e[6,1]*allele.prob[, "AAAA"]
##		p2 <- e[6,1]*allele.prob[, "AAAB"]
##		p3 <- e[6,1]*allele.prob[, "AABB"]
##		p4 <- e[6,1]*allele.prob[, "ABBB"]
##		p5 <- e[6,1]*allele.prob[, "BBBB"]
##		p6 <- e[6,2]
##		denom <- (p1*d[, "A"] + p2*d[, "AAAB"] + p3*d[, "AB"] + p4*d[, "ABBB"] + p5*d[, "B"] + p6)
##		mAAAA <- p1*d[,"A"]/denom
##		mAAAB <- p2*d[,"AAAB"]/denom
##		mAABB <- p3*d[,"AB"]/denom
##		mABBB <- p4*d[,"ABBB"]/denom
##		mBBBB <- p5*d[,"B"]/denom
##		m4out <- 1-mAAAA-mAAAB-mAABB-mABBB-mBBBB
##		##.index=subjectHits(findOverlaps(ir[1,], featureData(trioSet)))
##		##h[.index, ]
##		g0 <- h[, 1]*d0
##		g1A <- h[, 2]*mA
##		g1B <- h[, 2]*mB
##		g1out <- h[, 2]*m1out
##		g2AA <- h[, 3]*mAA
##		g2AB <- h[, 3]*mAB
##		g2BB <- h[, 3]*mBB
##		g2out <- h[,3]*m2out
##		g2rohAA <- h[, 4]*mA
##		g2rohBB <- h[, 4]*mB
##		g3AAA <- h[,5]*mAAA
##		g3AAB <- pmax(h[,5]*mAAB,1e-9)
##		g3ABB <- pmax(h[,5]*mABB,1e-9)
##		g3BBB <- h[,5]*mBBB
##		g3out <- h[,5]*m3out
##		g4AAAA <- h[,6]*mAAAA
##		g4AAAB <- pmax(h[,6]*mAABB, 1e-9)
##		g4AABB <- pmax(h[,6]*mAABB, 1e-9)
##		g4ABBB <- pmax(h[,6]*mABBB, 1e-9)
##		g4BBBB <- h[,6]*mBBBB
##		g4out <- h[, 6]*m4out
##		sigma.new <- mu.new <- rep(NA, length(mus))
##		names(sigma.new) <- names(mu.new) <- names(mus)
##		mu.new["A"] <- sum(g1A*x, na.rm=TRUE)/sum(g1A, na.rm=TRUE)
##		mu.new["B"] <- sum(g1B*x, na.rm=TRUE)/sum(g1B, na.rm=TRUE)
##		mu.new["AB"] <- sum(g2AB*x,na.rm=TRUE)/sum(g2AB, na.rm=TRUE)
##		mu.new["AAB"] <- sum(g3AAB*x,na.rm=TRUE)/sum(g3AAB, na.rm=TRUE)
##		mu.new["ABB"] <- sum(g3ABB*x, na.rm=TRUE)/sum(g3ABB, na.rm=TRUE)
##		mu.new["AAAB"] <- sum(g4AAAB*x,na.rm=TRUE)/sum(g4AAAB, na.rm=TRUE)
##		mu.new["ABBB"] <- sum(g4ABBB*x,na.rm=TRUE)/sum(g4ABBB, na.rm=TRUE)
##		sigma.new["A"] <- sqrt(sum(g1A*(x-mu.new["A"])^2, na.rm=TRUE)/sum(g1A, na.rm=TRUE))
##		sigma.new["B"] <- sqrt(sum(g1B*(x-mu.new["B"])^2, na.rm=TRUE)/sum(g1B, na.rm=TRUE))
##		sigma.new["AB"] <- sqrt(sum(g2AB*(x-mu.new["AB"])^2, na.rm=TRUE)/sum(g2AB, na.rm=TRUE))
##		sigma.new["AAB"] <- sqrt(sum(g3AAB*(x-mu.new["AAB"])^2, na.rm=TRUE)/sum(g3AAB, na.rm=TRUE))
##		sigma.new["ABB"] <- sqrt(sum(g3ABB*(x-mu.new["ABB"])^2, na.rm=TRUE)/sum(g3ABB, na.rm=TRUE))
##		sigma.new["AAAB"] <- sqrt(sum(g4AAAB*(x-mu.new["AAAB"])^2, na.rm=TRUE)/sum(g4AAAB, na.rm=TRUE))
##		sigma.new["ABBB"] <- sqrt(sum(g4ABBB*(x-mu.new["ABBB"])^2, na.rm=TRUE)/sum(g4ABBB, na.rm=TRUE))
##		mus <- mu.new
##		mus <- constrainMuBaf(mus)
##		mus <- makeMusBafNondecreasing(mus)
##		sigmas <- sigma.new
##		if(!is.null(constrainSdBaf)){
##			sigmas <- constrainSdBaf(mus=mus, sds=sigmas)
##		}
##		## ensures that .05, 0.95 are at least the .975 or 0.025 quantile
##		sigmas <- pmax(sigmas, 0.025)
##		p1bar <- sum(g1out,na.rm=TRUE)/sum(g1A+g1B+g1out, na.rm=TRUE)
##		e[4, 2] <- e[2, 2] <- p1bar
##		p2bar <- sum(g2out,na.rm=TRUE)/sum(g2AA+g2AB+g2BB+g2out, na.rm=TRUE)
##		e[3, 2] <- p2bar
##		p3bar <- sum(g3out,na.rm=TRUE)/sum(g3AAA+g3AAB+g3ABB+g3BBB+g3out, na.rm=TRUE)
##		e[5, 2] <- p3bar
##		p4bar <- sum(g4out,na.rm=TRUE)/sum(g4AAAA+g4AAAB+g4AABB+g4ABBB+g4BBBB+g4out, na.rm=TRUE)
##		e[6, 2] <- p4bar
##		pr.min <- max(e[, 2], 1e-6)
##		e[, 2] <- pmax(e[, 2], pr.min)
##		e[, 2] <- pmin(e[, 2], pOutlierMax)
##		e[2, 2] <- min(e[2, 2], pOutlierMaxROH)
##		e[, 1] <- 1-e[, 2]
##	}
##	names(mus) <- names(sigmas) <- nms
##	return(list(mu=mus, sigma=sigmas, p=e[, 2]))
##}
##
####constrainSdBaf <- function(mus, sds){
####	sds["AAB"] <- min((0.5 - mus["AAB"])/2, sds["AAB"])
####	sds["ABB"] <- min(abs((0.5 - mus["AAB"]))/2, sds["ABB"])
####	sds["AAAB"] <- min(abs((0.5 - mus["AAAB"]))/3, sds["AAAB"])
####	sds["ABBB"] <- min(abs((0.5 - mus["ABBB"]))/3, sds["ABBB"])
####	return(sds)
####}
##constrainSdBaf <- function(mus, sds) pmax(sds, sds["AB"])
##
##
constrainMuBaf <- function(mus){
	mus["A"] <- ifelse(mus["A"] > 0.001, 0.001, mus["A"])
	mus["B"] <- ifelse(mus["B"] < 0.999, 0.999, mus["B"])
	mus["AAAB"] <- ifelse(mus["AAAB"] > 0.35, 0.35, mus["AAAB"])
	mus["AAAB"] <- ifelse(mus["AAAB"] < 0.05, 0.05, mus["AAAB"])
	mus["AAB"] <- ifelse(mus["AAB"] < 0.1, 0.1, mus["AAB"])
	mus["AAB"] <- ifelse(mus["AAB"] > 0.45, 0.45, mus["AAB"])
	mus["AB"] <- ifelse(mus["AB"] < 0.45, 0.45, mus["AB"])
	mus["AB"] <- ifelse(mus["AB"] > 0.55, 0.55, mus["AB"])
	mus["ABBB"] <- ifelse(mus["ABBB"] > 0.95, 0.95, mus["ABBB"])
	mus["ABBB"] <- ifelse(mus["ABBB"] < 0.65, 0.65, mus["ABBB"])
	mus["ABB"] <- ifelse(mus["ABB"] > 0.9, 0.9, mus["ABB"])
	mus["ABB"] <- ifelse(mus["ABB"] < 0.55, 0.55, mus["ABB"])
##	## symmetry around 0.5
	##delta1 <- mus["AB"] - mus["AAB"]
	delta1 <- 0.5-mus["AAB"]
	delta2 <- mus["ABB"] - 0.5##mus["AB"]
	avg <- (delta1+delta2)/2
	mus["AAB"] <- 0.5-avg
	mus["ABB"] <- 0.5+avg
	delta1 <- 0.5 - mus["AAAB"]
	delta2 <- mus["ABBB"] - 0.5
	avg <- (delta1+delta2)/2
	mus["AAAB"] <- 0.5 - avg
	mus["ABBB"] <- 0.5 + avg
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
	##mus <- mus2[match(nms, names(mus2))]
	mus2
}

isNonDecreasing <- function(mus){
	nms <- names(mus)
	index <- match(c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B"), names(mus))
	all(diff(mus[index]) > 0)
}

