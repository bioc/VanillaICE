simulateSingleDupBaf <- function(b, is.snp, from, to, ...){
	stopifnot(is(b, "numeric"))
	names(b) <- NULL
	b.all <- b
	b <- b[is.snp]
	sds <- VanillaICE:::updateSigma(b, is.snp=rep(TRUE, length(b)), sigma0=c(0.02, 0.04, 0.02))

	p <- 0.5
	q3 <- (1-p)^3
	p2q <- p^2*(1-p)
	pq2 <- p*(1-p)^2
	p3 <- p^3
	##tmp <- cbind(q3, 3*pq2, 3*p2q, p3)
	##stopifnot(all(rowSums(tmp) == 1))
	index <- seq(from, to, by=1)

	z <- sample(1:4, size=length(index), replace=TRUE, prob=c(q3, 3*pq2, 3*p2q, p3))
	nZ <- table(z)
	d1 <- rtnorm(length(index), 0, sds[1], lower=0, upper=1)
	d2 <- rtnorm(length(index), 1/3, sds[1], lower=0, upper=1)
	d3 <- rtnorm(length(index), 2/3, sds[1], lower=0, upper=1)
	d4 <- rtnorm(length(index), 1, sds[3], lower=0, upper=1)
	simB <- rep(NA, length(index))
	simB[z==1] <- sample(d1, nZ[1])
	simB[z==2] <- sample(d2, nZ[2])
	simB[z==3] <- sample(d3, nZ[3])
	simB[z==4] <- sample(d4, nZ[4])
	b.all[index] <- simB
	##het <- rtnorm(length(index), 0.5, sds[2], lower=0, upper=1)
	return(b.all)
}

simulateDoubleDupBaf <- function(b, is.snp, from, to, ...){
	stopifnot(is(b, "numeric"))
	names(b) <- NULL
	b.all <- b
	b <- b[is.snp]
	sds <- VanillaICE:::updateSigma(b, is.snp=rep(TRUE, length(b)), sigma0=c(0.02, 0.04, 0.02))
	p <- 0.5
	q4 <- (1-p)^4
	pq3 <- p*(1-p)^3
	p2q2 <- p^2*(1-p)^2
	p3q <- p^3*(1-p)
	p4 <- p^4
	##tmp <- cbind(q3, 3*pq2, 3*p2q, p3)
	##stopifnot(all(rowSums(tmp) == 1))
	index <- seq(from, to, by=1)
	z <- sample(1:5, size=length(index), replace=TRUE, prob=c(q4, pq3, p2q2, p3q, p4))
	nZ <- table(z)
	d1 <- rtnorm(length(index), 0, sds[1], lower=0, upper=1)
	d2 <- rtnorm(length(index), 1/4, sds[1], lower=0, upper=1)
	d3 <- rtnorm(length(index), 1/2, sds[2], lower=0, upper=1)
	d4 <- rtnorm(length(index), 3/4, sds[1], lower=0, upper=1)
	d5 <- rtnorm(length(index), 1, sds[3], lower=0, upper=1)
	simB <- rep(NA, length(index))
	simB[z==1] <- sample(d1, nZ[1])
	simB[z==2] <- sample(d2, nZ[2])
	simB[z==3] <- sample(d3, nZ[3])
	simB[z==4] <- sample(d4, nZ[4])
	simB[z==5] <- sample(d5, nZ[5])
	b.all[index] <- simB
	return(b.all)
}

simulateSingleDelBaf <- function(b, is.snp, from, to, ...){
	index <- seq(from, to, by=1)
	stopifnot(all(diff(index) > 0))
	stopifnot(length(index) > 1)
	stopifnot(is(b, "numeric"))
	names(b) <- NULL
	b.all <- b
	b <- b[is.snp]
	sds <- VanillaICE:::updateSigma(b, is.snp=rep(TRUE, length(b)), sigma0=c(0.02, 0.04, 0.02))
	p <- 0.5
	q <- 1-p

	z <- sample(1:2, size=length(index), replace=TRUE, prob=c(p, q))
	nZ <- table(z)
	d1 <- rtnorm(length(index), 0, sds[1], lower=0, upper=1)
	d2 <- rtnorm(length(index), 1, sds[3], lower=0, upper=1)
	simB <- rep(NA, length(index))
	simB[z==1] <- sample(d1, nZ[1])
	simB[z==2] <- sample(d2, nZ[2])
	b.all[index] <- simB
	##het <- rtnorm(length(index), 0.5, sds[2], lower=0, upper=1)
	return(b.all)
}

artificialData <- function(states, nmarkers){
	data(oligoSetExample, package="oligoClasses")
	oligoSet <- chromosomePositionOrder(oligoSet)
	pos <- position(oligoSet)
	state.path <- rep(states, nmarkers)
	copynumber <- rep(2, length(state.path))
	copynumber[state.path==2] <- 1.5##bias
	copynumber[state.path==5] <- 2.5##bias
	genotypes <- rep(NA, length(copynumber))
	gt <- rmultinom(n=length(copynumber), size=1, prob=rep(1/3,3))
	genotypes[gt[1, ] == 1] <- 1L
	genotypes[gt[2, ] == 1] <- 2L
	genotypes[gt[3, ] == 1] <- 3L
	genotypes[state.path==4 | state.path==2] <- 1L
	genotypes <- as.matrix(genotypes)
	## make signal fairly obvious
	sigmas <- rgamma(length(copynumber), 4, scale=0.05)
	b <- rbaf(as.matrix(genotypes), sigma=0.01, epsilon=0.001, states=state.path)
	dat <- as.matrix(rnorm(length(state.path), mean=copynumber, sd=sigmas))
	i <- seq_along(state.path)
	rownames(dat) <- rownames(genotypes) <- featureNames(oligoSet)[i]
	##pos <- seq(1, by=3e3, length.out=length(copynumber))
	object <- new("oligoSnpSet",
		      copyNumber=integerMatrix(dat, scale=100),
		      call=as.matrix(genotypes),
		      callProbability=snpCallProbability(oligoSet)[i, , drop=FALSE])
	##baf(object) <- b
	assayDataElement(object, "baf") <- integerMatrix(b, scale=1000)
	##df <- data.frame(position=pos, chromosome=rep(1L, length(pos)), isSnp=
	fData(object)$position <- as.integer(pos[i])
	fData(object)$chromosome <- 1L
	fData(object)$isSnp <- TRUE
	return(object)
}

rbaf <- function(genotypes, sigma, epsilon, states){
	baf <- matrix(NA, nrow(genotypes), ncol(genotypes))
	Ns <- table(genotypes)
	a <- pnorm(0, mean=0, sd=sigma)
	b <- pnorm(1, mean=0, sd=sigma)
	I <- runif(Ns[1], 0, 1) > epsilon
	baf[genotypes==1] <- I*qnorm(a+runif(Ns[1], 0, b-a), mean=0, sd=sigma) + (1-I)*runif(Ns[1], 0, 1)
	I <- runif(Ns[2], 0, 1) > epsilon
	baf[genotypes==2] <- I*rnorm(Ns[2], mean=0.5, sd=sigma*2) + (1-I)*runif(Ns[2], 0, 1)
	a <- pnorm(0, mean=1, sd=sigma)
	b <- pnorm(1, mean=1, sd=sigma)
	I <- runif(Ns[3], 0, 1) > epsilon
	baf[genotypes==3] <- I*qnorm(a+runif(Ns[3], 0, b-a), mean=1, sd=sigma) + (1-I) * runif(Ns[3], 0, 1)

	## assume 1/2 are hets to make it easy
	ndup <- sum(states==5)
	ndup.het <- ceiling(ndup/2)
	ndup.hom <- ndup-ndup.het

	index5 <- which(states==5)
	index25 <- sample(index5, ceiling(ndup.het/2))
	index5 <- setdiff(index5, index25)

	index75 <- setdiff(index5, floor(ndup.het/2))
	indexhom <- setdiff(index5, index75)

	n25 <- length(index25)
	n75 <- length(index75)
	nhom <- length(index5)
	I <- runif(n25, 0, 1) > epsilon
	baf[index25] <- I*rnorm(n25, mean=1/3, sd=sigma*2) + (1-I)*runif(n25, 0, 1)
	I <- runif(n75, 0, 1) > epsilon
	baf[index75] <- I*rnorm(n75, mean=2/3, sd=sigma*2) + (1-I)*runif(n75, 0, 1)
	a <- pnorm(0, mean=1, sd=sigma)
	b <- pnorm(1, mean=1, sd=sigma)
	baf[indexhom] <- qnorm(a+runif(Ns[3], 0, b-a), mean=1, sd=sigma)
	rownames(baf) <- rownames(genotypes)
	return(baf)
}
