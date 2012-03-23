test_hmm_oligoSnpSetWithBAFs <- function(){
	## Results are sensitive to the seed, and differs between the unit test and my
	library(oligoClasses)
	states <- as.integer(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3))
	nmarkers <- as.integer(c(996, 102, 902, 50, 2467, 102, 76, 1822,
				 99, 900, 20, 160))
	statepath <- rep(states, nmarkers)
	if(FALSE){
		oligoset <- VanillaICE:::artificialData(states, nmarkers)
		save(oligoset, file="~/Software/VanillaICE/inst/extdata/oligosetForUnitTest.rda")
	} else {
		path <- system.file("extdata", package="VanillaICE")
		load(file.path(path, "oligosetForUnitTest.rda"))
	}
	## produces an error -- can't find replacement method for baf
	copyNumber(oligoset)[c(5, 6), ] <- NA
	##	baf(oligoset)[c(5,6), ] <- NA
	##because of the centromere, there's an extra normal
	##trace(VanillaICE:::hmmBeadStudioSet, browser)
	fit <- hmm(oligoset, is.log=FALSE, p.hom=1, cnStates=c(0.5, 1.5, 2, 2, 2.5, 3.2))
	checkIdentical(state(fit), states)
	if(FALSE){
		rect2 <- function(object){
			object <- object[state(object) !=3 , ]
			object <- object[order(width(object), decreasing=TRUE), ]
			rect(xleft=start(object)/1e6,xright=end(object)/1e6,
			     ybottom=rep(0.7,nrow(object)),
			     ytop=rep(1,nrow(object)),
			     col=(1:3)[as.integer(as.factor(state(object)))],
			     border=(1:3)[as.integer(as.factor(state(object)))])
		}
		plot(position(oligoset)/1e6, copyNumber(oligoset)/100, pch=".", col="grey")
		rect2(fit)
	}
	checkEquals(coverage2(fit), nmarkers, tolerance=0.02)
	## do not call copy-neutral ROH
	res2 <- hmm(object=oligoset, is.log=FALSE, p.hom=0)
	checkIdentical(state(res2), states[-c(1,2)])
	nmarkers2 <- nmarkers[-c(1,2)]
	nmarkers2[1] <- sum(nmarkers[1:3])
	checkEquals(coverage2(res2), nmarkers2, tolerance=0.02)
}

test_hmm_genotypesOnly <- function(){
	states <- as.integer(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3))
	nmarkers <- as.integer(c(996, 102, 902, 50, 2467, 102, 76, 1822,
				 99, 900, 20, 160))
	statepath <- rep(states, nmarkers)
	states2 <- ifelse(states==c(2,4), 2L, 1L)
	states2 <- c(1,2,1,2,1,1,2,1)
	path <- system.file("extdata", package="VanillaICE")
	load(file.path(path, "oligosetForUnitTest.rda"))
	snpSet <- as(oligoset, "SnpSet")
	res <- hmm(snpSet, S=2L, prGtHom=c(0.7, 0.99), normalIndex=1L)

	checkIdentical(state(res), as.integer(states2))

	## test ICE hmm
	## annotation package not supported
	library(Biobase)
	checkException(hmm(snpSet, S=2L, ICE=TRUE, normalIndex=1L))
	annotation(snpSet) <- "genomewidesnp6Crlmm"
	snpCallProbability(snpSet)[c(1049,1050)] <- p2i(0.5)
	res2 <- hmm(snpSet, S=2L, ICE=TRUE, normalIndex=1L)
}

test_hmm_cnset <- function(){
	library(oligoClasses)
	library2(crlmm)
	data(cnSetExample, package="crlmm")
	oligoset <- as(cnSetExample, "oligoSnpSet")
	## this object is not ordered by physical position
	## make sure the right answer is returned even though
	## its not ordered
	oligoset <- chromosomePositionOrder(oligoset)
	##trace(hmmBeadStudioSet, browser)
	res <- hmm(oligoset, p.hom=0)
	rd <- res[state(res)!=3, ]
	if(FALSE){
		SNPchip:::xyplotLrrBaf(rd, oligoset,
				       frame=200e3,
				       panel=SNPchip:::xypanelBaf,
				       scales=list(x="free"))
				       ##cex=0.2)

		i <- subjectHits(findOverlaps(rd[6, ], oligoset))
		b <- baf(oligoset)[i, 2]
		b <- b/1000
		hist(b, breaks=100)
	}
	checkIdentical(state(rd), as.integer(c(5,2,4, 2)))
	checkEquals(coverage2(rd), as.integer(c(775, 36, 45, 4)))
}
