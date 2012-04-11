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
	##trace(VanillaICE:::updateBaf, browser)
	##trace(VanillaICE:::bafEmissionFromMatrix2, browser)
	##
##	arm <- VanillaICE:::.getArm(chromosome(oligoset), position(oligoset))
##	index <- split(seq_along(arm), arm)
##	tmp <- oligoset[index[[2]], ]
##
##	trace(VanillaICE:::updateBaf, browser)
##
##	.index <- subjectHits(findOverlaps(fit[2, ], featureData(tmp)))
##	trace(VanillaICE::viterbi2Wrapper, browser)
##	fit <- hmm(tmp, is.log=FALSE, p.hom=1, cnStates=c(0.5, 1.5, 2, 2, 2.5, 3.2))
##	trace(VanillaICE::viterbi2Wrapper, browser)
	fit <- hmm(oligoset, is.log=FALSE, p.hom=1, cnStates=c(0.5, 1.5, 2, 2, 2.5, 3.2))
##	fit <- hmm(oligoset[, 2], is.log=FALSE, p.hom=1, cnStates=c(0.5, 1.5, 2, 2, 2.5, 3.2))
	checkIdentical(state(fit), states)
	checkEquals(coverage2(fit), nmarkers, tolerance=0.02)
	if(FALSE){
		library(IRanges)
		library(Biobase)
		rect2 <- function(object){
			object <- object[state(object) !=3 , ]
			object <- object[order(width(object), decreasing=TRUE), ]
			rect(xleft=start(object)/1e6,xright=end(object)/1e6,
			     ybottom=rep(0.7,nrow(object)),
			     ytop=rep(1,nrow(object)),
			     col=(1:3)[as.integer(as.factor(state(object)))],
			     border=(1:3)[as.integer(as.factor(state(object)))])
		}
		plot(position(oligoset)/1e6, copyNumber(oligoset)/100, pch=".", col="black")
		rect2(fit)
		o <- subjectHits(findOverlaps(fit[4, ], featureData(oligoset)))
		plot(position(oligoset)/1e6, baf(oligoset)/1000, pch=".", col="black")
		rect2(fit)
	}
	## do not call copy-neutral ROH
##	arm <- VanillaICE:::.getArm(chromosome(oligoset), position(oligoset))
##	index <- split(seq_along(arm), arm)
##	tmp <- oligoset[index[[2]], ]
##
##	trace(VanillaICE:::updateBaf, browser)
##	fit <- hmm(tmp, is.log=FALSE, p.hom=1, cnStates=c(0.5, 1.5, 2, 2, 2.5, 3.2))
##	res2 <- hmm(object=tmp, is.log=FALSE, p.hom=0)
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
##	trace(VanillaICE::viterbi2Wrapper, browser)
	res <- hmm(oligoset, p.hom=0)
	rd <- res[state(res)!=3 & coverage2(res) >= 5, ]
	##baf(oligoset)[.index, ]/1000
	if(FALSE){
		library(Biobase)
		library(IRanges)
		SNPchip:::xyplotLrrBaf(rd, oligoset,
				       frame=200e3,
				       panel=SNPchip:::xypanelBaf,
				       cex=0.5,
				       scales=list(x=list(relation="free"),
				       y=list(alternating=1,
				       at=c(-1, 0, log2(3/2), log2(4/2)),
				       labels=expression(-1, 0, log[2](3/2), log[2](4/2)))),
				       par.strip.text=list(cex=0.7),
				       ylim=c(-3,1),
				       col.hom="grey50",
				       col.het="grey50",
				       col.np="grey20",
				       key=list(text=list(c(expression(log[2]("R ratios")), expression("B allele freqencies")),
						col=c("grey", "blue")), columns=2))
		i <- subjectHits(findOverlaps(rd[6, ], oligoset))
		b <- baf(oligoset)[i, 2]
		b <- b/1000
		hist(b, breaks=100)
	}
	##checkIdentical(state(rd), as.integer(c(5,2,4, 2)))
	##checkEquals(coverage2(rd), as.integer(c(775, 36, 45, 4)), tolerance=5)
	##checkIdentical(state(rd), as.integer(c(5,2)))
	##checkEquals(coverage2(rd), as.integer(c(778,25)), tolerance=5)
	checkIdentical(state(rd), as.integer(c(5,2,5)))
	checkEquals(coverage2(rd), as.integer(c(779,34,36)), tolerance=5)
}
