test_hmm_oligoSnpSetWithBAFs <- function(){
	## Results are sensitive to the seed, and differs between the unit test and my
	library(oligoClasses);library(VanillaICE);library(RUnit)
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
	fit <- hmm(oligoset, is.log=FALSE,
		   p.hom=1,
		   cnStates=c(0.5, 1.5, 2, 2, 2.5, 3.2),
		   prOutlierBAF=list(initial=1e-4, max=1e-3, maxROH=1e-5))
	checkEquals(state(fit[[1]]), states)
	checkEquals(coverage2(fit[[1]]), nmarkers, tolerance=0.02)
	if(FALSE){
		library(IRanges)
		library(Biobase)
		rect2 <- function(object){
			object <- object[state(object) !=3 , ]
			object <- object[order(width(object), decreasing=TRUE), ]
			rect(xleft=start(object)/1e6,
			     xright=end(object)/1e6,
			     ybottom=rep(0.7,length(object)),
			     ytop=rep(1,length(object)),
			     col=(1:3)[as.integer(as.factor(state(object)))],
			     border=(1:3)[as.integer(as.factor(state(object)))])
		}
		franges <- makeFeatureGRanges(oligoset)
		olaps <- lapply(fit, findOverlaps, subject=franges)
		plot(position(oligoset)/1e6, copyNumber(oligoset)/100, pch=".", col="black")
		##trace(rect2, browser)
		rect2(fit[[1]])
		rect2(res2[[1]])
		plot(position(oligoset)/1e6, baf(oligoset)/1000, pch=".", col="black")
		rect2(fit)
	}
	res2 <- hmm(object=oligoset, is.log=FALSE, p.hom=0,
		    cnStates=c(0.5, 1, 2, 2, 3, 3.5),
		    prOutlierBAF=list(initial=1e-4, max=1e-3, maxROH=1e-5))
	nmarkers2 <- nmarkers[-c(1,2)]
	nmarkers2[1] <- sum(nmarkers[1:3])
	checkEquals(coverage2(res2[[1]]), nmarkers2, tolerance=0.02)
}

test_hmm_genotypesOnly <- function(){
	library(IRanges)
	states <- as.integer(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3))
	nmarkers <- as.integer(c(996, 102, 902, 50, 2467, 102, 76, 1822,
				 99, 900, 20, 160))
	statepath <- rep(states, nmarkers)
	states2 <- ifelse(states==c(2,4), 2L, 1L)
	states2 <- c(1,2,1,2,1,1,2,1)
	path <- system.file("extdata", package="VanillaICE")
	load(file.path(path, "oligosetForUnitTest.rda"))
	snpSet <- as(oligoset, "SnpSet2")
	##trace(hmmSnpSet, browser)
	res <- hmm(snpSet, S=2L, prGtHom=c(0.7, 0.99), normalIndex=1L)
	checkEquals(state(res[[1]]), states2)

	## test ICE hmm
	## annotation package not supported
	library(Biobase)
	checkException(hmm(snpSet, S=2L, ICE=TRUE, normalIndex=1L))
	annotation(snpSet) <- "genomewidesnp6Crlmm"
	snpCallProbability(snpSet)[c(1049,1050)] <- p2i(0.5)
	hmm(snpSet, S=2L, ICE=TRUE, normalIndex=1L)
}

test_hmm_cnset <- function(){
	library(oligoClasses)
	library(GenomicRanges)
	library(Biobase)
	library2(crlmm)
	data(cnSetExample, package="crlmm")
	if(FALSE){
		cnSetExample@genome <- "hg19"
		save(cnSetExample, file="~/Software/crlmm/data/cnSetExample.rda")
	}
	oligoset2 <- constructOligoSetListFrom(cnSetExample)
	library(foreach);registerDoSEQ()
	res <- hmm(oligoset2, p.hom=0)
	checkEquals(names(res), sampleNames(oligoset2))
	res <- res[[1]]
	rd <- res[state(res)!=3 & coverage2(res) >= 5, ]
	if(FALSE){
		xyplotLrrBaf(rd, oligoset2,
			     frame=200e3,
			     panel=xypanelBaf,
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
		i <- subjectHits(findOverlaps(rd[6, ], oligoset2))
		b <- baf(oligoset2)[i, 2]
		b <- b/1000
		hist(b, breaks=100)
	}
	query <- GRanges("chr8", IRanges(3.7e6,3.8e6))
	j <- subjectHits(findOverlaps(query, rd))
	checkTrue(state(rd)[j]==5)
	checkEquals(coverage2(rd)[j], 776L, tolerance=5)
}

test_oligoSnpSetGT <- function(){
	library(VanillaICE);library(RUnit)
	library(oligoClasses)
	library(GenomicRanges)
	library(Biobase)
	data(oligoSetExample, package="oligoClasses")
	oligoSet <- oligoSet[chromosome(oligoSet) == 1, ]
	hmmResults <- hmm(oligoSet)
	if(FALSE){
		library(IRanges)
		library(Biobase)
		rect2 <- function(object){
			object <- object[order(width(object), decreasing=TRUE), ]
			rect(xleft=start(object)/1e6,
			     xright=end(object)/1e6,
			     ybottom=rep(-1,length(object)),
			     ytop=rep(-0.8,length(object)),
			     col=c("grey0","grey30", "white", "lightblue", "blue", "blue")[state(object)],
			     border=c("grey0","grey30", "white", "lightblue", "blue")[state(object)])
		}
		franges <- makeFeatureGRanges(oligoSet)
		olaps <- lapply(hmmResults, findOverlaps, subject=franges)
		par(las=1)
		plot(position(oligoSet)/1e6, copyNumber(oligoSet)/100, pch=21, cex=0.6,
		     col=c("lightblue", "red", "lightblue")[calls(oligoSet)],
		     bg=c("lightblue", "red", "lightblue")[calls(oligoSet)],
		     xlab="position (Mb)", ylab="log2 copy number")
		rect2(hmmResults[[1]])
	}
	index <- subjectHits(findOverlaps(GRanges("chr1", IRanges(70.1e6, 73e6)), hmmResults[[1]]))
	checkTrue(state(hmmResults[[1]])[index] >= 5)
	options(warn=2)
	## bad starting values, but works out ok
	res2 <- hmm(oligoSet, cnStates=c(-3, 2, 4, 4, 6, 7))
	checkIdentical(state(res2[[1]]), state(hmmResults[[1]]))
}


temp_test_armSplit <- function(){
	object <- bset
	arm <- getArm(object)
	##arm <- oligoClasses:::.getArm(chromosome(object), position(object), genomeBuild(object))
	arm <- factor(arm, unique(arm))
	index <- split(seq_len(nrow(object)), arm)
	l <- elementLengths(index)
	index <- index[l > 10]
}
