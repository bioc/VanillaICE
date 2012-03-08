test_Viterbi2 <- function(){
	checkTrue(validObject(new("Viterbi")))
	checkTrue(validObject(new("Viterbi", numberFeatures=100L,
				  numberStates=2L,
				  normalIndex=1L)))
	checkTrue(validObject(new("Viterbi2")))
	checkTrue(validObject(new("Viterbi2", numberFeatures=100L,
				  numberStates=6L)))
}

test_BeadStudioSetList <- function(){
	validObject(new("BeadStudioSetList"))
	path <- system.file("extdata", package="VanillaICE")
	fname <- list.files(path, pattern="LRRand", full.names=TRUE)
	obj1 <- BeadStudioSetList(fnames=fname,
				  ##annotationPkg="gw6crlmm",
				  annotationPkg="genomewidesnp6Crlmm",
				  genomeBuild="")
	checkTrue(validObject(obj1))

	library(ff)
	ldPath(tempdir())
	registerDoSEQ()
	obj2 <- BeadStudioSetList(fnames=fname,
				  ##annotationPkg="gw6crlmm",
				  annotationPkg="genomewidesnp6Crlmm",
				  genomeBuild="")
	checkTrue(validObject(obj2))
	checkTrue(identical(as.numeric(baf(obj1)[[1]]),
			    as.numeric(baf(obj2)[[1]][,1, drop=FALSE])))
}
