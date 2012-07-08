setAs("CNSet", "oligoSetList", function(from, to){
	constructOligoSetListFrom(from)
})

constructOligoSetListFrom <- function(object, ...){
	##row.index <- seq_len(nrow(object))
	##col.index <- seq_len(ncol(object))
	is.lds <- ifelse(is(calls(object), "ff_matrix") | is(calls(object), "ffdf"), TRUE, FALSE)
	if(is.lds) stopifnot(isPackageLoaded("ff"))
	b.r <- calculateRBaf(object, ...)
	b <- b.r[["baf"]]
	r <- b.r[["lrr"]]
	j <- match(colnames(r[[1]]), sampleNames(object))
	rns <- lapply(r, rownames)
	fDList <- foreach(featureid=rns) %do%{
		featureData(object)[match(featureid, featureNames(object)), ]
	}
	names(fDList) <- sapply(fDList, function(x) chromosome(x)[1])
	gtPlist <- gtlist <- vector("list", length(r))
	for(i in seq_along(r)){
		gtlist[[i]] <- initializeBigMatrix("call", nr=nrow(r[[i]]), nc=length(j), vmode="integer")
		gtPlist[[i]] <- initializeBigMatrix("callPr", nr=nrow(r[[i]]), nc=length(j), vmode="integer")
		featureid <- rownames(r[[i]])
		ix <- match(featureid, featureNames(object))
		rownames(gtPlist[[i]]) <- rownames(gtlist[[i]]) <- featureid
		colnames(gtPlist[[i]]) <- colnames(gtlist[[i]]) <- colnames(r[[i]])
		for(k in seq_along(j)){
			gtlist[[i]][, k] <- calls(object)[ix, j[k]]
			gtPlist[[i]][, k] <- snpCallProbability(object)[ix, j[k]]
		}
	}
	ad <- AssayDataList(baf=b, copyNumber=r, call=gtlist, callProbability=gtPlist)
	object <- new("oligoSetList",
		      assayDataList=ad,
		      featureDataList=fDList,
		      chromosome=names(fDList),
		      phenoData=phenoData(object)[j, ],
		      annotation=annotation(object),
		      genome=genomeBuild(object))
	return(object)
}

constructBafLrrSetListFrom <- function(object, ...){
	is.lds <- ifelse(is(calls(object), "ff_matrix") | is(calls(object), "ffdf"), TRUE, FALSE)
	if(is.lds) stopifnot(isPackageLoaded("ff"))
	b.r <- calculateRBaf(object, ...)
	b <- b.r[["baf"]]
	r <- b.r[["lrr"]]
	j <- match(colnames(r[[1]]), sampleNames(object))
	rns <- lapply(r, rownames)
	featureid <- NULL
	fDList <- foreach(featureid=rns) %do%{
		featureData(object)[match(featureid, featureNames(object)), ]
	}
	names(fDList) <- sapply(fDList, function(x) chromosome(x)[1])
	ad <- AssayDataList(baf=b, lrr=r)
	object <- new("BafLrrSetList",
		      assayDataList=ad,
		      featureDataList=fDList,
		      chromosome=names(fDList),
		      phenoData=phenoData(object)[j, ],
		      annotation=annotation(object),
		      genome=genomeBuild(object))
	return(object)
}

hmmCNSet <- function(object, ...){
	object <- as(object, "oligoSnpSet")
	hmmOligoSnpSet(object, ...)
}

