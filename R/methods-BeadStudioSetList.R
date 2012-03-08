setMethod("initialize", signature(.Object="BeadStudioSetList"),
	  function(.Object,
		   assayDataList=AssayDataList(BAF=BAF, logRRatio=logRRatio),
		   logRRatio=new("matrix"),
		   BAF=matrix(NA, nrow(logRRatio), ncol(logRRatio)),
		   featureDataList=GenomeAnnotatedDataFrameFromList(assayDataList),
		   chromosome=integer(),
		   phenoData=annotatedDataFrameFrom(assayDataList, byrow=FALSE),
		   annotation=character(),
		   genomeBuild=character(),
		   ...){
		  callNextMethod(.Object,
				 assayDataList=assayDataList,
				 featureDataList=featureDataList,
				 phenoData=phenoData,
				 chromosome=chromosome,
				 annotation=annotation,
				 genomeBuild=genomeBuild,
				 ...)
	  })

setMethod("updateObject", signature(object="BeadStudioSetList"),
          function(object, ..., verbose=FALSE) {
		  if (verbose) message("updateObject(object = 'BeadStudioSetList')")
		  obj <- tryCatch(callNextMethod(object), error=function(e) NULL)
		  if(is.null(obj)){
			  obj <- new("BeadStudioSetList",
				     assayDataList = assayDataList(object),
				     phenoData = phenoData(object),
				     annotation = updateObject(annotation(object),
				     ..., verbose=verbose),
				     featureDataList=featureDataList(object),
				     chromosome=chromosome(object),
				     ...)
		  }
		  obj
          })

initializeLrrAndBafMatrix <- function(nrow, ncol, col.names, outdir){
	ldPath(outdir)
	bafs <- initializeBigMatrix("baf", nr=nrow, nc=ncol, vmode="integer")
	lrrs <- initializeBigMatrix("lrr", nr=nrow, nc=ncol, vmode="integer")
	colnames(bafs) <- colnames(lrrs) <- col.names
	res <- list(baf=bafs, lrr=lrrs)
	return(res)
}

BeadStudioSetList <- function(fnames,
			      annotationPkg,
			      genomeBuild="hg19",
			      outdir=ldPath(),
			      sampleIds,
			      phenoData,
			      ...){
	if(is.null(getCluster())) registerDoSEQ()
	dat <- read.bsfiles(filenames=fnames[1])
	##
	## build featureDataList
	##
	##lrrs <- matrix(dat[,1,], dimnames=list(rownames(dat),NULL))
	fD <- GenomeAnnotatedDataFrameFromArray(dat, annotationPkg, universe=genomeBuild)
	fD <- fD[chromosome(fD) <= 23 & !is.na(chromosome(fD)), ]
	fD <- fD[order(chromosome(fD), position(fD)), ]
	dat <- dat[match(featureNames(fD), rownames(dat)), , , drop=FALSE]
	uchrom <- unique(chromosome(fD))
	featureDataList <- vector("list", length(uchrom))
	marker.index <- split(seq_len(nrow(fD)), chromosome(fD))
	featureDataList <- foreach(i=marker.index) %do% fD[i, ]
	##
	## build assayDataList
	##
	if(isPackageLoaded("ff")){
		message("Initializing ff arrays for BAFs and LRRs")
	} else {
		message("ff package not loaded -- initializing arrays for BAFs and LRRs ...")
	}
	if(missing(sampleIds)) sampleIds <- basename(fnames)
	if(isPackageLoaded("ff")){
		outdir <- ldPath()
		bafAndLrrList <- foreach(i=marker.index, .packages=c("ff", "VanillaICE")) %dopar% {
			VanillaICE:::initializeLrrAndBafMatrix(nrow=length(i),
							       ncol=length(fnames),
							       outdir=outdir,
							       col.names=sampleIds)
		}
		baflist <- lapply(bafAndLrrList, "[[", 1)
		lrrlist <- lapply(bafAndLrrList, "[[", 2)
	} else {
		baflist <- foreach(i=marker.index) %do% matrix(NA, nrow=length(i), ncol=length(sampleIds),
				   dimnames=list(featureNames(fD)[i], sampleIds))
		lrrlist <- foreach(i=marker.index) %do% matrix(NA, nrow=length(i), ncol=length(sampleIds),
				   dimnames=list(featureNames(fD)[i], sampleIds))
	}
	readAndWrite <- function(filename, marker.index, baflist, lrrlist, sampleIds, featureIds){
		for(i in seq_along(filename)){
			j <- match(sampleIds[i], colnames(baflist[[1]]))
			dat <- tryCatch(read.bsfiles(filenames=filename[i]), error=function(e) NULL)
			if(is.null(dat)) {
				warning("some files were not read")
				next()
			}
			i <- match(featureIds, rownames(dat))
			dat <- dat[i, , , drop=FALSE]
			lrrs <- as.integer(dat[, 1, ]*100)
			bafs <- as.integer(dat[, 2, ]*1000)
			rm(dat)
			for(m in seq_along(marker.index)){
				M <- marker.index[[m]]
				baflist[[m]][, j] <- bafs[M]
				lrrlist[[m]][, j] <- lrrs[M]
			}
		}
		if(isPackageLoaded("ff")) return(TRUE)
		return(list(baflist, lrrlist))
	}
	if(isPackageLoaded("ff")){
		message("Reading ", length(fnames),
			" files and writing ff files to directory ", ldPath())
		fileIndices <- splitIndicesByNode(seq_along(fnames))
		res <- foreach(i=fileIndices, .packages=c("ff", "VanillaICE")) %dopar% {
			readAndWrite(filename=fnames[i],
				     marker.index=marker.index,
				     baflist=baflist,
				     lrrlist=lrrlist,
				     sampleIds=sampleIds[i],
				     featureIds=featureNames(fD))
		}
	} else {
		res <- readAndWrite(filename=fnames,
				    marker.index=marker.index,
				    baflist=baflist,
				    lrrlist=lrrlist,
				    sampleIds=sampleIds,
				    featureIds=featureNames(fD))
		baflist <- res[[1]]
		lrrlist <- res[[2]]
		rm(res); gc()
	}
	ad <- AssayDataList(logRRatio=lrrlist, BAF=baflist)
	if(missing(phenoData)) phenoData <- annotatedDataFrameFrom(lrrlist[[1]], byrow=FALSE)
	object <- new("BeadStudioSetList",
		      assayDataList=ad,
		      featureDataList=featureDataList,
		      phenoData=phenoData,
		      chromosome=uchrom,
		      annotation=annotationPkg,
		      genomeBuild=genomeBuild, ...)
	return(object)
}


setReplaceMethod("assayData", signature=signature(object="BeadStudioSetList",
			      value="AssayData"),
                 function(object, value) {
			 object@assayDataList <- value
			 object
                 })

setMethod("[", signature(x="BeadStudioSetList"),
	  function(x, i, j, ..., drop=FALSE){
		  ## using 'i' to subset markers does not really make
		  ## sense
		  ##
		  ## Use i to subset the list. example, x[1] is still a BeadStudioSetList, but is one chromosome
		  ##
		  if(!missing(i) & !missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- ad[[elt]][i]
				  tmp <- lapply(tmp, function(x, j) {
					  x[, j, drop=FALSE]
				  }, j=j)
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@chromosome <- chromosome(x)[i]
			  x@featureDataList <- featureDataList(x)[i]
			  x@phenoData <- phenoData(x)[j, ]
		  }
		  if(!missing(i) & missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- ad[[elt]][i]
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@chromosome <- chromosome(x)[i]
			  x@featureDataList <- featureDataList(x)[i]
		  }
		  if(missing(i) & !missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- lapply(ad[[elt]], function(x, j) {
					  x[, j, drop=FALSE]
				  }, j=j)
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@phenoData <- phenoData(x)[j, ]
		  }
		  return(x)
	  })

setMethod("ncol", signature(x="BeadStudioSetList"),
	  function(x) ncol(x[[1]]))

setMethod("[[", signature(x="BeadStudioSetList"),
	  function(x, i, j, ..., exact=TRUE){
		  if(missing(i)) return(x)
		  if(length(i) == 1){
			  lrrs <- lrr(x)[[i]]
			  bafs <- baf(x)[[i]]
			  fdlist <- featureDataList(x)[[i]]
			  x <- new("BeadStudioSet",
				   lrr=lrrs,
				   baf=bafs,
				   phenoData=phenoData(x),
				   featureData=fdlist,
				   ##is.integer=isInteger(x),
				   genomeBuild=genomeBuild(x),
				   annotation=annotation(x))
		  } else {
			  stop("subscript out of bounds")
		  }
	  })

setMethod("show", signature(object="BeadStudioSetList"),
	  function(object){
		  lo <- length(lrr(object))
		  cat(class(object), " of length ", lo, "\n", sep="")
	  })

setMethod("assayDataList", signature(object="BeadStudioSetList"),
	  function(object)  object@assayDataList)

setMethod("featureDataList", signature(object="BeadStudioSetList"),
	  function(object)  object@featureDataList)

setMethod("featureNames", signature(object="BeadStudioSetList"),
	  function(object){
		  lapply(featureDataList(object), featureNames)
	  })

setMethod("position", signature(object="BeadStudioSetList"),
	  function(object){
		  lapply(featureDataList(object), position)
	  })

setMethod("lrr", signature(object="BeadStudioSetList"),
	  function(object){
		  ##lapply(object, lrr)
		  assayDataList(object)[["logRRatio"]]
	  })

setMethod("baf", signature(object="BeadStudioSetList"),
	  function(object){
		  ##lapply(object, baf)
		  assayDataList(object)[["BAF"]]
	  })

setMethod("chromosome", signature(object="BeadStudioSetList"),
	  function(object, as.list=FALSE, ...){
		  ##lapply(object, chromosome)
		  if(!as.list) object@chromosome else chromosomeList(object)
	  })

setMethod("varLabels", signature(object="BeadStudioSetList"),
	  function(object) varLabels(phenoData(object)))

setMethod("pData", signature(object="BeadStudioSetList"),
	  function(object) pData(phenoData(object)))

setMethod("$", signature(x="BeadStudioSetList"),
	  function(x, name){
		  eval(substitute(phenoData(x)$NAME_ARG, list(NAME_ARG=name)))
	  })
setMethod("order", signature(...="BeadStudioSetList"),
	  function(..., na.last=TRUE,decreasing=FALSE){
		  x <- list(...)[[1]]
		  for(i in seq_along(x)){
			  x[[i]] <- chromosomePositionOrder(x[[i]])
		  }
		  return(x)
	  })

setMethod("annotation", signature(object="BeadStudioSetList"), function(object) object@annotation)
setMethod("genomeBuild", signature(object="BeadStudioSetList"), function(object) object@genomeBuild)

setMethod("dims", signature(object="BeadStudioSetList"), function(object){
	nchr <- length(chromosome(object))
	ntrios <- ncol(baf(object)[[1]])
	dm <- c(nchr, ntrios)
	names(dm) <- c("chromosomes", "trios")
	return(dm)
})

setMethod("sampleNames", signature(object="BeadStudioSetList"),
	  function(object) sampleNames(phenoData(object)))

setMethod("phenoData", signature(object="BeadStudioSetList"),
	  function(object) object@phenoData)

AssayDataList <- function(storage.mode = c("lockedEnvironment", "environment", "list"), ...) {
	storage.mode <- match.arg(storage.mode) ## defaults to "lockedEnvironment"
	assayData <- switch(storage.mode,
			    lockedEnvironment =,
			    environment = new.env(parent=emptyenv()),
			    list = list())
	arglist <- list(...)
	for (nm in names(arglist)) assayData[[nm]] <- arglist[[nm]]
	if (storage.mode == "lockedEnvironment") Biobase:::assayDataEnvLock(assayData)
	assayData
}


setMethod("assayData", signature(object="BeadStudioSetList"),
	  function(object) assayDataList(object))

setMethod("storageMode", "BeadStudioSetList", function(object) storageMode(assayData(object)))

setMethod("length", signature(x="BeadStudioSetList"), function(x) length(x@chromosome))


assayDataListDims <- function(object) {
	nms <- if (assayDataStorageMode(object) == "list") names(object) else ls(object)
	if (length(nms) == 0)
		return(matrix(integer(0), nrow = 2, ncol = 0,
			      dimnames = list(c("Features", "Samples"), character(0))))
	d <- lapply(nms, function(i) lapply(object[[i]], dim)) ##dim(object[[i]]))
	##rownames(d) <- c("Features", "Samples", rep("...", nrow(d)-2))
	names(d) <- nms
	##colnames(d) <- nms
	##d[,order(colnames(d)), drop=FALSE]
	return(d)
}

validAssayDataDims <- function(object){
	msg <- NULL
	d <- assayDataListDims(object)
	firstElement <- d[[1]]
	d <- d[-1]
	res <- sapply(d, function(i) identical(i, firstElement))
	## check that the 3rd dimension is 3
	if(!all(res)){
		msg <- "Assay data elements must have the same dimension"
	}
	if(is.null(msg)) return(TRUE) else return(msg)
}

assayDataStorageMode <- Biobase:::assayDataStorageMode

##setMethod("isInteger", signature(object="BeadStudioSetList"),
##	  function(object) object@is.integer)

setMethod("stack", signature(x="BeadStudioSetList"),
	  function(x, ...){
		  b <- baf(x)
		  Rs <- sapply(b, nrow)
		  Cs <- ncol(b[[1]])
		  logRR <- bf <- matrix(NA, sum(Rs), Cs)
		  chrom <- factor(rep(chromosome(x), Rs), ordered=TRUE, levels=chromosome(x))
		  index <- split(seq_len(sum(Rs)), chrom)
		  for(i in seq_along(x)){
			  j <- index[[i]]
			  bf[j, ] <- baf(x[[i]])[,]
			  logRR[j, ] <- lrr(x[[i]])[,]
		  }
		  fns <- unlist(featureNames(x))
		  dimnames(bf) <- dimnames(logRR) <- list(fns, sampleNames(x[[1]]))
		  pos <- unlist(position(x))
		  issnp <- unlist(lapply(x@featureDataList, isSnp))
		  featureData <- new("GenomeAnnotatedDataFrame",
				     position=pos,
				     chromosome=as.integer(rep(chromosome(x), Rs)),
				     isSnp=issnp,
				     row.names=fns)
		  obj <- new("BeadStudioSet",
			     baf=bf,
			     lrr=logRR,
			     featureData=featureData,
			     phenoData=phenoData(x),
			     annotation=annotation(x),
			     genomeBuild=genomeBuild(x))
			     ##is.integer=isInteger(x))
		  return(obj)
	  })

stackFeatureDataList <- function(x){
	fns <- unlist(lapply(x, featureNames))
	pos <- unlist(lapply(x, position))
	issnp <- unlist(lapply(x, isSnp))
	chrom <- unlist(lapply(x, chromosome))
	featureData <- new("GenomeAnnotatedDataFrame",
			   position=pos,
			   chromosome=chrom,
			   isSnp=issnp,
			   row.names=fns)
}

setReplaceMethod("$", "BeadStudioSetList", function(x, name, value) {
	phenoData(x)[[name]] = value
	x
})

setReplaceMethod("phenoData",
                 signature=signature(
		 object="BeadStudioSetList",
                   value="AnnotatedDataFrame"),
                 function(object, value) {
			 object@phenoData <- value
			 object
                 })

setReplaceMethod("pData",
                 signature=signature(
		 object="BeadStudioSetList",
                   value="data.frame"),
                 function(object, value) {
			 pd <- phenoData(object)
			 pData(pd) <- value
			 phenoData(object) <- pd
			 object
                 })
