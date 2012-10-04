setMethod("initialize", signature(.Object="BeadStudioSetList"),
	  function(.Object,
		   assayDataList=AssayDataList(baf=baf, lrr=lrr),
		   lrr=list(),
		   baf=lapply(lrr, function(x) matrix(nrow=nrow(x), ncol=ncol(x))),
		   featureDataList=GenomeAnnotatedDataFrameFrom(assayDataList, annotation, genome),
		   chromosome=vector("list", length(lrr)),
		   phenoData,
		   annotation=character(),
		   genome=character(),
		   ...){
		  if(missing(phenoData)){
			  if(length(lrr) > 0){
				  phenoData <- annotatedDataFrameFrom(lrr[[1]], byrow=FALSE)
			  } else {
				  phenoData <- new("AnnotatedDataFrame")
			  }
		  }
		  callNextMethod(.Object,
				 assayDataList=assayDataList,
				 featureDataList=featureDataList,
				 phenoData=phenoData,
				 chromosome=chromosome,
				 annotation=annotation,
				 genome=genome,
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
				     genome=genomeBuild(object),
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
			      genome=c("hg19","hg18"),
			      outdir=ldPath(),
			      sampleIds,
			      phenoData,
			      byArm=FALSE,
			      centromere,
			      ...){
	if(!getDoParRegistered()) registerDoSEQ()
	if(!missing(sampleIds)) if(length(sampleIds) != length(fnames)) stop("sampleIds and fnames not the same length")
	genome <- match.arg(genome)
	dat <- read.bsfiles(filenames=fnames[1])
	##fD <- GenomeAnnotatedDataFrameFromArray(dat, annotationPkg, genome=genome)
	fD <- GenomeAnnotatedDataFrameFrom(dat, annotationPkg, genome=genome)
	fD <- fD[chromosome(fD) <= 23 & !is.na(chromosome(fD)), ]
	fD <- fD[order(chromosome(fD), position(fD)), ]
	dat <- dat[match(featureNames(fD), rownames(dat)), , , drop=FALSE]
	uchrom <- unique(chromosome(fD))
	if(byArm){
		if(missing(centromere)) stop("byArm is TRUE but centromere coordinates not provided")
		marker.index <- split(seq_len(nrow(fD)), chromosome(fD))
		poslist <- lapply(marker.index, function(i, object) position(object)[i], object=fD)
		assignArm <- function(x, chrom, centromere){
			index <- match(paste("chr", chrom, sep=""), centromere$chrom)
			start.centromere <- centromere[index, "chromStart"]
			ifelse(x <= start.centromere, "p", "q")
		}
		arm <- foreach(i=seq_along(poslist)) %do% {
			assignArm(x=poslist[[i]], chrom=names(poslist)[i], centromere=centromere)
		}
		arm <- unlist(arm)
		arm <- as.integer(factor(arm, levels=c("p", "q")))
		fD <- new("GenomeAnnotatedDataFrame",
			  position=position(fD),
			  chromosome=chromosome(fD),
			  isSnp=isSnp(fD),
			  arm=arm,
			  row.names=featureNames(fD))
		featureNames(fD) <- featureNames(fD)
		chrarm <- paste(chromosome(fD), fD$arm, sep="_")
		levels <- unique(chrarm)
		chrarm <- factor(chrarm, levels=levels)
		marker.index <- split(seq_len(nrow(fD)), chrarm)
		minlength <- 1000
		L <- sapply(marker.index, length)
		if(any(L < minlength)){
			droplevels <- levels[L < minlength]
			newlevel <- strsplit(droplevels, "_")[[1]][[1]]
			armdrop <- strsplit(droplevels, "_")[[1]][[2]]
			i <- grep(droplevels, chrarm)
			chrarm <- as.character(chrarm)
			chrarm[i] <- newlevel
			otherlevel <- paste(newlevel, (1:2)[-as.integer(armdrop)], sep="_")
			chrarm[chrarm==otherlevel] <- newlevel
			levels <- unique(chrarm)
			chrarm <- factor(chrarm, levels=levels)
			marker.index <- split(seq_len(nrow(fD)), chrarm)
		}
		featureDataList <- foreach(i=marker.index) %do% fD[i, ]
	} else {
		##featureDataList <- vector("list", length(uchrom))
		##levels <- factor(chromosome(fD), unique(chromosome(fD)))
		levels <- factor(uchrom)
		marker.index <- split(seq_len(nrow(fD)), chromosome(fD))
		featureDataList <- foreach(i=marker.index) %do% fD[i, ]
	}
	names(featureDataList) <- as.character(levels)
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
		i <- NULL
		bafAndLrrList <- foreach(i=marker.index, .packages=c("oligoClasses", "ff", "VanillaICE")) %dopar% {
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
	for(i in seq_along(baflist)){
		fns <- featureNames(featureDataList[[i]])
		rownames(baflist[[i]]) <- fns
		rownames(lrrlist[[i]]) <- fns
	}
	readAndWrite <- function(filename, marker.index, baflist, lrrlist, sampleIds, featureIds){
		for(i in seq_along(filename)){
			j <- match(sampleIds[i], colnames(baflist[[1]]))
			dat <- tryCatch(read.bsfiles(filenames=filename[i]), error=function(e) NULL)
			if(is.null(dat)) {
				warning("some files were not read")
				next()
			}
			ii <- match(featureIds, rownames(dat))
			dat <- dat[ii, , , drop=FALSE]
			lrrs <- as.integer(dat[, 1, ]*100)
			bafs <- as.integer(dat[, 2, ]*1000)
			rm(dat)
			for(m in seq_along(marker.index)){
				M <- marker.index[[m]]
				baflist[[m]][, j] <- bafs[M]
				lrrlist[[m]][, j] <- lrrs[M]
			}
		}
		if(isPackageLoaded("ff")) TRUE
		return(list(baflist, lrrlist))
	}
	if(isPackageLoaded("ff")){
		message("Reading ", length(fnames),
			" files and writing ff files to directory ", ldPath())
		fileIndices <- splitIndicesByNode(seq_along(fnames))
		res <- foreach(i=fileIndices, .packages=c("oligoClasses", "Biobase", "ff", "VanillaICE")) %dopar% {
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
	}
	rm(res)
	ad <- AssayDataList(lrr=lrrlist, baf=baflist)
	if(missing(phenoData)) phenoData <- annotatedDataFrameFrom(lrrlist[[1]][,,drop=FALSE], byrow=FALSE)
	object <- new("BeadStudioSetList",
		      assayDataList=ad,
		      featureDataList=featureDataList,
		      phenoData=phenoData,
		      chromosome=names(featureDataList),
		      annotation=annotationPkg,
		      genome=genome, ...)
	return(object)
}


setReplaceMethod("assayData", signature=signature(object="BeadStudioSetList",
			      value="AssayData"),
                 function(object, value) {
			 object@assayDataList <- value
			 object
                 })

setMethod("ncol", signature(x="BeadStudioSetList"),
	  function(x) ncol(x[[1]]))

setMethod("[[", signature(x="BeadStudioSetList"),
	  function(x, i, j, ..., exact=TRUE){
		  if(missing(i)) return(x)
		  ad <- assayDataList(x)
		  fdlist <- featureData(x)
		  adnew <- switch(storage.mode(ad),
			  lockedEnvironment =,
				  environment = new.env(parent=emptyenv()),
				  list = list())
		  nms <- ls(ad)
		  if(length(i) == 1){
			  for (nm in ls(ad)){
				  elt <- ad[[nm]][[i]]
				  dimnames(elt) <- lapply(dimnames(elt), unname)
				  adnew[[nm]] <- elt
			  }
		  }
		  x <- new("BeadStudioSet",
			   assayData=adnew,
			   phenoData=phenoData(x),
			   featureData=fdlist[[i]],
			   genome=genomeBuild(x),
			   annotation=annotation(x))
	  })

setMethod("[[", signature(x="BafLrrSetList"),
	  function(x, i, j, ..., exact=TRUE){
		 x <- callNextMethod()
		 new("BafLrrSet",
		     assayData=assayData(x),
		     phenoData=phenoData(x),
		     featureData=featureData(x),
		     genome=genomeBuild(x),
		     annotation=annotation(x))
	  })

setMethod("[", signature(x="gSetList"),
	  function(x, i, j, ..., drop=TRUE){
		  ##if(!missing(i)) return(x)
		  if(missing(i) && missing(j)) return(x)
		  ad <- assayDataList(x)
		  if(!missing(i)){
			  fdlist <- featureData(x)[i]
		  }
		  adnew <- switch(storage.mode(ad),
			  lockedEnvironment =,
				  environment = new.env(parent=emptyenv()),
				  list = list())
		  nms <- ls(ad)
		  if(!missing(i)){
			  for (nm in ls(ad)){
				  elt <- ad[[nm]][i]
				  adnew[[nm]] <- elt
			  }
		  }
		  if(missing(j)){
			  x@featureDataList <- fdlist
			  x@chromosome <- x@chromosome[i]
		  } else {
			  for (nm in ls(ad)){
				  elt <- lapply(ad[[nm]], function(y, j) y[, j, drop=FALSE], j=j)
				  adnew[[nm]] <- elt
			  }
			  phenoData(x) <- phenoData(x)[j, ]
			  x@protocolData <- x@protocolData[j, ]
		  }
		  x@assayDataList <- adnew
		  return(x)
	  })

setReplaceMethod("[[", signature(x="BafLrrSetList", value="BafLrrSet"),
		 function(x, i, j, ..., value){
			 fdl <- x@featureDataList
			 fdl[[i]] <- featureData(value)
			 adl <- x@assayDataList
			 r <- adl[["lrr"]]
			 r[[i]] <- lrr(value)
			 b <- adl[["baf"]]
			 b[[i]] <- baf(value)
			 adl <- AssayDataList(lrr=r, baf=b)
			 new("BafLrrSetList",
			     assayDataList=adl,
			     featureDataList=fdl,
			     phenoData=phenoData(x),
			     chromosome=chromosome(x),
			     annotation=annotation(x),
			     genome=genomeBuild(x))
		 })

setMethod("copyNumber", signature(object="oligoSetList"),
	  function(object) assayData(object)[["copyNumber"]])
setMethod("calls", signature(object="oligoSetList"),
	  function(object) assayData(object)[["call"]])
setMethod("snpCallProbability", signature(object="oligoSetList"),
	  function(object) assayData(object)[["callProbability"]])
setMethod("baf", signature(object="oligoSetList"),
	  function(object) assayData(object)[["baf"]])

setMethod("lrr", signature(object="BeadStudioSetList"),
	  function(object){
		  ##lapply(object, lrr)
		  assayDataList(object)[["lrr"]]
	  })

setMethod("baf", signature(object="BeadStudioSetList"),
	  function(object){
		  ##lapply(object, baf)
		  assayDataList(object)[["baf"]]
	  })

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
		  fns <- as.character(unlist(Biobase::featureNames(x)))
		  dimnames(bf) <- dimnames(logRR) <- list(fns, sampleNames(x[[1]]))
		  pos <- unlist(position(x))
		  issnp <- unlist(lapply(x@featureDataList, isSnp))
		  featureData <- new("GenomeAnnotatedDataFrame",
				     position=pos,
				     chromosome=as.integer(rep(chromosome(x), Rs)),
				     isSnp=issnp,
				     row.names=fns)
		  featureData <- featureData[match(fns, featureNames(featureData)), ]
		  ##all.equal(featureNames(featureData), rownames(bf))
		  obj <- new("BeadStudioSet",
			     baf=bf,
			     lrr=logRR,
			     featureData=featureData,
			     phenoData=phenoData(x),
			     annotation=annotation(x),
			     genome=genomeBuild(x))
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


