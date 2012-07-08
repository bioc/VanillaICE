hmmSnpSet <- function(object,
		      ICE=FALSE,
		      chromosome=1:22,
		      normalIndex=1L,
		      rohIndex=normalIndex+1L,
		      S=2L,
		      prGtHom=c(.95, 0.05),
		      prGtMis=rep(1/S,S),
		      prHetCalledHom=0.001,
		      prHetCalledHet=0.995,
		      prHomInNormal=0.8,
		      prHomInRoh=0.999, TAUP=1e8, ...){ ## additional arguments to Viterbi construction
	object <- object[chromosome(object) %in% chromosome, ]
	marker.index <- validChromosomeIndex(object)
	object <- object[marker.index, ]
	object <- chromosomePositionOrder(object)
	arm <- .getArm(chromosome(object), position(object), genomeBuild(object))
	index <- split(seq_len(nrow(object)), arm)
	chr <- sapply(index, function(i, chr) unique(chr[i]), chr=chromosome(object))
	if(ICE) checkAnnotationForICE(object)
	chrom <- i <- NULL
	rd <- foreach(i=index, chrom=chr) %do% {
		viterbiForSnpSet(gt=calls(object)[i, , drop=FALSE],
				 is.snp=isSnp(object)[i],
				 pos=position(object)[i],
				 gt.conf=confs(object)[i, , drop=FALSE],
				 chrom=chrom,
				 ICE=ICE,
				 annotationPkg=annotation(object),
				 prGtHom=prGtHom,
				 prGtMis=prGtMis,
				 prHetCalledHom=prHetCalledHom,
				 prHetCalledHet=prHetCalledHet,
				 prHomInNormal=prHomInNormal,
				 prHomInRoh=prHomInRoh,
				 TAUP=TAUP,...)
	}
	rd <- stackGRangesList(rd, genomeBuild(object))
	return(rd)
}

viterbiForSnpSet <- function(gt, S=2L,
			     is.snp, pos, gt.conf,
			     chrom,
			     ICE=FALSE,
			     rohIndex=2L,
			     annotationPkg,
			     prGtHom,
			     prGtMis=rep(1/S,S),
			     prHetCalledHom=0.001,
			     prHetCalledHet=0.995,
			     prHomInNormal=0.8,
			     prHomInRoh=0.999,
			     normalIndex=1L,
			     TAUP=1e8,
			     ...){ ## additional arguments passed to Viterbi (e.g., TAUP, normal2altered...)
	log.beta.gt <- gtEmissionFromMatrix(object=gt,
					    S=S,
					    is.snp=is.snp,
					    pos=pos,
					    gt.conf=gt.conf,
					    ICE=ICE,
					    annotationPkg=annotationPkg,
					    prGtHom=prGtHom,
					    prGtMis=prGtMis,
					    prHetCalledHom=prHetCalledHom,
					    prHetCalledHet=prHetCalledHet,
					    prHomInNormal=prHomInNormal,
					    prHomInRoh=prHomInRoh,
					    rohIndex=2L)
	transitionProb <- computeTransitionProb(pos, S=S, TAUP=TAUP)
	## construct ViterbiSet
	J <- ncol(gt)
	viterbiList <- vector("list", J)
	for(j in seq_len(J)){
		tmp <- new("Viterbi",
			   numberFeatures=nrow(gt),
			   numberStates=S,
			   emission=as.matrix(as.numeric(log.beta.gt[, j, ])),
			   normalIndex=normalIndex,
			   transitionProb=transitionProb,
			   ...)
		viterbiList[[j]] <- fit(tmp)
	}
	id <- object <- NULL
	rdlist <- foreach(object = viterbiList,
			  id=colnames(gt)) %do% {
				  computeLoglikFromViterbi2(object=object,
							    chrom=chrom,
							    id=id,
							    pos=pos,
							    log.it=FALSE)
			  }
	res <- GRangesList(rdlist)
	names(res) <- colnames(gt)
	return(res)
}


setMethod("gtEmission", signature(object="SnpSet2"),
	  function(object, hmm.params, ...){
		  is.ordered <- checkOrder(object)
		  stopifnot(is.ordered)
		  log.emit <- gtEmission(calls(object), hmm.params,
					 is.snp=isSnp(object),
					 gt.conf=confs(object),
					 cdfName=annotation(object), ...)
		  return(log.emit)
	  })

setMethod("bafEmission", signature(object="SnpSet2"),
	  function(object, is.snp, cdfName,
		   prOutlier=1e-3, p.hom=0.95, ...){
		  is.ordered <- checkOrder(object)
		  stopifnot(is.ordered)
		  log.emit <- bafEmission(baf(object),
					  is.snp=isSnp(object),
					  cdfName=annotation(object),
					  prOutlier=prOutlier,
					  p.hom=p.hom, ...)
		  return(log.emit)
	  })

setMethod("baf", signature(object="SnpSet2"), function(object) assayDataElement(object, "baf"))

