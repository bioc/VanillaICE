setMethod("GenomeAnnotatedDataFrameFrom", signature(object="array"),
	  function(object, annotationPkg, ...){
		  GenomeAnnotatedDataFrameFromArray(object, annotationPkg, ...)
	  })

GenomeAnnotatedDataFrameFromArray <- function(object, annotationPkg, ...){
	## coerce to matrix
	dims <- dim(object)
	is.array <- length(dims) == 3
	if(is.array){
		res <- oligoClasses:::GenomeAnnotatedDataFrameFromMatrix(object[, , 1], annotationPkg, ...)
	} else {
		##dim(object) <- dim(object)[c(1,2)]
		res <- oligoClasses:::GenomeAnnotatedDataFrameFromMatrix(object, annotationPkg, ...)
	}
	res
}

GenomeAnnotatedDataFrameFromList <- function(object, annotationPkg){
	nms <- ls(object)
	elt <- object[[nms[1]]]
	fdlist <- vector("list", length(elt))
	for(i in seq_along(elt)){
		fdlist[[i]] <- GenomeAnnotatedDataFrameFrom(elt[[i]], annotationPkg)
	}
	return(fdlist)
}
