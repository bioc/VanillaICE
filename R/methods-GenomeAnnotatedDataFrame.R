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
