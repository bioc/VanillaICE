assayDataListDims <- function(object) {
	nms <- if (assayDataStorageMode(object) == "list") names(object) else ls(object)
	if (length(nms) == 0)
		return(matrix(integer(0), nrow = 2, ncol = 0,
			      dimnames = list(c("Features", "Samples"), character(0))))
	d <- lapply(nms, function(i) lapply(object[[i]], dim)) ##dim(object[[i]]))
	names(d) <- nms
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
