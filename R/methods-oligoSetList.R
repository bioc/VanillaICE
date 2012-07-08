setMethod("[[", signature(x="oligoSetList"),
	  function(x, i, j, ..., exact=TRUE){
		  if(missing(i)) return(x)
		  if(length(i) == 1){
			  r <- copyNumber(x)[[i]]
			  b <- baf(x)[[i]]
			  gt <- calls(x)[[i]]
			  gtp <- snpCallProbability(x)[[i]]
			  colnames(r) <- colnames(b) <- colnames(gt) <- colnames(gtp) <- sampleNames(phenoData(x))
			  fd <- featureDataList(x)[[i]]
			  x <- new("oligoSnpSet",
				   call=gt,
				   callProbability=gtp,
				   copyNumber=r,
				   baf=b,
				   phenoData=phenoData(x),
				   featureData=fd,
				   genome=genomeBuild(x),
				   annotation=annotation(x))
		  } else {
			  stop("subscript out of bounds")
		  }
	  })

hmmOligoSetList <- function(object, sampleIds, ...){
	message("Fitting HMM to each chromosome")
	if(isPackageLoaded("ff")) pkgs <- c("ff", "Biobase", "VanillaICE") else pkgs <- c("Biobase", "VanillaICE")
	if(missing(sampleIds)) sampleIds <- sampleNames(object)
	if(ncol(object[[1]]) < ocSamples()){
		rdl <- foreach(obj=object, .packages=pkgs) %dopar% {
			hmmOligoSnpSet(object=obj, sampleIds=sampleIds, ...)
		}
		rd <- stackRangedData(rdl)
	} else {
		sample.index <- splitIndicesByLength(seq_along(sampleIds), ocSamples())
		## create a single stream of tasks that can be executed in parallel
		obj <- NULL; j <- NULL
		rdl <- foreach(j=sample.index, .packages=pkgs) %:%
			foreach(obj=object) %dopar%{
				hmmOligoSnpSet(object=obj, sampleIds=sampleNames(object)[j], ...)
			}
		rdl2 <- foreach(ranges=rdl) %do% stack(RangedDataList(ranges))
		rdl2 <- foreach(ranges=rdl2) %do% ranges[, -ncol(ranges)]
		rdl3 <- stack(RangedDataList(rdl2))
		rd <- rdl3[, -ncol(rdl3)]
	}
	return(rd)
}
