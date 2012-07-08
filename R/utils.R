getCrlmmReference <- function(x) paste(annotation(x), "Conf", sep="")
isLoaded <- function(dataset, environ=.vanillaIcePkgEnv)
	exists(dataset, envir=environ)
getVarInEnv <- function(dataset, environ=.vanillaIcePkgEnv){
	if (!isLoaded(dataset))
		stop("Variable ", dataset, " not found in .vanillaIcePkgEnv")
	environ[[dataset]]
}

loader <- function(theFile, envir, pkgname){
	theFile <- file.path(system.file(package=pkgname),
			     "extdata", theFile)
	if (!file.exists(theFile))
		stop("File ", theFile, "does not exist in ", pkgname)
	load(theFile, envir=envir)
}

isAffy <- function(pkg){
	if(pkg=="") {
		return(NULL)
	} else {
		is.affy <- any(c("pd.mapping50k.hind240", "pd.mapping50k.xba240",
				 "pd.mapping50k.hind240,pd.mapping50k.xba240",
				 "pd.mapping250k.nsp",
				 "pd.mapping250k.sty",
				 "pd.mapping250k.nsp,pd.mapping250k.sty",
				 "pd.genomewidesnp.5",
				 "pd.genomewidesnp.6",
				 "genomewidesnp6Crlmm",
				 "genomewidesnp5Crlmm", "gw6crlmm") %in% pkg)
	}
	return(is.affy)
}
