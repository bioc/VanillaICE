THISPKG <- "VanillaICE"

.First.lib <- function(libname, pkgname){
	library.dynam(THISPKG, pkgname, libname)
}

##.onAttach <- function(libname, pkgname) {
##	message("Welcome to VanillaICE version ", packageDescription(THISPKG, field="Version"))
##}

.onUnload <- function(libpath){
	library.dynam.unload(THISPKG, libpath)
}

.onLoad <- function(libname, pkgname){
	require("methods")
	.vanillaIcePkgEnv <- new.env(parent=emptyenv())	
}


