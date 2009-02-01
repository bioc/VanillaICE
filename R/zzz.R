THISPKG <- "VanillaICE"

.First.lib <- function(libname, pkgname){
	library.dynam(THISPKG, pkgname, libname)
}

.vanillaIcePkgEnv <- new.env(parent=emptyenv())

.onAttach <- function(libname, pkgname) {
	message("Welcome to VanillaICE version ", packageDescription(THISPKG, field="Version"))
}

.onUnload <- function(libpath){
	library.dynam.unload(THISPKG, libpath)
}
