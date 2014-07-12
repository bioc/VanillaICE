THISPKG <- "VanillaICE"

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to VanillaICE version ", packageDescription(THISPKG, fields="Version"))
}

.onUnload <- function(libpath){
  library.dynam.unload(THISPKG, libpath)
}

.vanillaIcePkgEnv <- new.env(parent=emptyenv())
