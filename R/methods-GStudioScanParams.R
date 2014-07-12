#' @export
GStudioScanParams <- function(cnvar="Log R Ratio",
                              bafvar="B Allele Freq",
                              gtvar=c("Allele1 - AB", "Allele2 - AB"),
                              index_genome=integer(),
                              select=integer(),
                              scale=1000){
  new("GStudioScanParams", cnvar=cnvar, bafvar=bafvar, gtvar=gtvar,
      index_genome=index_genome, scale=scale, select=select)
}

setValidity("GStudioScanParams", function(object){
  msg <- TRUE
  if(any(is.na(indexGenome(object)))) {
    msg <- "indexGenome values can not contain NAs"
    return(msg)
  }
  msg
})

setMethod("scale", "GStudioScanParams", function(x, center=TRUE, scale=TRUE) x@scale)
setMethod("indexGenome", "GStudioScanParams", function(object) object@index_genome)
setMethod("selectCols", "GStudioScanParams", function(object) object@select)

setMethod("show", "GStudioScanParams", function(object){
  cat("'GStudioScanParams' class:\n")
  cat("  o columns to select:", paste0(selectCols(object), collapse=","), "\n")
  cat("  o expected variable labels (as returned by fread):\n")
  cat("        ", cnvar(object), "\n")
  cat("        ", bafvar(object), "\n")
  cat("        ", paste(gtvar(object), collapse=", "), "\n")
  cat("  o numeric data will be multiplied by ", scale(object), " and written\n    to disk as integers\n")
  cat("  o Accessors: genomeIndex(), selectCols(), scale(), cnvar(), bafvar(), gtvar()\n")
})

cnvar <- function(object) object@cnvar
bafvar <- function(object) object@bafvar
gtvar <- function(object) object@gtvar
