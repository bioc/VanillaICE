#' @param cnvar length-one character vector providing name of variable for log R ratios
#' @param bafvar length-one character vector providing name of variable for B allele frequencies
#' @param gtvar length-one character vector providing name of variable for genotype calls
#' @param index_genome integer vector indicating which rows of the of
#' the source files (e.g., GenomeStudio) to keep.  By matching on a sorted GRanges
#' object containing the feature annotation (see example), the
#' information on the markers will also be sorted.
#' @param select integer vector specifying indicating which columns of the source files to import (see examples)
#' @param scale length-one numeric vector for rescaling the raw data
#' and coercing to class integer. By default, the low-level data will
#' be scaled and saved on disk as integers.
#' @param row.names length-one numeric vector indicating which column
#' the SNP names are in
#' @examples
#' CopyNumScanParams() ## empty container
#' @rdname CopyNumScanParams
#' @seealso \code{\linkS4class{ArrayViews}} \code{\link{parseSourceFile}}
#' @export
CopyNumScanParams <- function(cnvar="Log R Ratio",
                              bafvar="B Allele Freq",
                              gtvar=c("Allele1 - AB", "Allele2 - AB"),
                              index_genome=integer(),
                              select=integer(),
                              scale=1000,
                              row.names=1L){
  new("CopyNumScanParams", cnvar=cnvar, bafvar=bafvar, gtvar=gtvar,
      index_genome=index_genome, scale=scale, select=select, row.names=row.names)
}

setValidity("CopyNumScanParams", function(object){
  msg <- TRUE
  if(any(is.na(indexGenome(object)))) {
    msg <- "indexGenome values can not contain NAs"
    return(msg)
  }
  msg
})

setMethod("scale", "CopyNumScanParams", function(x, center=TRUE, scale=TRUE) x@scale)
setMethod("indexGenome", "CopyNumScanParams", function(object) object@index_genome)
setMethod("selectCols", "CopyNumScanParams", function(object) object@select)

#' @param object a \code{CopyNumScanParams} object
#' @aliases show,CopyNumScanParams-method
#' @rdname CopyNumScanParams
setMethod("show", "CopyNumScanParams", function(object){
  cat("'CopyNumScanParams' class:\n")
  cat("  o columns to select:", paste0(selectCols(object), collapse=","), "\n")
  cat("  o expected variable labels (as returned by fread):\n")
  cat("        ", cnvar(object), "\n")
  cat("        ", bafvar(object), "\n")
  cat("        ", paste(gtvar(object), collapse=", "), "\n")
  cat("  o numeric data will be multiplied by ", scale(object), " and written\n    to disk as integers\n")
  cat("  o Accessors: indexGenome(), selectCols(), scale(), cnvar(), bafvar(), gtvar()\n")
})

cnvar <- function(object) object@cnvar
bafvar <- function(object) object@bafvar
gtvar <- function(object) object@gtvar
