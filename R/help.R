#' A hidden markov model for detection of germline copy number variants from arrays
#'
#' @docType package
#' @name VanillaICE
#' @import oligoClasses
#' @importFrom oligoClasses baf lrr
#' @importFrom crlmm calculateRBaf
#' @importMethodsFrom oligoClasses state baf lrr copyNumber
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import methods
#' @import IRanges
#' @import S4Vectors
#' @import GenomeInfoDb
#' @import BiocGenerics
#' @import MatrixGenerics
#' @importMethodsFrom Biobase pData fData phenoData featureData featureNames sampleNames assayData snpCallProbability
#' @importFrom Biobase rowMax assayDataElement assayDataElementNames assayDataNew
#' @importFrom data.table fread
#' @import foreach
#' @import grid
#' @import lattice
#' @importFrom matrixStats rowMedians colMedians
#' @importFrom tools file_path_sans_ext
#' @importFrom stats  acf dnorm dunif pnorm rnorm runif setNames
#' @importFrom utils packageDescription read.table
#' @useDynLib VanillaICE
NULL
