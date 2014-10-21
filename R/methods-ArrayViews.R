#' @include help.R AllGenerics.R AllClasses.R
NULL

#' @param class  character string
#' @param colData DataFrame
#' @param rowData GRanges object
#' @param sourcePaths character string provide complete path to plain text source files (one file per sample) containing log R ratios and B allele frequencies
#' @param scale  log R ratios and B allele frequencies can be stored as integers on disk to increase IO speed.   If scale =1, the raw data is not transformed.  If scale = 1000 (default), the log R ratios and BAFs are multipled by 1000 and coerced to an integer.
#' @param sample_ids character vector indicating how to name samples
#' @param parsedPath character vector indicating where parsed files should be saved
#' @seealso \code{\link{CopyNumScanParams}} \code{\link{parseSourceFile}}
#' @aliases ArrayViews
#' @rdname ArrayViews-class
#' @examples
#' ArrayViews()
#' ## From unit test
#'   require(BSgenome.Hsapiens.UCSC.hg18)
#'   require(data.table)
#'   extdir <- system.file("extdata", package="VanillaICE", mustWork=TRUE)
#'   features <- suppressWarnings(fread(file.path(extdir, "SNP_info.csv")))
#'   fgr <- GRanges(paste0("chr", features$Chr), IRanges(features$Position, width=1),
#'                  isSnp=features[["Intensity Only"]]==0)
#'   fgr <- SnpGRanges(fgr)
#'   names(fgr) <- features[["Name"]]
#'   seqlevels(fgr) <- seqlevels(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(BSgenome.Hsapiens.UCSC.hg18) %in% seqlevels(fgr)]
#'   seqinfo(fgr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(fgr),]
#'   fgr <- sort(fgr)
#'   files <- list.files(extdir, full.names=TRUE, recursive=TRUE, pattern="FinalReport")
#'   ids <- gsub(".rds", "", gsub("FinalReport", "", basename(files)))
#'   views <- ArrayViews(rowData=fgr,
#'                       sourcePaths=files,
#'                       sample_ids=ids)
#' @export
ArrayViews <- function(class="ArrayViews",
                       colData,
                       rowData=GRanges(),
                       sourcePaths=character(),
                       scale=1000,
                       sample_ids,
                       parsedPath="./"){
  if(missing(colData)){
    if(!missing(sample_ids)) {
      colData <- DataFrame(row.names=sample_ids)
    } else colData <- DataFrame(row.names=basename(sourcePaths))
  }
  new(class,
      colData=colData,
      rowData=rowData,
      index=seq_len(length(rowData)),
      sourcePaths=sourcePaths,
      scale=scale,
      parsedPath=parsedPath)
}

##ArrayViews <- function(class="CNSet",
##                       colData,
##                       rowData=GRanges(),
##                       sourcePaths=character(),
##                       ##index=seq_along(rowData),
##                       scale=1000,
##                       sample_ids,
##                       parsedPath){
##  if(missing(colData)){
##    if(!missing(sample_ids)) {
##      colData <- DataFrame(row.names=sample_ids)
##    } else colData <- DataFrame(row.names=basename(sourcePaths))
##  }
##  if(missing(parsedPath)){
##    parsedPath <- if(!missing(sourcePaths)) dirname(sourcePaths)[1] else character()
##  }
##  new(class,
##      colData=colData,
##      rowData=rowData,
##      index=seq_len(length(rowData)),
##      sourcePaths=sourcePaths,
##      scale=scale,
##      parsedPath=parsedPath)
##}

setValidity("ArrayViews", function(object){
  msg <- TRUE
  if(length(sourcePaths(object)) > 0){
    if(!all(file.exists(sourcePaths(object)))){
      msg <- "Not all files in sourcePaths(object) exist"
      return(msg)
    }
  }
  if(length(parsedPath(object)) > 0){
    ddir <- parsedPath(object)
    if(!file.exists(ddir)){
      msg <- "Directory parsedPath(object) does not exist"
    }
  }
  return(msg)
})

#' @aliases ArrayViews,numeric,numeric-method "[",ArrayViews,ANY-method
#' @param i numeric vector or missing
#' @param j numeric vector or missing
#' @param drop ignored
#' @rdname ArrayViews-class
setMethod("[", signature(x="ArrayViews", i="ANY", j="ANY"), function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@rowData <- rowData(x)[i]
    x@index <- indexGenome(x)[i]
  }
  if(!missing(j)){
    if(is.character(j)) j <- match(j, colnames(x))
    x@colData <- colData(x)[j, ]
    x@sourcePaths <- x@sourcePaths[j]
  }
  x
})


#' @param value a character-string vector
#' @rdname ArrayViews-class
#' @aliases colnames<-,ArrayViews,character-method
#' @name colnames<-
#' @usage colnames(x) <- value
#' @export
setReplaceMethod("colnames", c("ArrayViews", "character"), function(x, value){
  coldata <- colData(x)
  rownames(coldata) <- value
  x@colData <- coldata
  x
})

setMethod("rowData", "ArrayViews", function(x, ...) x@rowData)
setMethod("colData", "ArrayViews", function(x, ...) x@colData)
setMethod("scale", "ArrayViews", function(x, center=TRUE, scale=TRUE) x@scale)
setMethod("rownames", "ArrayViews", function(x, do.NULL=TRUE, prefix="col") .rownames(x))
setMethod("colnames", "ArrayViews", function(x, do.NULL=TRUE, prefix="col") .colnames(x))
setMethod("indexGenome", "ArrayViews", function(object) object@index)
##setMethod("gstudioPaths", "GStudioViews", function(object) object@sourcePaths)
setMethod("sourcePaths", "ArrayViews", function(object) object@sourcePaths)

#' @param object a \code{ArrayViews} object
#' @aliases show,ArrayViews-method
#' @rdname ArrayViews-class
setMethod("show", "ArrayViews", function(object){
  cat(paste0("class '", class(object), "'\n"))
  cat("   No. files  :", ncol(object), "\n")
  cat("   No. markers:", nrow(object), "\n")
})

.snp_id_column <- function(object) object@row.names

.resolveIndex <- function(object, param){
  browser()
}

.parseSourceFile <- function(object, param){
  if(ncol(object) > 1) warning("Only parsing the first file in the views object")
  object <- object[,1]
  outfiles <- lowlevelFiles(object)
  if(all(file.exists(outfiles))) return(NULL)
  file <- sourcePaths(object)
  nms <- .rownames(object)
  is_gz <- length(grep(".gz", file)) > 0
  if(is_gz){
    ## unzip in a temporary directory using a system call (platform dependent)
    to <- paste0(tempfile(), ".gz")
    file.copy(file, to)
    system(paste("gunzip", to))
    file <- gsub(".gz", "", to)
  }
  dat <- fread(file[1], select=selectCols(param), showProgress=FALSE)
  dat <- dat[indexGenome(param), ]
  ##nms <- dat[["SNP Name"]]
  nms <- dat[[.snp_id_column(param)]]
##  if(!identical(nms, .rownames(object))){
##    dat <- .resolveIndex(dat, object)
##  }
  stopifnot(identical(nms, .rownames(object)))
  gtindex <- match(gtvar(param), colnames(dat))
  if(length(gtvar(param))==2){
    gt <- sapply(gtindex, function(i, x) x[[i]], x=dat)
    gt <- paste0(gt[,1], gt[,2])
    if(!all(gt %in% c("AA", "AB", "BB"))){
      msg <- which(!gt %in% c("AA", "AB", "BB"))
      gt[msg] <- NA
    }
  } else gt <- dat[[gtindex]]
  if(is.character(gt)){
    gt <- as.integer(factor(gt, levels=c("AA", "AB", "BB")))
  } else gt <- as.integer(gt)
  j <- match(cnvar(param), colnames(dat))
  k <- match(bafvar(param), colnames(dat))
  r <- scaleBy(dat[[j]], scale(param))
  b <- scaleBy(dat[[k]], scale(param))
  saveRDS(r, file=outfiles[1])
  saveRDS(b, file=outfiles[2])
  saveRDS(gt, file=outfiles[3])
  NULL
}

#' @aliases parseSourceFile,ArrayViews,CopyNumScanParams-method
#' @rdname parseSourceFile
#' @examples
#'   require(BSgenome.Hsapiens.UCSC.hg18)
#'   require(data.table)
#'   extdir <- system.file("extdata", package="VanillaICE", mustWork=TRUE)
#'   features <- suppressWarnings(fread(file.path(extdir, "SNP_info.csv")))
#'   fgr <- GRanges(paste0("chr", features$Chr), IRanges(features$Position, width=1),
#'                  isSnp=features[["Intensity Only"]]==0)
#'   fgr <- SnpGRanges(fgr)
#'   names(fgr) <- features[["Name"]]
#'   seqlevels(fgr) <- seqlevels(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(BSgenome.Hsapiens.UCSC.hg18) %in% seqlevels(fgr)]
#'   seqinfo(fgr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(fgr),]
#'   fgr <- sort(fgr)
#'   files <- list.files(extdir, full.names=TRUE, recursive=TRUE, pattern="FinalReport")
#'   views <- ArrayViews(rowData=fgr, sourcePaths=files, parsedPath=tempdir())
#'   show(views)
#'
#' ## read the first file
#' dat <- fread(files[1])
#' ## information to store on the markers
#' select <- match(c("SNP Name", "Allele1 - AB", "Allele2 - AB", "Log R Ratio", "B Allele Freq"), names(dat))
#' ##
#' ## which rows to keep in the MAP file. By matching on the sorted GRanges object
#' ## containing the feature annotation, the low-level data for the log R ratios/
#' ## B allele frequencies will also be sorted
#' ##
#' index_genome <- match(names(fgr), dat[["SNP Name"]])
#' scan_params <- CopyNumScanParams(index_genome=index_genome, select=select)
#' ##
#' ## parse the source files
#' ##
#' parseSourceFile(views, scan_params)
#' list.files(parsedPath(views))
#' ##
#' ##  Inspecting source data through accessors defined on the views object
#' ##
#' require(oligoClasses)
#' ## log R ratios
#' r <- head(lrr(views))
#' ## B allele frequencies
#' b <- head(baf(views))
#' g <- head(genotypes(views))
setMethod("parseSourceFile", c("ArrayViews", "CopyNumScanParams"),
          function(object, param) {
            message("Writing parsed files to ", parsedPath(object))
            invisible(sapply(object, .parseSourceFile, param))
          })

#' @export
#' @aliases sapply,ArrayViews-method
#' @param X a \code{ArrayViews} object
#' @param FUN a function to apply to each column of \code{X}
#' @param simplify logical indicating whether result should be simplied
#' @param USE.NAMES whether the output should be a named vector
#' @param ... additional arguments to \code{FUN}
#' @rdname ArrayViews-class
setMethod("sapply", "ArrayViews", function(X, FUN, ...,
                                             simplify=TRUE, USE.NAMES=TRUE){
  FUN <- match.fun(FUN)
  answer <- .lapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
    names(answer) <- X
  if (!identical(simplify, FALSE) && length(answer))
    simplify2array(answer, higher = (simplify == "array"))
  else answer
})

.lapply <- function(X, FUN, ..., simplify=FALSE, USE.NAMES=FALSE){
  FUN <- match.fun(FUN)
  J <- seq_len(ncol(X))
  j <- NULL
  ##answer <- foreach(j = J, .packages=c("Klein")) %dopar% {
  answer <- foreach(j=J, .packages="VanillaICE") %dopar% {
    FUN(X[, j], ...)
  }
  answer
}

setMethod("lapply", "ArrayViews", function(X, FUN, ...){
  ## Apply FUN to each element in X.  Assumes
  FUN <- match.fun(FUN)
  J <- seq_len(ncol(X))
  j <- NULL
  answer <- foreach(j = J, .packages=c("VanillaICE")) %dopar% {
    FUN(X[, j], ...)
  }
  answer
})

readAssays <- function(object, files){
  keep <- file.exists(files)
  files <- files[keep]
  result <- matrix(NA, nrow(object), length(files))
  i <- indexGenome(object)
  for(j in seq_along(files)) result[, j] <- readRDS(files[j])[i]
  dimnames(result) <- list(.rownames(object), .colnames(object)[keep])
  result
}


#' @aliases lrr,ArrayViews-method
#' @rdname LowLevelSummaries
setMethod("lrr", "ArrayViews", function(object){
  files <- lrrFile(object)
  x <- readAssays(object, files)
  scaleRead(x, scale(object))
})

#' @aliases baf,ArrayViews-method
#' @rdname LowLevelSummaries
setMethod("baf", "ArrayViews", function(object){
  files <- bafFile(object)
  x <- readAssays(object, files)
  scaleRead(x, scale(object))
})

#' @aliases genotypes,ArrayViews-method
#' @rdname LowLevelSummaries
setMethod("genotypes", "ArrayViews", function(object){
  files <- gtFile(object)
  readAssays(object, files)
})

#' @param x a \code{ArrayViews} object
#' @export
#' @rdname ArrayViews-class
setMethod("ncol", "ArrayViews", function(x) nrow(colData(x)))

#' @export
#' @rdname ArrayViews-class
setMethod("nrow", "ArrayViews", function(x) length(rowData(x)))

#' @export
#' @rdname ArrayViews-class
setMethod("dim", "ArrayViews", function(x) c(ncol(x), nrow(x)))

#' @rdname SnpExperiment
#' @aliases SnpExperiment,ArrayViews-method
#' @examples
#' view <- ArrayViews()
#' SnpExperiment(view)
setMethod("SnpExperiment", "ArrayViews", function(object){
  if(ncol(object) == 0) return(SnpArrayExperiment())
  view <- object
  r <- as.matrix(lrr(view))
  b <- as.matrix(baf(view))
  g <- as.matrix(genotypes(view))
  gr <- SnpGRanges(rowData(view), isSnp=rep(TRUE, nrow(view)))
  SnpArrayExperiment(cn=r, baf=b, rowData=gr, colData=colData(view))
})

writeHmm <- function(object){
  file <- hmmFile(.colnames(object))
  if(file.exists(file)) return(TRUE)
  gr <- hmm2(object)[[1]]
  saveRDS(gr, file=file)
  TRUE
}

#' @param tolerance length-one numeric vector.  When the difference in
#' the log-likelihood of the Viterbi state path between successive
#' models (updated by Baum Welch) is less than the tolerance, no
#' additional model updates are performed.
#' @param verbose logical.  Whether to display messages indicating progress.
#' @aliases hmm2,ArrayViews-method
#' @rdname hmm2
setMethod("hmm2", "ArrayViews", function(object, emission_param=EmissionParam(),
                                         transition_param=TransitionParam(),
                                         tolerance=2,
                                         verbose=FALSE, ...){
  se <- as(object, "SnpArrayExperiment")
  hmm2(se, emission_param=emission_param,
       transition_param=transition_param,
       tolerance=tolerance,
       verbose=verbose, ...)
})


setMethod("fileName", "ArrayViews", function(object, label){
  stable_file_identifiers <- make.unique(basename(sourcePaths(object)))
  file.path(parsedPath(object), paste0(stable_file_identifiers, "_", label, ".rds"))
})


#' @aliases parsedPath,ArrayViews-method
#' @rdname parsedPath
setMethod("parsedPath", "ArrayViews", function(object) object@parsedPath)


#' @aliases lrrFile,ArrayViews-method
#' @rdname IO
#' @examples
#' views <- ArrayViews(parsedPath=tempdir())
#' sourcePaths(views)
#' lrrFile(views)
#' bafFile(views)
#' gtFile(views)
setMethod("lrrFile", "ArrayViews", function(object, label="lrr"){
  fileName(object, label)
})

#' @aliases bafFile,ArrayViews-method
#' @rdname IO
setMethod("bafFile", "ArrayViews", function(object, label="baf"){
  fileName(object, label)
})

#' @aliases gtFile,ArrayViews-method
#' @rdname IO
setMethod("gtFile", "ArrayViews", function(object, label="gt"){
  fileName(object, label)
})



##lrrFile <- function(object, label="lrr") fileName(object, label)
##bafFile <- function(object, label="baf") fileName(object, label)
##gtFile <- function(object, label="gt") fileName(object, label)
hmmFile <- function(object, label="hmm") fileName(object, label)

## This creates filenames for storing log R ratios, etc.
lowlevelFiles <- function(views){
  c(lrrFile(views), bafFile(views), gtFile(views))
}

#' Filter sex chromosomes
#'
#' Removes markers on chromosomes X and Y.
#'
#' @param object an object for which the methods \code{seqnames} and \code{rowData} are defined.
#' @return an object of the same class as the input
#' @export
dropSexChrom <- function(object){
  chrom <- as.character(seqnames(rowData(object)))
  is_autosome <- chrom %in% paste0("chr", 1:22)
  if(all(is_autosome))  return(object)
  message("Dropping sex chromosomes...")
  object[is_autosome, ]
}

setMethod("seqnames", "ArrayViews", function(x) seqnames(rowData(x)))

#' @aliases start,ArrayViews-method
#' @rdname ArrayViews-class
setMethod("start", "ArrayViews", function(x) start(rowData(x)))

#' @aliases end,ArrayViews-method
setMethod("end", "ArrayViews", function(x) end(rowData(x)))


#' Drop markers on the same chromosome having the same genomic
#' coordinates
#'
#' If there are multiple markers on the same chromosome with the same
#' annotated position, only the first is kept.
#'
#' @param object a container for which the methods seqnames and start
#' are defined
#' @return an object of the same class with duplicated genomic positions removed
#' @examples
#' data(snp_exp)
#' g <- rowData(snp_exp)
#' ## duplicate the first row
#' g[length(g)] <- g[1]
#'  rowData(snp_exp) <- g
#'  snp_exp2 <- dropDuplicatedMapLocs(snp_exp)
#' @export
dropDuplicatedMapLocs <- function(object){
  starts <- paste0(as.character(seqnames(object)), start(object), sep="_")
  dups <- duplicated(starts)
  if(!any(dups)) return(object)
  object[!dups, ]
}

setMethod("sort", "ArrayViews", function(x, decreasing=FALSE, ...){
  index <- order(rowData(x))
  if(identical(index, seq_len(nrow(x)))) return(x)
  message("Sorting views object by genomic position...")
  x[index,]
})


setMethod("scaleBy", c("numeric", "numeric"), function(x, by) as.integer(x*by))
setMethod("scaleRead", c("numeric", "numeric"), function(x, params) x/params)
setMethod("scaleRead", c("matrix", "numeric"), function(x, params) x/params)

.rownames <- function(object) names(rowData(object))
.colnames <- function(object) rownames(colData(object))
##.nrow <- function(x) length(rowData(x))
## ncol <- function(x) length(sourcePaths(x))
.path <- function(object) object@path

setAs("ArrayViews", "SnpArrayExperiment", function(from, to){
  r <- lrr(from)
  b <- baf(from)
  g <- genotypes(from)
  SnpArrayExperiment(cn=r, baf=b, genotypes=g, rowData=SnpGRanges(rowData(from), isSnp=rep(TRUE, nrow(b))),
                     colData=colData(from))

})
