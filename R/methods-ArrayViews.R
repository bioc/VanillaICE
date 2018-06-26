#' @include help.R AllGenerics.R AllClasses.R
NULL

#' @param class  character string
#' @param colData DataFrame
#' @param rowRanges GRanges object
#' @param sourcePaths character string provide complete path to plain text source files (one file per sample) containing log R ratios and B allele frequencies
#' @param scale  log R ratios and B allele frequencies can be stored as integers on disk to increase IO speed.   If scale =1, the raw data is not transformed.  If scale = 1000 (default), the log R ratios and BAFs are multipled by 1000 and coerced to an integer.
#' @param sample_ids character vector indicating how to name samples.  Ignored if colData is specified.
#' @param parsedPath character vector indicating where parsed files
#' should be saved
#' @param lrrFiles character vector of file names for storing log R ratios
#' @param bafFiles character vector of file names for storing BAFs
#' @param gtFiles character vector of file names for storing genotypes
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
#'   bsgenome <- BSgenome.Hsapiens.UCSC.hg18
#'   seqlevels(fgr) <- seqlevels(bsgenome)[seqlevels(bsgenome) %in% seqlevels(fgr)]
#'   seqinfo(fgr) <- seqinfo(bsgenome)[seqlevels(fgr),]
#'   fgr <- sort(fgr)
#'   files <- list.files(extdir, full.names=TRUE, recursive=TRUE, pattern="FinalReport")
#'   ids <- gsub(".rds", "", gsub("FinalReport", "", basename(files)))
#'   views <- ArrayViews(rowRanges=fgr,
#'                       sourcePaths=files,
#'                       sample_ids=ids)
#'   lrrFile(views)
#'   ## view of first 10 markers and samples 3 and 5
#'   views <- views[1:10, c(3,5)]
#' @export
ArrayViews <- function(class="ArrayViews",
                       colData,
                       rowRanges=GRanges(),
                       sourcePaths=character(),
                       scale=1000,
                       sample_ids,
                       parsedPath=getwd(),
                       lrrFiles=character(),
                       bafFiles=character(),
                       gtFiles=character()){
  if(missing(colData)){
    if(!missing(sample_ids)) {
      colData <- DataFrame(row.names=sample_ids)
    } else colData <- DataFrame(row.names=basename(sourcePaths))
  }
  if(length(sourcePaths) > 0 && length(lrrFiles) == 0){
    stable_file_identifiers <- fileName(sourcePaths, "lrr")
    lrrFiles = file.path(parsedPath, paste0(stable_file_identifiers, "_lrr.rds"))
  }
  if(length(sourcePaths) > 0 && length(bafFiles) == 0){
    stable_file_identifiers <- fileName(sourcePaths, "baf")
    bafFiles = file.path(parsedPath, paste0(stable_file_identifiers, "_baf.rds"))
  }
  if(length(sourcePaths) > 0 && length(gtFiles) == 0){
    stable_file_identifiers <- fileName(sourcePaths, "gt")
    gtFiles = file.path(parsedPath, paste0(stable_file_identifiers, "_gt.rds"))
  }
  new(class,
      colData=colData,
      rowRanges=rowRanges,
      index=seq_len(length(rowRanges)),
      sourcePaths=sourcePaths,
      scale=scale,
      parsedPath=parsedPath,
      lrrFiles=lrrFiles,
      bafFiles=bafFiles,
      gtFiles=gtFiles)
}


setValidity("ArrayViews", function(object){
  msg <- TRUE
  if(length(sourcePaths(object)) > 0){
    if(!all(file.exists(sourcePaths(object)))){
      msg <- "Not all files in sourcePaths(object) exist"
      return(msg)
    }
  }
  ## should we check that files have .rds extension?
  if(length(bafFile(object)) != length(lrrFile(object))){
    msg <- "lrrFiles vector must be the same length as sourcePaths"
    return(msg)
  }
  if(length(bafFile(object)) > 0){
    if(length(sourcePaths(object)) != length(bafFile(object))){
      msg <- "bafFiles vector must be the same length as sourcePaths"
      return(msg)
    }
    if(length(sourcePaths(object)) != length(gtFile(object))){
      msg <- "gtFiles vector must be the same length as sourcePaths"
      return(msg)
    }
  }
  if(length(parsedPath(object)) > 0){
    ddir <- parsedPath(object)
    if(!file.exists(ddir)){
      msg <- "Directory parsedPath(object) does not exist"
      return(msg)
    }
  }
  if(length(object@index) != length(rowRanges(object))){
    msg <- "index slot should have same length as rowRanges"
    return(msg)
  }
  return(msg)
})

setMethod("seqinfo", "ArrayViews", function(x){
  seqinfo(rowRanges(x))
})

setReplaceMethod("seqlevels", "ArrayViews",
  function(x,
           pruning.mode=c("error", "coarse", "fine", "tidy"),
           value){
  i <- setNames(x@index, names(rowRanges(x)))
  rd <- rowRanges(x)
  seqlevels(rd, pruning.mode=pruning.mode) <- value
  rowRanges(x) <- rd
  x@index <- i[names(rd)]
  x
})

#' @aliases ArrayViews,numeric,numeric-method "[",ArrayViews,ANY-method
#' @param i numeric vector or missing
#' @param j numeric vector or missing
#' @param drop ignored
#' @rdname ArrayViews-class
setMethod("[", signature(x="ArrayViews", i="ANY", j="ANY"), function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@rowRanges <- rowRanges(x)[i]
    x@index <- indexGenome(x)[i]
  }
  if(!missing(j)){
    if(is.character(j)) j <- match(j, colnames(x))
    x@colData <- colData(x)[j, ]
    x@sourcePaths <- x@sourcePaths[j]
    ##
    ## We do not want to check whether this slot is a character string
    ##  -- should enforce character class and length of in the validity method
    x@lrrFiles <- x@lrrFiles[j]
    x@bafFiles <- x@bafFiles[j]
    x@gtFiles <- x@gtFiles[j]
  }
  x
})


#' @param value a character-string vector
#' @rdname ArrayViews-class
#' @aliases colnames<-,ArrayViews,character-method
#' @name colnames<-
#' @usage colnames(x) <- value
#' @rdname ArrayViews-class
#' @export
setReplaceMethod("colnames", c("ArrayViews", "character"), function(x, value){
  coldata <- colData(x)
  rownames(coldata) <- value
  x@colData <- coldata
  x
})

#' @aliases colnames colnames,ArrayViews-method
#' @rdname ArrayViews-class
#' @param do.NULL ignored
#' @param prefix ignored
#' @export
setMethod("colnames", "ArrayViews", function(x, do.NULL=TRUE, prefix="col") .colnames(x))

setMethod("rowRanges", "ArrayViews", function(x, ...) x@rowRanges)

setMethod("rowRanges<-", "ArrayViews", function(x, value) {
  x@rowRanges <- value
  x
})

setMethod("colData", "ArrayViews", function(x, ...) x@colData)

setReplaceMethod("colData", "ArrayViews", function(x, value){
  x@colData <- value
  x
})


#' @rdname ArrayViews-class
#' @aliases $,ArrayViews-method
#' @param name character string indicating name in colData slot of
#' ArrayViews object
#' @export
setMethod("$", "ArrayViews", function(x, name) {
  eval(substitute(colData(x)$NAME_ARG, list(NAME_ARG=name)))
})

#' @rdname ArrayViews-class
#' @export
setReplaceMethod("$", "ArrayViews", function(x, name, value) {
  colData(x)[[name]] <- value
  x
})

##setReplaceMethod("$", "ArrayViews", function(x, value){
##  x@colData <- value
##  x
##})

setMethod("scale", "ArrayViews", function(x, center=TRUE, scale=TRUE) x@scale)

setMethod("rownames", "ArrayViews", function(x, do.NULL=TRUE, prefix="col") .rownames(x))



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
  stop("not all files have markers in the same order, or some files are from a different platform")
}

.parseSourceFile <- function(object, param){
  if(ncol(object) > 1) warning("Only parsing the first file in the views object")
  object <- object[,1]
  outfiles <- lowlevelFiles(object)
  if(all(file.exists(outfiles))) return(NULL)
  file <- sourcePaths(object)
  nms <- .rownames(object)
  is_gz <- length(grep(".gz$", file)) > 0
  if(is_gz){
    ## unzip in a temporary directory using a system call (platform dependent)
    to <- paste0(tempfile(), ".gz")
    file.copy(file, to)
    system(paste("gunzip", to))
    file <- gsub(".gz", "", to)
  }
  dat <- fread(file[1], select=selectCols(param), showProgress=FALSE, skip="[Data]")
  dat <- dat[indexGenome(param), ]
  ##nms <- dat[["SNP Name"]]
  nms <- dat[[.snp_id_column(param)]]
  if(!identical(nms, rownames(object))){
    rownames(dat) <- nms
    dat <- .resolveIndex(dat, object)
  }
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
#'   bsgenome <- BSgenome.Hsapiens.UCSC.hg18
#'   require(data.table)
#'   extdir <- system.file("extdata", package="VanillaICE", mustWork=TRUE)
#'   features <- suppressWarnings(fread(file.path(extdir, "SNP_info.csv")))
#'   fgr <- GRanges(paste0("chr", features$Chr), IRanges(features$Position, width=1),
#'                  isSnp=features[["Intensity Only"]]==0)
#'   fgr <- SnpGRanges(fgr)
#'   names(fgr) <- features[["Name"]]
#'   seqlevels(fgr) <- seqlevels(bsgenome)[seqlevels(bsgenome) %in% seqlevels(fgr)]
#'   seqinfo(fgr) <- seqinfo(bsgenome)[seqlevels(fgr),]
#'   fgr <- sort(fgr)
#'   files <- list.files(extdir, full.names=TRUE, recursive=TRUE, pattern="FinalReport")
#'   views <- ArrayViews(rowRanges=fgr, sourcePaths=files, parsedPath=tempdir())
#'   show(views)
#'
#' ## read the first file
#' dat <- fread(files[1], skip="[Data]")
#' ## information to store on the markers
#' select <- match(c("SNP Name", "Allele1 - AB", "Allele2 - AB",
#'                   "Log R Ratio", "B Allele Freq"), names(dat))
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
setMethod("sapply", "ArrayViews", function(X, FUN, ..., simplify=TRUE,
                                           USE.NAMES=TRUE){
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
  x <- scaleRead(x, scale(object))
  x
})

#' @aliases baf,ArrayViews-method
#' @rdname LowLevelSummaries
setMethod("baf", "ArrayViews", function(object){
  files <- bafFile(object)
  x <- readAssays(object, files)
  x <- scaleRead(x, scale(object))
  x
})

#' @aliases genotypes,ArrayViews-method
#' @rdname LowLevelSummaries
setMethod("genotypes", "ArrayViews", function(object){
  files <- gtFile(object)
  gt <- readAssays(object, files)
  gt
})

#' @param x a \code{ArrayViews} object
#' @export
#' @rdname ArrayViews-class
setMethod("ncol", "ArrayViews", function(x) nrow(colData(x)))

#' @export
#' @rdname ArrayViews-class
setMethod("nrow", "ArrayViews", function(x) length(rowRanges(x)))

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
  if(is(rowRanges(view), "SnpGRanges")){
    gr <- rowRanges(view)
  } else {
    gr <- SnpGRanges(rowRanges(view), isSnp=rep(TRUE, nrow(view)))
  }
  SnpArrayExperiment(cn=r, baf=b, rowRanges=gr, colData=colData(view))
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

setMethod("fileName", "character", function(object, label){
  source_paths <- file_path_sans_ext(file_path_sans_ext(basename(object)))
  stable_file_identifiers <- make.unique(source_paths)
})

setMethod("fileName", "ArrayViews", function(object, label){
  ## strip ending
  source_paths <- file_path_sans_ext(file_path_sans_ext(basename(sourcePaths(object))))
  stable_file_identifiers <- make.unique(source_paths)
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
setMethod("lrrFile", "ArrayViews", function(object) object@lrrFiles)

#' @param value a character vector of filenames for the log R ratios
#' @aliases lrrFile<-,ArrayViews-method
#' @rdname IO
setReplaceMethod("lrrFile", "ArrayViews", function(object, value){
  object@lrrFiles <- value
  object
})


#' @aliases bafFile,ArrayViews-method
#' @rdname IO
setMethod("bafFile", "ArrayViews", function(object) object@bafFiles)
##  fileName(object, label)
##})

#' @aliases gtFile,ArrayViews-method
#' @rdname IO
setMethod("gtFile", "ArrayViews", function(object) object@gtFiles)
##  fileName(object, label)
##})

hmmFile <- function(object, label="hmm") fileName(object, label)

## This creates filenames for storing log R ratios, etc.
lowlevelFiles <- function(views){
  files <- c(lrrFile(views), bafFile(views), gtFile(views))
  if(any(is.na(files))) stop("low level file name is invalid")
  files
}

#' Filter sex chromosomes
#'
#' Removes markers on chromosomes X and Y.
#'
#' @param object an object for which the methods \code{seqnames} and \code{rowRanges} are defined.
#' @return an object of the same class as the input
#' @export
dropSexChrom <- function(object){
  chrom <- as.character(seqnames(rowRanges(object)))
  is_autosome <- chrom %in% paste0("chr", 1:22)
  if(all(is_autosome))  return(object)
  message("Dropping sex chromosomes...")
  object[is_autosome, ]
}

setMethod("seqnames", "ArrayViews", function(x) seqnames(rowRanges(x)))

#' @aliases start,ArrayViews-method
#' @rdname ArrayViews-class
setMethod("start", "ArrayViews", function(x) start(rowRanges(x)))

#' @aliases end,ArrayViews-method
setMethod("end", "ArrayViews", function(x) end(rowRanges(x)))


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
#' g <- rowRanges(snp_exp)
#' ## duplicate the first row
#' g[length(g)] <- g[1]
#'  rowRanges(snp_exp) <- g
#'  snp_exp2 <- dropDuplicatedMapLocs(snp_exp)
#' @export
dropDuplicatedMapLocs <- function(object){
  starts <- paste0(as.character(seqnames(object)), start(object), sep="_")
  dups <- duplicated(starts)
  if(!any(dups)) return(object)
  object[!dups, ]
}

setMethod("sort", "ArrayViews", function(x, decreasing=FALSE, ...){
  index <- order(rowRanges(x))
  if(identical(index, seq_len(nrow(x)))) return(x)
  message("Sorting views object by genomic position...")
  x[index,]
})


setMethod("scaleBy", c("numeric", "numeric"), function(x, by) as.integer(round(x*by, 0)))
setMethod("scaleRead", c("numeric", "numeric"), function(x, params) x/params)
setMethod("scaleRead", c("matrix", "numeric"), function(x, params) x/params)

.rownames <- function(object) names(rowRanges(object))
.colnames <- function(object) rownames(colData(object))
##.nrow <- function(x) length(rowRanges(x))
## ncol <- function(x) length(sourcePaths(x))
.path <- function(object) object@path

setAs("ArrayViews", "SnpArrayExperiment", function(from, to){
  r <- lrr(from)
  b <- baf(from)
  g <- genotypes(from)
  SnpArrayExperiment(cn=r, baf=b, genotypes=g, rowRanges=SnpGRanges(rowRanges(from), isSnp=rep(TRUE, nrow(b))),
                     colData=colData(from))

})

setMethod("isAutosome", "ArrayViews", function(object){
  isAutosome(rowRanges(object))
})

setMethod("chromosome", "ArrayViews", function(object) as.character(seqnames(object)))

#' @aliases isHeterozygous,ArrayViews-method
#' @rdname isHeterozygous
setMethod("isHeterozygous", "ArrayViews", function(object, cutoff){
  isHeterozygous(baf(object), cutoff)
})
