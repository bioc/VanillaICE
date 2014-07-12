#' @export
ArrayViews <- function(class="ArrayViews",
                       colData,
                       rowData=GRanges(),
                       filePaths=character(),
                       index=seq_along(rowData),
                       scale=1000,
                       sample_ids,
                       datadir){
  if(missing(colData)){
    if(!missing(sample_ids)) {
      colData <- DataFrame(row.names=sample_ids)
    } else colData <- DataFrame(row.names=basename(filePaths))
  }
  if(missing(datadir) && !missing(filePaths)) datadir <- dirname(filePaths)[1] else datadir <- character()
  new(class,
      colData=colData,
      rowData=rowData,
      index=index,
      filePaths=filePaths,
      scale=scale,
      datadir=datadir)
}

GStudioViews <- function(...)  ArrayViews(class="GStudioViews", ...)

#' @export "["
setMethod("[", signature(x="ArrayViews", i="ANY", j="ANY"), function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@rowData <- rowData(x)[i]
    x@index <- indexGenome(x)[i]
  }
  if(!missing(j)){
    x@colData <- colData(x)[j, ]
    x@filePaths <- x@filePaths[j]
  }
  x
})

setMethod("rowData", "ArrayViews", function(x, ...) x@rowData)
setMethod("colData", "ArrayViews", function(x, ...) x@colData)
setMethod("scale", "ArrayViews", function(x, center=TRUE, scale=TRUE) x@scale)
setMethod("rownames", "ArrayViews", function(x, do.NULL=TRUE, prefix="col") .rownames(x))
setMethod("colnames", "ArrayViews", function(x, do.NULL=TRUE, prefix="col") .colnames(x))
setMethod("indexGenome", "ArrayViews", function(object) object@index)
setMethod("gstudioPaths", "GStudioViews", function(object) object@filePaths)
setMethod("filePaths", "ArrayViews", function(object) object@filePaths)

setMethod("show", "ArrayViews", function(object){
  cat(paste0("class '", class(object), "'\n"))
  cat("   No. files  :", .ncol(object), "\n")
  cat("   No. markers:", .nrow(object), "\n")
})

setMethod("show", "GStudioViews", function(object){
  callNextMethod()
})

setMethod("parseGStudio", c("GStudioViews", "GStudioScanParams"),
          function(object, param){
  outfiles <- lowlevelFiles(object)
  if(all(file.exists(outfiles))) return(NULL)
  file <- gstudioPaths(object)
  nms <- .rownames(object)
  is_gz <- length(grep(".gz", file)) > 0
  if(is_gz){
    ## unzip in a temporary directory using a system call (platform dependent)
    to <- paste0(tempfile(), ".gz")
    file.copy(file, to)
    system(paste("gunzip", to))
    file <- gsub(".gz", "", to)
  }
  ##dat <- fread(to, select=c(1, 10, 11, 18, 19), showProgress=FALSE)
  dat <- fread(file, select=selectCols(param), showProgress=FALSE)
  dat <- dat[indexGenome(param), ]
  nms <- dat[["SNP Name"]]
  stopifnot(identical(nms, .rownames(object)))
  gtindex <- match(gtvar(param), colnames(dat))
  gt <- sapply(gtindex, function(i, x) x[[i]], x=dat)
  gt <- paste0(gt[,1], gt[,2])
  if(!all(gt %in% c("AA", "AB", "BB"))){
    msg <- which(!gt %in% c("AA", "AB", "BB"))
    gt[msg] <- NA
  }
  gt <- as.integer(factor(gt, levels=c("AA", "AB", "BB")))
  r <- scaleBy(dat[["Log R Ratio"]], scale(param))
  b <- scaleBy(dat[["B Allele Freq"]], scale(param))
  saveRDS(r, file=outfiles[1])
  saveRDS(b, file=outfiles[2])
  saveRDS(gt, file=outfiles[3])
  TRUE
})

#' @export
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
  J <- seq_len(.ncol(X))
  j <- NULL
  ##answer <- foreach(j = J, .packages=c("Klein")) %dopar% {
  answer <- foreach(j=J, .packages="VanillaICE") %dopar% {
    FUN(X[, j], ...)
  }
  answer
}

setMethod("lapply", "ArrayViews", function(X, FUN, ...){
  ## Apply FUN to each element in X.  Assumes
  browser()
  FUN <- match.fun(FUN)
  J <- seq_len(.ncol(X))
  j <- NULL
  answer <- foreach(j = J, .packages=c("VanillaICE")) %dopar% {
    FUN(X[, j], ...)
  }
  answer
})

readAssays <- function(object, files){
  keep <- file.exists(files)
  files <- files[keep]
  result <- matrix(NA, .nrow(object), length(files))
  i <- indexGenome(object)
  for(j in seq_along(files)) result[, j] <- readRDS(files[j])[i]
  dimnames(result) <- list(.rownames(object), .colnames(object)[keep])
  result
}


setMethod("lrr", "ArrayViews", function(object){
  files <- lrrFile(object)
  x <- readAssays(object, files)
  scaleRead(x, scale(object))
})

setMethod("baf", "ArrayViews", function(object){
  files <- bafFile(object)
  x <- readAssays(object, files)
  scaleRead(x, scale(object))
})

setMethod("genotypes", "ArrayViews", function(object){
  files <- gtFile(object)
  readAssays(object, files)
})

setMethod("SnpExperiment", "ArrayViews", function(object){
  view <- object[, 1]
  r <- as.matrix(lrr(view))
  b <- as.matrix(baf(view))
  g <- genotypes(view)
  gr <- SnpGRanges(rowData(view), isSnp=rep(TRUE, .nrow(view)))
  SnpArrayExperiment(cn=r, baf=b, rowData=gr, colData=colData(view))
})

#' @export
writeHmm <- function(object){
  file <- hmmFile(.colnames(object))
  if(file.exists(file)) return(TRUE)
  gr <- hmm2(object)[[1]]
  saveRDS(gr, file=file)
  TRUE
}


setMethod("fileName", "ArrayViews", function(object, label){
  file.path(datadir(object), paste0(.colnames(object), "_", label, ".rds"))
})



setMethod("datadir", "ArrayViews", function(object) object@datadir)


setMethod("lrrFile", "ArrayViews", function(object, label="lrr"){
  fileName(object, label)
})

setMethod("bafFile", "ArrayViews", function(object, label="baf"){
  fileName(object, label)
})

setMethod("gtFile", "ArrayViews", function(object, label="gt"){
  fileName(object, label)
})



##lrrFile <- function(object, label="lrr") fileName(object, label)
##bafFile <- function(object, label="baf") fileName(object, label)
##gtFile <- function(object, label="gt") fileName(object, label)
hmmFile <- function(object, label="hmm") fileName(object, label)

lowlevelFiles <- function(views){
  c(lrrFile(views),
    bafFile(views),
    gtFile(views))
}

#' @export
dropSexChrom <- function(object){
  chrom <- as.character(seqnames(rowData(object)))
  is_autosome <- chrom %in% paste0("chr", 1:22)
  if(all(is_autosome))  return(object)
  message("Dropping sex chromosomes...")
  object[is_autosome, ]
}

setMethod("seqnames", "ArrayViews", function(x) seqnames(rowData(x)))
setMethod("start", "ArrayViews", function(x) start(rowData(x)))
setMethod("end", "ArrayViews", function(x) end(rowData(x)))

#' @export
dropDuplicatedMapLocs <- function(object){
  starts <- paste0(as.character(seqnames(object)), start(object), sep="_")
  dups <- duplicated(starts)
  if(!any(dups)) return(object)
  object[!dups, ]
}

setMethod("sort", "ArrayViews", function(x, decreasing=FALSE, ...){
  index <- order(rowData(x))
  if(identical(index, seq_len(.nrow(x)))) return(x)
  message("Sorting views object by genomic position...")
  x[index,]
})


setMethod("scaleBy", c("numeric", "numeric"), function(x, by) as.integer(x*by))
setMethod("scaleRead", c("numeric", "numeric"), function(x, params) x/params)
setMethod("scaleRead", c("matrix", "numeric"), function(x, params) x/params)

.rownames <- function(object) names(rowData(object))
.colnames <- function(object) rownames(colData(object))
.nrow <- function(x) length(rowData(x))
.ncol <- function(x) length(filePaths(x))
.path <- function(object) object@path

setAs("ArrayViews", "SnpArrayExperiment", function(from, to){
  r <- lrr(from)
  b <- baf(from)
  g <- genotypes(from)
  SnpArrayExperiment(cn=r, baf=b, genotypes=g, rowData=SnpGRanges(rowData(from), isSnp=rep(TRUE, nrow(b))),
                     colData=colData(from))

})
