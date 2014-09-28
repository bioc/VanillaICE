#' @examples
#' SnpArrayExperiment()
#' @aliases   SnpArrayExperiment,missing-method
#' @rdname SnpArrayExperiment-class
setMethod(SnpArrayExperiment, "missing",
          function(cn, baf, rowData=GRanges(),
                   colData=DataFrame(), isSnp=logical(), ...){
            se <- new("SnpArrayExperiment", assays=snpArrayAssays(...),
                      rowData=SnpGRanges(), ...)
            se
          })


#' @aliases SnpArrayExperiment,matrix-method
#' @rdname SnpArrayExperiment-class
setMethod(SnpArrayExperiment, "matrix",
          function(cn, baf, rowData=GRanges(),
                   colData=DataFrame(row.names=colnames(cn)),
                   isSnp=logical(), ...){
            assays <- snpArrayAssays(cn=cn, baf=baf, ...)
            if("isSnp" %in% colnames(mcols(rowData))){
              rowData <- SnpGRanges(rowData)
            } else {
              if(length(isSnp) != length(rowData)) stop(" isSnp must be the same length as rowData")
              rowData <- SnpGRanges(rowData, isSnp)
            }
            se <- new("SnpArrayExperiment", assays=assays,
                      rowData=rowData, colData=colData)
            se
          })


setAs("CNSet", "SnpArrayExperiment",
      function(from){
        b.r <- crlmm:::calculateRBaf(from)
        baflist <- b.r[["baf"]]
        not_matrix <- !is.matrix(baflist[[1]])
        if(not_matrix) baflist <- lapply(baflist, "[")
        lrrlist <- b.r[["lrr"]]
        if(not_matrix) lrrlist <- lapply(lrrlist, "[")
        bafs <- do.call(rbind, baflist)
        lrrs <- do.call(rbind, lrrlist)
        nms <- rownames(lrrs)
        i <- match(nms, featureNames(from))
        g <- GRanges(paste0("chr", integer2chromosome(chromosome(from)[i])),
                     IRanges(position(from)[i], width=1))
        rowdat <- SnpGRanges(g, isSnp=isSnp(from)[i])
        rowdat <- SnpGRanges(rowdat)
        names(rowdat) <- nms
        ##tmp=snpArrayAssays(cn=lrrs/100, baf=bafs/100)
        ##new("SnpArrayExperiment", assays=tmp, rowData=rowdat, colData=DataFrame(pData(from)))
        se <- SnpArrayExperiment(cn=lrrs/100, baf=bafs/1000,
                                 rowData=SnpGRanges(rowdat),
                                 colData=DataFrame(pData(from)))
        sort(se)
      })

setAs("oligoSnpSet", "SnpArrayExperiment",
      function(from){
        rowdata <- as(featureData(from), "SnpGRanges")
        coldata <- as(as(phenoData(from), "data.frame"), "DataFrame")
        cn_assays <- snpArrayAssays(cn=copyNumber(from),
                                    baf=baf(from),
                                    gt=calls(from))
        new("SnpArrayExperiment", assays=cn_assays, rowData=rowdata, colData=coldata)
      })



setValidity("SnpArrayExperiment", function(object){
  msg <- validObject(rowData(object))
  msg <- msg && (requiredAssays() %in% names(assays(object)))
  msg
})

#' @export
#' @aliases baf,SnpArrayExperiment-method
#' @rdname LowLevelSummaries
setMethod("baf", "SnpArrayExperiment",
          function(object) assays(object)[["baf"]])

setMethod("chromosome", "SnpArrayExperiment",
          function(object) as.character(seqnames(object)))

#' @export
#' @aliases copyNumber,SnpArrayExperiment-method
#' @rdname LowLevelSummaries
setMethod("copyNumber", "SnpArrayExperiment",
          function(object) assays(object)[["cn"]])

setMethod("isSnp", "SnpArrayExperiment",
          function(object) rowData(object)$isSnp)

#' @export
#' @aliases lrr,SnpArrayExperiment-method
#' @rdname LowLevelSummaries
#' @seealso \code{copyNumber}
setMethod("lrr", "SnpArrayExperiment",
          function(object) assays(object)[["cn"]])

#' @export
#' @aliases genotypes,SnpArrayExperiment-method
#' @rdname LowLevelSummaries
setMethod(genotypes, "SnpArrayExperiment",
          function(object) assays(object)[["genotypes"]])

setMethod(chromosome, "SnpArrayExperiment",
          function(object) as.character(seqnames(object)))

setMethod(position, "SnpArrayExperiment",
          function(object) start(object))

setMethod(genomeBuild, "SnpArrayExperiment",
          function(object) metadata(rowData(object))[["genome"]])


setMethod(distance, "SnpArrayExperiment", function(x) diff(start(x)))

setMethod(NA_index, "SnpArrayExperiment", function(x){
  b <- baf(x)
  r <- copyNumber(x)
  p <- start(x)
  bs <- rowSums(is.na(b))
  rs <- rowSums(is.na(r))
  i <- which(bs > 0 | rs > 0 | is.na(p))
  ##xx <- list(baf(x)[,1], copyNumber(x)[,1], start(x))
  ##i <- NA_index(xx)
})


setMethod(NA_filter, "SnpArrayExperiment", function(x, i){
  if(missing(i)) i <- NA_index(x)
  if(length(i) > 0) x <- x[-i, ]
  x
})

assayNames <- function(x) names(assays(x))
requiredAssays <- function() c("cn", "baf")
isCnAssays <- function(x) all(requiredAssays() %in% assayNames(x))

#' Create an assays object from log R ratios and B allele frequencies
#'
#'
#' @param cn matrix of log R ratios
#' @param baf matrix of B allele frequencies
#' @param ... one can pass additional matrices to the SnpExperiment
#' constructor, such as the SNP genotypes.
#' @export
snpArrayAssays <- function(cn=new("matrix"),
                           baf=new("matrix"), ...){
  assays <- GenomicRanges:::.ShallowSimpleListAssays(
    data = SimpleList(cn=cn,
      baf=baf, ...))
  assays
}



.calculateTransitionProbs <- function(x, param){
  p <- calculateTransitionProbability(start(x), param=param)
  ## transition probabilities between chromosomes
  chrom <- as.integer(as.factor(chromosome(x)))
  index <- which(diff(chrom)!=0)
  p[index] <- 1/6 ## any of states equally likely
  p
}

setMethod(calculateTransitionProbability, "SnpArrayExperiment",
          function(x, param=TransitionParam()){
            .calculateTransitionProbs(x, param)
          })

setMethod(calculateTransitionProbability, "oligoSnpSet",
          function(x, param=TransitionParam()){
            .calculateTransitionProbs(x, param)
          })



setMethod(isSnp, "SnpArrayExperiment",
          function(object) rowData(object)$isSnp)

setMethod(copyNumber, "SnpArrayExperiment",
          function(object) assays(object)[["cn"]])

setReplaceMethod("copyNumber", c("SnpArrayExperiment", "matrix"),
                 function(object, value) {
                   assays(object)[["cn"]] <- value
                 object
               })

setReplaceMethod("baf", c("SnpArrayExperiment", "matrix"),
                 function(object, value) {
                   assays(object)[["baf"]] <- value
                 object
               })

##  @describeIn SnpArrayExperiment-class
##  @aliases sweepMode,SnpArrayExperiment-method
setMethod("sweepMode", "SnpArrayExperiment",
          function(x, MARGIN){
            se <- x[chromosome(x) %in% autosomes(chromosome(x)), ]
            cn <- copyNumber(se)
            if(MARGIN==1){
              m <- rowModes(cn)
            }
            if(MARGIN==2){
              m <- colModes(cn)
            }
            ## subtract the autosomal model from all chromosomes
            cn <- sweep(copyNumber(x), MARGIN, m)
            copyNumber(x) <- cn
            x
          })


##setMethod(SnpArrayExperiment, "ShallowSimpleListAssays",
##          function(cn, baf, assays, rowData=GRanges(),
##                   colData=DataFrame(rownames=colnames(assays$data[[1]])),
##                   isSnp, ...){
##            if("isSnp" %in% colnames(mcols(rowData))){
##              rowData <- SnpGRanges(rowData)
##            } else rowData <- SnpGRanges(rowData, isSnp)
##            new("SnpArrayExperiment",
##                assays=assays,
##                rowData=rowData,
##                colData=colData, ...)
##          })
