setMethod(SnpArrayExperiment, "missing",
          function(cn, baf, rowData=GRanges(),
                   colData=DataFrame(), isSnp=logical(), ...){
            se <- new("SnpArrayExperiment", assays=snpArrayAssays(...),
                      rowData=SnpGRanges(), ...)
            se
          })

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

setValidity("SnpArrayExperiment", function(object){
  msg <- validObject(rowData(object))
  msg <- msg && (requiredAssays() %in% names(assays(object)))
  msg
})


setMethod(baf, "SnpArrayExperiment",
          function(object) assays(object)[["baf"]])

setMethod(chromosome, "SnpArrayExperiment",
          function(object) as.character(seqnames(object)))

setMethod(copyNumber, "SnpArrayExperiment",
          function(object) assays(object)[["cn"]])

setMethod(isSnp, "SnpArrayExperiment",
          function(object) rowData(object)$isSnp)

setMethod(lrr, "SnpArrayExperiment",
          function(object) assays(object)[["cn"]])

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
isCnAssays <- function(x) all(requiredAssays() %in% assayNames(object))

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
  p[index] <- 0.801 ## lowest value possible
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

##
## setMethod(SnpArrayExperiment, "ShallowSimpleListAssays",
##           function(cn, baf, assays, rowData=GRanges(),
##                    colData=DataFrame(rownames=colnames(assays$data[[1]])),
##                    isSnp, ...){
##             if("isSnp" %in% colnames(mcols(rowData))){
##               rowData <- SnpGRanges(rowData)
##             } else rowData <- SnpGRanges(rowData, isSnp)
##             new("SnpArrayExperiment",
##                 assays=assays,
##                 rowData=rowData,
##                 colData=colData, ...)
##           })

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
