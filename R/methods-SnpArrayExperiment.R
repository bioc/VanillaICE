#' @examples
#' ## empty container
#' SnpArrayExperiment()
#'
#' data(snp_exp) # example
#'
#' SnpArrayExperiment(cn=lrr(snp_exp), baf=baf(snp_exp),
#'                    rowRanges=rowRanges(snp_exp))
#' @aliases SnpArrayExperiment,missing-method
#' @rdname SnpArrayExperiment-class
setMethod(SnpArrayExperiment, "missing",
          function(cn, baf, rowRanges=GRanges(),
                   colData=DataFrame(), isSnp=logical(), ...){
            new("SnpArrayExperiment", SummarizedExperiment(
                assays=snpArrayAssays(...)$data,
                rowRanges=SnpGRanges(),
                ...
            ))
          })


#' @aliases SnpArrayExperiment,matrix-method
#' @rdname SnpArrayExperiment-class
setMethod(SnpArrayExperiment, "matrix",
          function(cn, baf, rowRanges=GRanges(),
                   colData=DataFrame(row.names=colnames(cn)),
                   isSnp=logical(), ...){
            assays <- snpArrayAssays(cn=cn, baf=baf, ...)
            if("isSnp" %in% colnames(mcols(rowRanges))){
              rowRanges <- SnpGRanges(rowRanges)
            } else {
              if(length(isSnp) != length(rowRanges)) stop(" isSnp must be the same length as rowRanges")
              rowRanges <- SnpGRanges(rowRanges, isSnp)
            }
            new("SnpArrayExperiment", SummarizedExperiment(
                assays=assays$data,
                rowRanges=rowRanges,
                colData=colData
            ))
          })


setAs("CNSet", "SnpArrayExperiment",
      function(from){
        b.r <- calculateRBaf(from)
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
        rowranges <- SnpGRanges(g, isSnp=isSnp(from)[i])
        rowranges <- SnpGRanges(rowranges)
        names(rowranges) <- nms
        ##tmp=snpArrayAssays(cn=lrrs/100, baf=bafs/100)
        ##new("SnpArrayExperiment", assays=tmp, rowRanges=rowranges, colData=DataFrame(pData(from)))
        se <- SnpArrayExperiment(cn=lrrs/100, baf=bafs/1000,
                                 rowRanges=SnpGRanges(rowranges),
                                 colData=DataFrame(pData(from)))
        sort(se)
      })

setAs("oligoSnpSet", "SnpArrayExperiment",
      function(from){
        rowranges <- as(featureData(from), "SnpGRanges")
        coldata <- as(as(phenoData(from), "data.frame"), "DataFrame")
        if(!"baf" %in% assayDataElementNames(from)) stop("BAFs not available")
        cn_assays <- snpArrayAssays(cn=copyNumber(from),
                                    baf=baf(from),
                                    gt=calls(from))
        new("SnpArrayExperiment", SummarizedExperiment(
            assays=cn_assays$data,
            rowRanges=rowranges,
            colData=coldata
        ))
      })



setValidity("SnpArrayExperiment", function(object){
  msg <- validObject(rowRanges(object))
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
          function(object) rowRanges(object)$isSnp)

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
          function(object) metadata(rowRanges(object))[["genome"]])


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

.ShallowSimpleListAssays <-
    function(..., data = SimpleList())
{
    xx <- SummarizedExperiment:::.ShallowSimpleListAssays0(...)
    xx$data <- data
    xx
}

#' Create an assays object from log R ratios and B allele frequencies
#'
#'
#' This function is exported primarily for internal use by other BioC
#' packages.
#' @param cn matrix of log R ratios
#' @param baf matrix of B allele frequencies
#' @param ... additional matrices of the same dimension, such as SNP genotypes.
#' @examples
#' data(snp_exp)
#' r <- lrr(snp_exp)
#' b <- baf(snp_exp)
#' sl <- snpArrayAssays(cn=r, baf=b)
#' @export
snpArrayAssays <- function(cn=new("matrix"),
                           baf=new("matrix"), ...){
  assays <- .ShallowSimpleListAssays(
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
          function(object) rowRanges(object)$isSnp)

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


#' @aliases sweepMode sweepMode,SnpArrayExperiment-method
#' @rdname sweepMode
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

## Inclusion of this function means that BSgenome must be imported and not
## listed in Suggests

#' Create an example SnpArrayExperiment from source files containing
#' marker-level genomic data that are provided in this package
#'
#' @return A \code{\link{SnpArrayExperiment}}
#' @param bsgenome a \code{BSgenome} object
#' @export
#' @examples
#' \dontrun{
#'    if(require("BSgenome.Hsapiens.UCSC.hg18")){
#'      genome <- BSgenome.Hsapiens.UCSC.hg18
#'      snp_exp <- getExampleSnpExperiment(genome)
#'    }
#' }
getExampleSnpExperiment <- function(bsgenome){
  extdir <- system.file("extdata", package="VanillaICE", mustWork=TRUE)
  features <- suppressWarnings(fread(file.path(extdir, "SNP_info.csv")))
  fgr <- GRanges(paste0("chr", features$Chr), IRanges(features$Position, width=1),
                 isSnp=features[["Intensity Only"]]==0)
  fgr <- SnpGRanges(fgr)
  names(fgr) <- features[["Name"]]
  seqlevels(fgr) <- seqlevels(bsgenome)[seqlevels(bsgenome) %in% seqlevels(fgr)]
  seqinfo(fgr) <- seqinfo(bsgenome)[seqlevels(fgr),]
  fgr <- sort(fgr)
  file <- list.files(extdir, full.names=TRUE, recursive=TRUE, pattern="FinalReport")[5]
  dat <- fread(file)
  select <- match(c("SNP Name", "Allele1 - AB", "Allele2 - AB", "Log R Ratio", "B Allele Freq"), names(dat))
  index_genome <- match(names(fgr), dat[["SNP Name"]])
  scan_params <- CopyNumScanParams(index_genome=index_genome, select=select)
  views <- ArrayViews(rowRanges=fgr, sourcePaths=file)
  parseSourceFile(views, param=scan_params)
  SnpExperiment(views)
}


setMethod("isAutosome", "SnpArrayExperiment", function(object){
  isAutosome(rowRanges(object))
})


#' @aliases isHeterozygous,SnpArrayExperiment-method
#' @rdname isHeterozygous
setMethod("isHeterozygous", "SnpArrayExperiment", function(object, cutoff){
  isHeterozygous(baf(object), cutoff)
})

#' @aliases isHeterozygous,numeric-method
#' @rdname isHeterozygous
setMethod("isHeterozygous", "numeric", function(object, cutoff){
  object > cutoff[1] & object < cutoff[2] & !is.na(object)
})

#' @aliases isHeterozygous,matrix-method
#' @rdname isHeterozygous
setMethod("isHeterozygous", "matrix", function(object, cutoff){
  x <- as.numeric(object)
  ishet <- isHeterozygous(x, cutoff)
  matrix(ishet, nrow(object), ncol(object))
})
