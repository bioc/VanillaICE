setClass("HmmOptionList", contains="list")


setClass("Viterbi",
         representation(state="integer",
                        loglik="numeric",
                        forward_backward="matrix"))


setClass("EmissionParam",
         representation=representation(
           cn_means="numeric",
           cn_sds="numeric",
           baf_means="numeric",
           baf_sds="numeric",
           initial="numeric",
           EMupdates="integer",
           CN_range="numeric",
           proportionOutlier="numeric",
           temper="numeric",
           p_outlier="numeric",
           modelHomozygousRegions="logical"))



setClass("TransitionParam",
         representation(
           taup="numeric",
           taumax="numeric"))

#' Classes and methods for storing/getting log-likelihoods from Viterbi algorithm
#'
#' Exported for internal use by other BioC packages
#' @slot loglik a numeric vector
#' @slot tolerance a numeric vector
#' @docType class
#' @rdname LogLik-class
#' @seealso \code{\link{LogLik}}
#' @export
setClass("LogLik", representation(loglik="numeric", tolerance="numeric"))


#' A class allowing matrix or NULL objects
#'
#' Exported for internal use by other BioC packages
#'
#' @rdname matrixOrNULL-class
#' @aliases matrixOrNULL-class
#' @name matrixOrNULL
#' @docType class
NULL
#' @export
setClassUnion("matrixOrNULL", c("matrix", "NULL"))

setClass("HmmParam",
         representation(
           emission="matrixOrNULL",
           emission_param="EmissionParam",
           transition="numeric",
           chromosome="character",
           loglik="LogLik",
           viterbi="Viterbi",
           compute_posteriors="logical",
           verbose="logical"))


setClass("SnpDataFrame", contains="DataFrame")

#' An extension to GRanges for representing SNPs
#'
#' @slot elementMetadata a \code{SnpDataFrame}
#' @rdname SnpGRanges
#' @export
setClass("SnpGRanges", contains="GRanges",
         representation(elementMetadata="SnpDataFrame"))


#' A RangedSummarizedExperiment-derived class of marker-level SNP array data
#' for copy number inference
#'
#' @rdname SnpArrayExperiment-class
#' @export
setClass("SnpArrayExperiment",
         contains = "RangedSummarizedExperiment",
         representation(rowRanges="SnpGRanges"))

# # @aliases HmmGRanges-class
# # @rdname HmmGRanges-class
# # @export
setClass("HmmGRanges", contains="GRanges", representation(emission_param="EmissionParam"))

#' Container for the common criteria used to  filtering genomic ranges
#'
#' The maximum a posteriori estimate of the trio copy number state for
#' each genomic range is represented in a
#' \code{\link{GRanges}}-derived class.  Ultimately, these ranges will
#' be filtered based on the trio copy number state (e.g., denovo
#' deletions), size, number of features (SNPs), or chromosome.
#' \code{FilterParam} is a container for the parameters commmonly used
#' to filter the genomic ranges.
#' @slot probability a length-one numeric vector indicating the
#' minimum posterior probability for the called state. Genomic
#' intervals with posterior probabilities below \code{probability}
#' will be filtered.
#' @slot numberFeatures a positive integer indicating the minimum
#' number of features in a segment
#' @slot seqnames a character vector of \code{seqnames} to select
#' (i.e., 'chr1' for only those intervals on chromosome 1)
#' @slot width positive integer indicating the minimal width of
#' genomic intervals
#' @slot state character string indicating which hidden Markov model states to select
#' @examples
#' fp <- FilterParam()
#' width(fp)
#' numberFeatures(fp)
#' seqnames(fp)
#' @rdname FilterParam-class
#' @export
setClass("FilterParam", representation(probability="numeric",
                                       numberFeatures="numeric",
                                       seqnames="character",
                                       width="integer",
                                       state="character"))

#' Container for the segmented data and the 6-state HMM model parameters
#'
#' @slot granges a \code{GRanges} object
#' @slot param a \code{HmmParam} object
#' @slot posterior a matrix of posterior probabilities
#' @slot filters a \code{FilterParam} object
#' @seealso \code{\link{hmm2}}
#' @rdname HMM
#' @examples
#' data(snp_exp)
#' hmm_list <- hmm2(snp_exp[,1])
#' resultsFirstSample <- hmm_list[[1]]
#' resultsFirstSample
#' @export
setClass("HMM", representation(granges="GRanges",
                               param="HmmParam",
                               posterior="matrix",
                               filters="FilterParam"))

#' Class, constructor, and methods for representing HMM results from
#' multiple samples
#'
#' Each element of the HMMList contains the genomic intervals of the
#' HMM segmentation (GRanges-derived object), parameters from the
#' Baum-Welch, and a \code{FilterParam} object.
#'
#' @slot .Data a list. Each element of the list should be a \code{HMM} object.
#' @rdname HMMList-class
#' @seealso \code{\linkS4class{HMM}}
#' @examples
#' data(snp_exp)
#' fit <- hmm2(snp_exp)
#' class(fit)
#' identical(length(fit), ncol(snp_exp))
#' unlist(fit)
#' @export
setClass("HMMList", contains="list")

#' Parameters for parsing source files containing SNP-array processed
#' data, such as GenomeStudio files for the Illumina platform
#'
#' Raw SNP array processed files have headers and variable labels that
#' may depend the software, how the output files was saved, the
#' software version, and other factors.  The purpose of this container
#' is to collect the parameters relevant for reading in the source
#' files for a particular project in a single container.  This may
#' require some experimentation as the example illustrates.  The
#' function \code{\link{fread}} in the \code{data.table} package
#' greatly simplifies this process.
#' 
#' @slot index_genome an integer vector
#' @slot cnvar the column label for the log R ratios
#' @slot bafvar the column label for the B allele frequencies
#' @slot gtvar the column label(s) for the genotypes
#' @slot scale length-one numeric vector indicating how the low-level data
#'   should be scaled prior to saving on disk
#' @slot select numeric vector indicating which columns to read
#' @slot row.names length-one numeric vector indicating which column
#' the SNP names are in
#'
#' @rdname CopyNumScanParams
#' @export
setClass("CopyNumScanParams", representation(index_genome="integer",
                                             cnvar="character",
                                             bafvar="character",
                                             gtvar="character",
                                             scale="numeric",
                                             select="integer",
                                             row.names="integer"))

#' ArrayViews class, constructor, and methods
#'
#'
#' ArrayViews provides views to the low-level data -- log R ratios, B
#' allele frequencies, and genotypes that are stored in parsed files
#' on disk, often scaled and coerced to an integer.  Accessors to the
#' low-level data are provided that extract the marker-level summaries
#' from disk, rescaling when appropriate.
#'
#' @slot colData  A character string
#' @slot rowRanges A \code{DataFrame}. WARNING: The accessor for this slot is \code{rowRanges}, not \code{rowRanges}!
#' @slot index A \code{GRanges} object
#' @slot sourcePaths A character string providing complete path to source files (one file per sample) containing low-level summaries (Log R ratios, B allele frequencies, genotypes)
#' @slot scale A length-one numeric vector
#' @slot parsedPath A character string providing full path to where
#' parsed files should be saved
#' @slot lrrFiles character vector of filenames for log R ratios
#' @slot bafFiles character vector of filenames for BAFs
#' @slot gtFiles character vector of filenames for genotypes
#' @aliases '[',ArrayViews,ANY-method baf_means,ArrayViews-method
#' @rdname ArrayViews-class
#' @export
setClass("ArrayViews",
         representation(colData="DataFrame",
                        rowRanges="GRanges",
                        index="integer",
                        sourcePaths="character",
                        scale="numeric",
                        parsedPath="character",
                        lrrFiles="character",
                        bafFiles="character",
                        gtFiles="character"))

setClass("HmmTrellisParam", representation(expandfun="function",
                                           ylimits="list"))

#' Paramater class for plotting idiograms
#'
#' @slot seqnames length-one character vector providing chromosome name
#' @slot seqlengths length-one numeric vector indicating size of chromosome
#' @slot unit character string indicating unit for genomic position (default is 'kb')
#' @slot genome character string indicating genome build
#' @slot box a list of parameters for plotting the box around the part of the idiogram that is plotted.
#' @export
#' @rdname IdiogramParams-class
#' @examples
#' if(require(BSgenome.Hsapiens.UCSC.hg18) && require(grid)){
#'    si <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)
#'    iparam <- IdiogramParams(seqnames="chr1",
#'                             genome="hg18",
#'                             seqlengths=seqlengths(si)["chr1"],
#'                             box=list(xlim=c(20e6L, 25e6L), color="blue", lwd=2))
#'    iparam
#'    idiogram <- plot(iparam)
#'    vp <- viewport(x=0.05, y=0.8, width=unit(0.9, "npc"), height=unit(0.2, "npc"),
#'                   name="vp1", just=c("left", "bottom"))
#'    grid.newpage()
#'    pushViewport(vp)
#'    print(idiogram, vp=vp, newpage=FALSE)
#' }
#'
setClass("IdiogramParams", representation(seqnames="character",
                                          seqlengths="numeric",
                                          unit="character",
                                          genome="character",
                                          box="list"))



setGeneric("plot")
