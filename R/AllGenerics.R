#' Fit a 6-state HMM to log R ratios and B allele frequencies
#' estimated from SNP arrays
#'
#' This function is intended for estimating the integer copy number
#' from germline or DNA of clonal origin using a 6-state HMM.  The
#' states are homozygous deletion, hemizygous deletion, diploid copy
#' number, diploid region of homozygosity, single copy gain, and two+
#' copy gain.  Because heterozygous markers are more informative for
#' copy number than homozygous markers and regions of homozgosity are
#' common in normal genomes, we currently computed a weighted average
#' of the BAF emission matrix with a uniform 0,1 distribution by the
#' probability that the marker is heterozygous, thereby downweighting
#' the contribution of homozygous SNPs to the likelihood.  In addition
#' to making the detection of copy-neutral regions of homozgosity less
#' likely, it also helps prevent confusing hemizygous deletions with
#' copy neutral regions of homozygosity -- the former would be driven
#' mostly by the log R ratios.  This is experimental and subject to
#' change.
#'
#' The \code{hmm2} method allows parallelization across samples using
#' the foreach paradigm.  Parallelization is automatic when enabled
#' via packages such as snow/doSNOW.
#'
#'
#' @param object A \code{\link{SnpArrayExperiment}}
#' @param emission_param A \code{\link{EmissionParam}} object
#' @param transition_param A \code{\link{TransitionParam}} object
#' @param ... currently ignored
#' @examples
#' tp <- TransitionParam()
#' TransitionParam(taup=1e12)
#' @export
setGeneric("hmm2", function(object, emission_param=EmissionParam(),
                            transition_param=TransitionParam(), ...) standardGeneric("hmm2"))

setGeneric("loglik", function(object) standardGeneric("loglik"))

#' Calculate the emission probabilities for the 6-state HMM
#'
#' Given the data and an object containing parameters for the HMM,
#' this function computes emission probabilities.  This function is
#' not intended to be called by the user and is exported for internal
#' use by other BioC packages.
#' @return A matrix of emission probabilities. Column correspond to
#' the HMM states and rows correspond to markers on the array (SNPs
#' and nonpolymorphic markers)
#'
#' @param x list of low-level data with two elements: a numeric vector of log
#' R ratios and a numeric vector of B allele frequencies
#' @param param  parameters for the 6-state HMM
#' @seealso baumWelchUpdate
#' @export
#' @aliases calculateEmission,RangedSummarizedExperiment-method calculateEmission,list-method calculateEmission,numeric-method
setGeneric("calculateEmission", function(x, param=EmissionParam()) standardGeneric("calculateEmission"))

setGeneric("forward_backward", function(x) standardGeneric("forward_backward"))

setGeneric("Viterbi",
           function(state=Rle(), loglik=numeric(1),
                    forward_backward=matrix(0.0, length(state), 6))
           standardGeneric("Viterbi"))

setGeneric("transition", function(object) standardGeneric("transition"))

setGeneric("initial", function(object) standardGeneric("initial"))

setGeneric("statei", function(object) standardGeneric("statei"))

setGeneric("state<-", function(object,value) standardGeneric("state<-"))

setGeneric("statef", function(object) standardGeneric("statef"))

#' A parameter class for computing Emission probabilities
#'
#' Parameters for computing emission probabilities for a 6-state HMM,
#' including starting values for the mean and standard deviations for
#' log R ratios (assumed to be Gaussian) and B allele frequencies
#' (truncated Gaussian), and initial state probabilities.
#'
#' The log R ratios are assumed to be emitted from a normal
#' distribution with a mean and standard deviation that depend on the
#' latent copy number.  Similarly, the BAFs are assumed to be emitted
#' from a truncated normal distribution with a mean and standard
#' deviation that depends on the latent number of B alleles relative
#' to the total number of alleles (A+B).
#'
#' @return numeric vector
#' @examples
#' ep <- EmissionParam()
#' cn_means(ep)
#' @rdname EmissionParam-methods
#' @aliases cn_means,EmissionParam-method cn_means,HmmParam-method cn_means<-,EmissionParam,numeric-method
#' @export
setGeneric("cn_means", function(object) standardGeneric("cn_means"))

#' @examples
#' ep <- EmissionParam()
#' cn_sds(ep)
#' @rdname EmissionParam-methods
#' @aliases cn_sds,EmissionParam-method cn_sds,HmmParam-method cn_sds<-,EmissionParam,numeric-method
#' @export
setGeneric("cn_sds", function(object) standardGeneric("cn_sds"))


#' @examples
#' ep <- EmissionParam()
#' baf_means(ep)
#' @rdname EmissionParam-methods
#' @aliases baf_means,EmissionParam-method baf_means,HmmParam-method
#' @export
setGeneric("baf_means", function(object) standardGeneric("baf_means"))

#' @examples
#' ep <- EmissionParam()
#' baf_sds(ep)
#' @rdname EmissionParam-methods
#' @aliases baf_sds,EmissionParam-method baf_sds,HmmParam-method
#' @export
setGeneric("baf_sds", function(object) standardGeneric("baf_sds"))

#' @examples
#' ep <- EmissionParam()
#' baf_means(ep) <- baf_means(ep)
#' @rdname EmissionParam-methods
#' @aliases baf_means<-,EmissionParam,numeric-method
#' @export
setGeneric("baf_means<-", function(object, value) standardGeneric("baf_means<-"))

#' @examples
#' ep <- EmissionParam()
#' baf_sds(ep) <- baf_sds(ep)
#' @param value numeric vector
#' @rdname EmissionParam-methods
#' @aliases baf_sds<-,EmissionParam,numeric-method
#' @export
setGeneric("baf_sds<-", function(object, value) standardGeneric("baf_sds<-"))

#' @examples
#' ep <- EmissionParam()
#' cn_sds(ep) <- cn_sds(ep)
#' @rdname EmissionParam-methods
#' @export
setGeneric("cn_sds<-", function(object, value) standardGeneric("cn_sds<-"))

#' @examples
#' ep <- EmissionParam()
#' cn_means(ep) <- cn_means(ep)
#' @rdname EmissionParam-methods
#' @export
setGeneric("cn_means<-", function(object, value) standardGeneric("cn_means<-"))

setGeneric("taup", function(object) standardGeneric("taup"))

setGeneric("taumax", function(object) standardGeneric("taumax"))

#' Constructor for EmissionParam class
#'
#'
#' @param cn_means numeric vector of starting values for log R ratio means (order is by copy number state)
#' @param cn_sds numeric vector of starting values for log R ratio standard deviations (order is by copy number state)
#' @param baf_means numeric vector of starting values for BAF means ordered.  See example for details on how these are ordered.
#' @param baf_sds numeric vector of starting values for BAF means ordered.  See example for details on how these are ordered.
#' @param initial numeric vector of intial state probabilities
#' @param EMupdates number of EM updates
#' @param CN_range the allowable range of log R ratios.  Log R ratios outside this range are thresholded.
#' @param temper Emission probabilities can be tempered by emit^temper. This is highly experimental.
#' @param p_outlier probability that an observation is an outlier (assumed to be the same for all markers)
#' @param modelHomozygousRegions logical.   If FALSE (default), the
#' emission probabilities for BAFs are modeled from a mixture of
#' truncated normals and a Unif(0,1) where the mixture probabilities
#' are given by the probability that the SNP is heterozygous. See
#' Details below for a discussion of the implications.
#'
#' @section Details:
#' When \code{modelHomozygousRegions} is FALSE (the default in
#' versions >= 1.28.0), emission probabilities for B allele frequences
#' are calculated from a mixture of a truncated normal densities and a
#' Unif(0,1) density with the mixture probabilities given by the
#' probability that a SNP is homozygous.  In particular, let \code{p}
#' denote a 6 dimensional vector of density estimates from a truncated
#' normal distribution for the latent genotypes 'A', 'B', 'AB', 'AAB',
#' 'ABB', 'AAAB', and 'ABBB'.  The probability that a genotype is
#' homozygous is estimated as
#'
#' \deqn{prHom=(p["A"]  + p["B"])/sum(p)}
#'
#' and the probability that the genotype is heterozygous (any latent
#' genotype that is not 'A' or 'B') is given by
#'
#' \deqn{prHet = 1-prHom}
#'
#' Since the density of a Unif(0,1) is 1, the 6-dimensional vector of
#' emission probability at a SNP is given by
#'
#' \deqn{emit = prHet * p + (1-prHet)}
#'
#' The above has the effect of minimizing the influence of BAFs near 0
#' and 1 on the state path estimated by the Viterbi algorithm. In
#' particular, the emission probability at homozygous SNPs will be
#' virtually the same for states 3 and 4, but at heterozygous SNPs the
#' emission probability for state 3 will be an order of magnitude
#' greater for state 3 (diploid) compared to state 4 (diploid region
#' of homozygosity).  The advantage of this parameterization are fewer
#' false positive hemizygous deletion calls.  [ Log R ratios tend to
#' be more sensitive to technical sources of variation than the
#' corresponding BAFs/ genotypes.  Regions in which the log R ratios
#' are low due to technical sources of variation will be less likely
#' to be interpreted as evidence of copy number loss if heterozygous
#' genotypes have more 'weight' in the emission estimates than
#' homozgous genotypes.  ]  The trade-off is that only states
#' estimated by the HMM are those with copy number alterations.  In
#' particular, copy-neutral regions of homozygosity will not be
#' called.
#'
#' By setting \code{modelHomozygousRegions = TRUE}, the emission
#' probabilities at a SNP are given simply by the \code{p} vector
#' described above and copy-neutral regions of homozygosity will be
#' called.#'
#'
#' @examples
#' ep <- EmissionParam()
#' show(ep)
#' cn_means(ep)
#' cn_sds(ep)
#' baf_means(ep)
#' baf_sds(ep)
#' @family EmissionParam-methods
#' @rdname EmissionParam-methods
#' @aliases EmissionParam EmissionParam,missing-method EmissionParam,numeric-method
#' @export
setGeneric("EmissionParam", function(cn_means=CN_MEANS(),
                                     cn_sds=CN_SDS(),
                                     baf_means=BAF_MEANS(),
                                     baf_sds=BAF_SDS(),
                                     initial=rep(1/6, 6),
                                     EMupdates=5L,
                                     CN_range=c(-5, 3),
                                     temper=1,
                                     p_outlier=1/100,
                                     modelHomozygousRegions=FALSE)
           standardGeneric("EmissionParam"))

#' Accessor for the maximum number of Baum-Welch updates
#'
#' This function is exported primarily for internal use by other BioC
#' packages.
#'
#' @param object see \code{showMethods("EMupdates")}
#' @aliases EMupdates,EmissionParam-method EMupdates,HmmParam-method
#' @rdname EmissionParam-methods
#' @export
setGeneric("EMupdates", function(object) standardGeneric("EMupdates"))

setGeneric("EMupdates<-", function(object,value) standardGeneric("EMupdates<-"))

setGeneric("probOutlier", function(object) standardGeneric("probOutlier"))

#' Methods to set and get emission probabilities
#'
#' Get or set a matrix of emission probabilities. This function is
#' exported primarily for internal use by other BioC packages.
#'
#' @return matrix
#' @param object  see \code{showMethods(emission)}
#' @export
#' @rdname emission
#' @aliases emission,HmmParam-method
setGeneric("emission", function(object) standardGeneric("emission"))

## used by MinimumDistance

#' @param value a matrix of emission probabilities
#' @rdname emission
#' @aliases emission<-,HmmParam-method emission<-,HMM-method
#' @export
setGeneric("emission<-", function(object, value) standardGeneric("emission<-"))

#' Constructor for TransitionParam class
#'
#' Contains parameters for computing transition probabilities
#'
#' @examples
#' TransitionParam()
#' ## higher values of taup make transitions between states less likely
#' TransitionParam(taup=1e12)
#' @param taup length-one numeric vector
#' @param taumax The maximum probability that the current state is the
#' same as the preceding state. See details
#'
#' @details Diagonal elements of the transition probability matrix are
#' computed as e^{-2*d/taup}, where d is the distance between markers
#' i and i-1 and \code{taup} is typically in the range of 1xe10.  This
#' probability is constrained to be no larger than \code{taumax}. The
#' probabilities on the off-diagonal elements are the same and are
#' subject to the constraint that the rows of the transition
#' probability matrix sum to 1.
#' @aliases TransitionParam TransitionParam,missing-method TransitionParam,numeric-method
#' @export
setGeneric("TransitionParam",
           function(taup=1e10,
                    taumax=1-5e6)
           standardGeneric("TransitionParam"))

#' Constructor for HmmParam class
#'
#' Contains emission probabilities, parameters for emission
#' probabilities, and transition probabilities required for computing
#' the most likely state path via the Viterbi algorithm
#'
#' @param emission A matrix of emission probabilities
#' @param emission_param an object of class \code{EmissionParam}
#' @param transition vector of transition probabilities whose length
#' is N-1, where N is the number of markers.  User should provide the
#' probability that the state at marker j is the same as the state at
#' marker j-1.  It is assumed that the probability of transitioning to
#' state_j from state_j-1 is the same for all states != state_j-1.
#' @param chromosome character vector
#' @param loglik an object of class \code{LogLik}
#' @param viterbi an object of class \code{Viterbi}
#' @param compute_posteriors logical
#' @param verbose logical
#' @examples
#' HmmParam()
#' @aliases HmmParam HmmParam,missing-method HmmParam,matrix-method
#' @rdname HmmParam
#' @export
setGeneric("HmmParam",
           function(emission=matrix(0.0, 0, 0),
                    ##initial=rep(1/6, 6),
                    emission_param=EmissionParam(),
                    transition=rep(0.99, nrow(emission)),
                    chromosome=character(nrow(emission)),
                    loglik=LogLik(),
                    viterbi=Viterbi(),
                    compute_posteriors=TRUE,
                    verbose=FALSE)
           standardGeneric("HmmParam"))

setGeneric("calculateTransitionProbability",
           function(x, param=TransitionParam())
           standardGeneric("calculateTransitionProbability"))

setGeneric("modev", function(x) standardGeneric("modev"))

#' Constructor for SnpGRanges class
#'
#' @param object A \code{GRanges} object
#' @param isSnp A logical vector.  Each genomic interval in the
#' \code{GRanges} container corresponds to a marker on the genotyping
#' array.  \code{isSnp} is FALSE for nonpolymorphic markers such as
#' those included on the Affymetrix 6.0 chips.
#' @param ... ignored
#' @examples
#' SnpGRanges()
#' g <- GRanges("chr1", IRanges(15L, 15L))
#' SnpGRanges(g, isSnp=TRUE)
#' @export
#' @rdname SnpGRanges
#' @aliases SnpGRanges
setGeneric("SnpGRanges", function(object=GRanges(),
                                  isSnp, ...)
           standardGeneric("SnpGRanges"))


#' Remove SNPs with NAs in any of the low-level estimates
#'
#' @return An object of the same class
#' @param x a container for SNP data (\code{\link{SnpArrayExperiment}})
#' @param i integer vector to subset
#' @aliases NA_filter NA_filter,SnpArrayExperiment-method NA_filter,character-method NA_filter,list-method NA_filter,numeric-method NA_filter,oligoSnpSet-method
#' @export
setGeneric("NA_filter", function(x, i)
           standardGeneric("NA_filter"))

setGeneric("distance", function(x) standardGeneric("distance"))


#' SnpArrayExperiment container for log R ratios and B allele frequencies
#'
#' Constructor for SnpArrayExperiment
#'
#' @param cn matrix of copy number estimates (e.g., log R ratios)
#' @param baf  matrix of B allele frequencies
#' @param rowRanges GRanges object for SNPs/nonpolymorphic markers
#' @param colData DataFrame containing sample-level covariates
#' @param isSnp  logical vector indicating whether marker is a SNP
#' @param ... additional arguments passed to \code{SummarizedExperiment()} constructor function
#' @docType methods
#' @rdname SnpArrayExperiment-class
#' @export
setGeneric("SnpArrayExperiment", function(cn,
                                          baf,
                                          rowRanges=GRanges(),
                                          colData=DataFrame(),
                                          isSnp=logical(), ...)
           standardGeneric("SnpArrayExperiment"))

#' Sweep the modal log R ratio (by row or column) from a matrix of log
#' R ratios
#'
#' This function simplifies the process of sweeping the modal log R
#' ratio from the rows or columns of a \code{SnpArrayExperiment}
#' object.  It is most useful when a large number of samples (more
#' than 10) are available and the dataset is a collection of germline
#' samples.  We assume that the samples are from a single batch and
#' that the modal value will be a robust estimate of the mean log R
#' ratio for diploid copy number.  Variation in the modal estimates
#' between markers is presumed to be attributable to probe effects
#' (e.g., differences hybridization efficiency/PCR do to sequence
#' composition). For sex chromosomes, one should apply this function
#' separately to men and women and then recenter the resulting matrix
#' according to the expected copy number.
#'
#' @param x see \code{showMethods(sweepMode)}
#' @param MARGIN integer indicating which margin (1=rows, 2=columns)
#' to sweep the mode
#' @examples
#' data(snp_exp)
#' snp_exp_rowcentered <- sweepMode(snp_exp, 1)
#' snp_exp_colcentered <- sweepMode(snp_exp, 2)
#' x <- lrr(snp_exp)
#' x_rowcentered <- sweep(x, 1, rowModes(x))
#' all.equal(lrr(snp_exp_rowcentered), x_rowcentered)
#' @export
#' @return an object of the same class as \code{x}
setGeneric("sweepMode", function(x, MARGIN) standardGeneric("sweepMode"))

setGeneric("NA_index", function(x) standardGeneric("NA_index"))

## #' HmmGRanges container
## #
## #' @param states  copy number number state inferred by HMM
## #' @param feature_starts start location in reference genome [basepairs]
## #' @param feature_chrom  end location in reference genome [basepairs]
## #' @param loglik  the log likelihood
## #' @param emission_param an instance of \code{EmissionParam} class
## #' @examples
## #'  library(oligoClasses)
## #'  library(IRanges)
## #'  load(system.file("extdata", "snp_exp.rda", package="VanillaICE"))
## #'  states <- Rle(factor(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3)),
## #'                as.integer(c(996, 102, 902, 50, 2467, 102, 76, 1822,
## #'                             99, 900, 20, 160)))
## #'  hgr <- HmmGRanges(states=states, feature_starts=start(snp_exp),
## #'                    feature_chrom=chromosome(snp_exp), loglik=15.3)
## #'
## # @aliases HmmGRangeso
## # @rdname HmmGRanges-class
## # @export
setGeneric("HmmGRanges", function(states, feature_starts,
                                  feature_chrom, loglik, emission_param=EmissionParam())
           standardGeneric("HmmGRanges"))

##setGeneric("isSnp", function(object) standardGeneric("isSnp"))

setGeneric("loglikRatio", function(object) standardGeneric("loglikRatio"))

setGeneric("tolerance", function(object) standardGeneric("tolerance"))

setGeneric("loglik", function(object) standardGeneric("loglik"))
setGeneric("loglikLast", function(object) standardGeneric("loglikLast"))

setGeneric("loglik<-", function(object,value) standardGeneric("loglik<-"))

setGeneric("exceedsTolerance", function(object) standardGeneric("exceedsTolerance"))

#' Accessor for parameters used to compute emission probabilities
#'
#' Parameters for computing emission probabilities include the
#' starting values for the Baum Welch update and initial state
#' probabilities.
#' @examples
#' hparam <- HmmParam()
#' emissionParam(hparam)
#' ep <- EmissionParam()
#' cn_means(ep) <- log2(c(.1/2, 1/2, 2/2, 2/2, 3/2, 4/2))
#' emissionParam(hparam) <- ep
#' @param object an object of class \code{EmissionParam}
#' @return \code{\link{EmissionParam}} instance
#' @rdname emissionParam
#' @aliases emissionParam,HmmGRanges-method emissionParam,HmmParam-method emissionParam<-,HmmGRanges,EmissionParam-method emissionParam<-,HmmParam,EmissionParam-method emissionParam,HMM-method
#' @export
setGeneric("emissionParam", function(object) standardGeneric("emissionParam"))

#' @param value an object of class \code{EmissionParam}
#' @rdname emissionParam
#' @export
setGeneric("emissionParam<-", function(object,value) standardGeneric("emissionParam<-"))

setGeneric("viterbi", function(object) standardGeneric("viterbi"))

setGeneric("viterbi<-", function(object,value) standardGeneric("viterbi<-"))

setGeneric("verbose", function(object) standardGeneric("verbose"))

setGeneric("CN_range", function(object) standardGeneric("CN_range"))
setGeneric("proportionOutlier", function(object) standardGeneric("proportionOutlier"))
setGeneric("temper", function(object) standardGeneric("temper"))

setGeneric("scale")

setGeneric("gstudioPaths", function(object) standardGeneric("gstudioPaths"))
setGeneric("indexGenome", function(object) standardGeneric("indexGenome"))
setGeneric("read", function(object) standardGeneric("read"))

#' Accessor for SNP genotypes
#'
#' Extract SNP genotypes. Genotypes are assumed to be represented as
#' integers: 1=AA, 2=AB, 3=BB.
#' @rdname LowLevelSummaries
#' @param object see \code{showMethods("genotypes")}
setGeneric("genotypes", function(object) standardGeneric("genotypes"))

#' Constructor for SnpArrayExperiment
#'
#' A single-argument generic function to construct a SnpArrayExperiment.
#'
#' @param object see \code{showMethods('SnpExperiment')} for a list of supported objects
#' @export
setGeneric("SnpExperiment", function(object) standardGeneric("SnpExperiment"))
setGeneric("scaleBy", function(x, by) standardGeneric("scaleBy"))
setGeneric("scaleRead", function(x, params) standardGeneric("scaleRead"))
setGeneric("selectCols", function(object) standardGeneric("selectCols"))

#' Complete path to directory for keeping parsed files
#'
#' A character string indicating the complete path for storing parsed
#' files.
#' @seealso \code{\link{parseSourceFile}} \code{\linkS4class{ArrayViews}}
#' @export
#' @param object a \code{ArrayViews} object
#' @seealso \code{\linkS4class{ArrayViews}}
setGeneric("parsedPath", function(object) standardGeneric("parsedPath"))

#' Accessor for file paths containing SNP-level summaries
#'
#'  Files containing SNP-level summaries for log R ratios, B allele
#' frequencies, and genotypes -- one sample per subject -- are
#' required.
#'
#' @examples
#' sourcePaths(ArrayViews())
#' @param object an \code{ArrayViews} object
#' @export
#' @aliases sourcePaths,ArrayViews-method
setGeneric("sourcePaths", function(object) standardGeneric("sourcePaths"))

#' Accessors for objects of class ArrayViews
#'
#' @param object see showMethods("lrrFile")
#' @rdname IO
#' @export
setGeneric("lrrFile", function(object) standardGeneric("lrrFile"))

#' @rdname IO
#' @export
setGeneric("lrrFile<-", function(object, value) standardGeneric("lrrFile<-"))


#' @rdname IO
#' @aliases bafFile,ArrayViews-method
#' @export
setGeneric("bafFile", function(object) standardGeneric("bafFile"))

#' @rdname IO
#' @export
#' @aliases gtFile,ArrayViews-method
setGeneric("gtFile", function(object) standardGeneric("gtFile"))

#' Function for parsing GenomeStudio files
#'
#' This function parses genome studio files, writing the low-level
#' data for log R ratios, B allele frequencies, and genotypes to disk
#' as integers (1 file per subject per data type).
#' @param object An \code{\linkS4class{ArrayViews}} object
#' @param param  An object of class \code{\link{CopyNumScanParams}}
#' @seealso \code{\link{ArrayViews}} \code{\link{ArrayViews}} \code{\link{CopyNumScanParams}}
#' @return  NULL
#' @export
setGeneric("parseSourceFile", function(object, param) standardGeneric("parseSourceFile"))

setGeneric("SnpDataFrame",
           function(x, row.names=NULL, check.names=TRUE,
                    isSnp)
           standardGeneric("SnpDataFrame"))

setGeneric("doPosterior", function(object) standardGeneric("doPosterior"))


setGeneric("posterior", function(object) standardGeneric("posterior"))

#' Accessor for HMM model parameters
#'
#' @param object see \code{showMethods(HmmParam)}
#' @aliases getHmmParams,HmmParam-method getHmmParams,HMM-method
#' @export
#' @examples
#' hmm_object <- HMM()
#' getHmmParams(hmm_object)
setGeneric("getHmmParams", function(object) standardGeneric("getHmmParams"))

setGeneric("granges<-", function(x, value) standardGeneric("granges<-"))

#' Filter the HMM-derived genomic ranges for copy number variants
#'
#' The HMM-derived genomic ranges are represented as a
#' \code{GRanges}-derived object.  \code{cnvFilter} returns a
#' \code{GRanges} object using the filters stipulated in the
#' \code{filters} argument.
#' @examples
#' data(snp_exp)
#' fit <- hmm2(snp_exp)
#' segs(fit) ## all intervals
#' cnvSegs(fit)
#' filter_param <- FilterParam(probability=0.95, numberFeatures=10, state=c("1", "2"))
#' cnvSegs(fit, filter_param)
#' filter_param <- FilterParam(probability=0.5, numberFeatures=2, state=c("1", "2"))
#' cnvSegs(fit, filter_param)
#' hemizygous(fit)
#' homozygous(fit)
#' duplication(fit)
#' @param object see \code{showMethods(cnvFilter)}
#' @param filters a \code{\link{FilterParam}} object
#' @seealso \code{\link{FilterParam}}
#' @rdname cnvFilter
#' @export
#' @aliases cnvFilter,HMM-method cnvFilter,GRanges-method
setGeneric("cnvFilter", function(object, filters=FilterParam()) standardGeneric("cnvFilter"))

#' Accessor for the HMM segments
#'
#' Accessor to obtain all segments from the HMM.
#'
#' @return a \code{GRanges}-derived object
#' @aliases segs segs,HMM-method
#' @param object see \code{showMethods(segs)}
#' @export
setGeneric("segs", function(object) standardGeneric("segs"))


#' @rdname cnvFilter
#' @export
setGeneric("cnvSegs", function(object, filters=FilterParam(state=c("1", "2", "5", "6"))) standardGeneric("cnvSegs"))

#' @export
#' @rdname cnvFilter
#' @aliases duplication,HMM-method
setGeneric("duplication", function(object, filters=FilterParam(state=c("5","6"))) standardGeneric("duplication"))

#' @export
#' @rdname cnvFilter
#' @aliases deletion,HMM-method
setGeneric("deletion", function(object, filters=FilterParam(state=c("1","2"))) standardGeneric("deletion"))

#' @aliases hemizygous,HMM-method
#' @export
#' @rdname cnvFilter
setGeneric("hemizygous", function(object, filters=FilterParam(state="2")) standardGeneric("hemizygous"))

#' @export
#' @aliases homozygous,HMM-method
#' @rdname cnvFilter
setGeneric("homozygous", function(object, filters=FilterParam(state="1")) standardGeneric("homozygous"))

#' The number of SNP/nonpolymorphic probes contained in a genomic interval
#'
#' @param object see \code{showMethods(numberFeatures)}
#' @aliases numberFeatures,FilterParam-method numberFeatures,HMM-method numberFeatures,HmmGRanges-method
#' @export
setGeneric("numberFeatures", function(object) standardGeneric("numberFeatures"))

#' Accessor for HMM filter parameters
#'
#' @param object see \code{showMethods(filters)}
#' @aliases filters,HmmParam-method filters,HMM-method
#' @export
setGeneric("filters", function(object) standardGeneric("filters"))

#' Accessor for probability filter
#'
#' @param object a \code{FilterParam} object
#' @export
setGeneric("probability", function(object) standardGeneric("probability"))

setGeneric("state")
setGeneric("lrr")
setGeneric("baf")
setGeneric("copyNumber")


#' Lattice-style plots for granges and SnpArrayExperiment objects
#'
#' Data for the graphic is generated by a call to \code{grangesData}.
#' @param granges a \code{HmmGRanges} object
#' @param se a \code{SnpArrayExperiment}
#' @param param trellis parameters for plotting HMM
#' @rdname plotting
#' @examples
#' if(require("BSgenome.Hsapiens.UCSC.hg18")){
#'   bsgenome <- BSgenome.Hsapiens.UCSC.hg18
#'   snp_exp <- getExampleSnpExperiment(bsgenome)
#'   seqlevels(snp_exp, pruning.mode="coarse") <- "chr22"
#'   fit <- hmm2(snp_exp)
#'   g <- reduce(hemizygous(fit), min.gapwidth=500e3)
#'   trellis_param <- HmmTrellisParam()
#'   fig <- xyplotList(g, snp_exp, trellis_param)
#'   vps <- viewports()
#'   xygrid(fig[[1]], vps, g)
#' }
#' @rdname plotting
#' @export
setGeneric("xyplotList", function(granges, se, param=HmmTrellisParam()) standardGeneric("xyplotList"))



setGeneric("isAutosome", function(object) standardGeneric("isAutosome"))

#' Assess whether genotype is heterozygous based on BAFs
#'
#' @param object a SnpArrayExperiment or ArrayViews object containing
#' BAFs, a matrix of BAFs, or a numeric vector of BAFs.
#' vector of BAFs
#' @param cutoff a length-two numeric vector providing the range of
#' BAFs consistent with allelic  heterozygosity
#' @examples
#' if(require("BSgenome.Hsapiens.UCSC.hg18")){
#'   bsgenome <- BSgenome.Hsapiens.UCSC.hg18
#'   snp_exp <- getExampleSnpExperiment(bsgenome)
#'   is_het <- isHeterozygous(snp_exp[, 1], c(0.4, 0.6))
#'   table(is_het)
#' }
#' @rdname isHeterozygous
#' @export
setGeneric("isHeterozygous", function(object, cutoff) standardGeneric("isHeterozygous"))
