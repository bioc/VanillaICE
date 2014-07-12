#' @export
setGeneric("hmm2", function(object, emission_param=EmissionParam(),
                            transition_param=TransitionParam(),
                            tolerance=2,
                            verbose=FALSE, ...) standardGeneric("hmm2"))

#' @export
setGeneric("loglik", function(object) standardGeneric("loglik"))

#' @export
setGeneric("calculateEmission", function(x, param=EmissionParam()) standardGeneric("calculateEmission"))

setGeneric("forward_backward", function(x) standardGeneric("forward_backward"))

#' @export
setGeneric("Viterbi",
           function(state=Rle(), loglik=numeric(1),
                    forward_backward=matrix(0.0, length(state), 6))
           standardGeneric("Viterbi"))

#' @export
setGeneric("transition", function(object) standardGeneric("transition"))

#' @export
setGeneric("initial", function(object) standardGeneric("initial"))

#' @export
setGeneric("statei", function(object) standardGeneric("statei"))

setGeneric("state<-", function(object,value) standardGeneric("state<-"))

#' @export
setGeneric("statef", function(object) standardGeneric("statef"))

#' @export
setGeneric("cn_means", function(object) standardGeneric("cn_means"))

#' @export
setGeneric("cn_sds", function(object) standardGeneric("cn_sds"))

#' @export
setGeneric("baf_means", function(object) standardGeneric("baf_means"))

#' @export
setGeneric("baf_sds", function(object) standardGeneric("baf_sds"))

#' @export
setGeneric("baf_means<-", function(object, value) standardGeneric("baf_means<-"))

#' @export
setGeneric("baf_sds<-", function(object, value) standardGeneric("baf_sds<-"))

#' @export
setGeneric("cn_sds<-", function(object, value) standardGeneric("cn_sds<-"))

#' @export
setGeneric("cn_means<-", function(object, value) standardGeneric("cn_means<-"))

#' @export
setGeneric("taup", function(object) standardGeneric("taup"))

#' @export
setGeneric("taumax", function(object) standardGeneric("taumax"))

#' @export
setGeneric("EmissionParam", function(cn_means=CN_MEANS(),
                                     cn_sds=CN_SDS(),
                                     baf_means=BAF_MEANS(),
                                     baf_sds=BAF_SDS(),
                                     initial=rep(1/6, 6),
                                     EMupdates=5L,
                                     CN_range=c(-5, 3),
                                     proportionOutlier=0.01,
                                     temper=25)
           standardGeneric("EmissionParam"))

#' @export
setGeneric("EMupdates", function(object) standardGeneric("EMupdates"))

#' @export
setGeneric("EMupdates<-", function(object,value) standardGeneric("EMupdates<-"))

#' @export
setGeneric("emission", function(object) standardGeneric("emission"))

#' @export
setGeneric("emission<-", function(object, value) standardGeneric("emission<-"))

#' @export
setGeneric("TransitionParam",
           function(taup=1e10,
                    taumax=1-5e6)
           standardGeneric("TransitionParam"))

#' @export
setGeneric("HmmParam",
           function(emission=matrix(0.0, 0, 0),
                    ##initial=rep(1/6, 6),
                    emission_param=EmissionParam(),
                    transition=rep(0.99, nrow(emission)),
                    chromosome=character(nrow(emission)),
                    loglik=LogLik(),
                    viterbi=Viterbi(),
                    verbose=FALSE)
           standardGeneric("HmmParam"))

#' @export
setGeneric("calculateTransitionProbability",
           function(x, param=TransitionParam())
           standardGeneric("calculateTransitionProbability"))

#' @export
setGeneric("modev", function(x) standardGeneric("modev"))



#' @export
setGeneric("SnpGRanges", function(object=GRanges(),
                                  isSnp, ...)
           standardGeneric("SnpGRanges"))

#' @export
setGeneric("BeadStudioSetList", function(x, ...)
           standardGeneric("BeadStudioSetList"))

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
#' @param rowData GRanges object for SNPs/nonpolymorphic markers
#' @param colData DataFrame containing sample-level covariates
#' @param isSnp  logical vector indicating whether marker is a SNP
#' @aliases SnpArrayExperiment
#' @docType methods
#' @rdname SnpArrayExperiment-class
#' @export
setGeneric("SnpArrayExperiment", function(cn,
                                          baf,
                                          rowData=GRanges(),
                                          colData=DataFrame(),
                                          isSnp=logical(), ...)
           standardGeneric("SnpArrayExperiment"))

#' @export
setGeneric("sweepMode", function(x, MARGIN) standardGeneric("sweepMode"))

#' @export
setGeneric("NA_index", function(x) standardGeneric("NA_index"))

#' HmmGRanges container
#
#' @param states  copy number number state inferred by HMM
#' @param feature_starts start location in reference genome [basepairs]
#' @param feature_chrom  end location in reference genome [basepairs]
#' @param loglik  the log likelihood
#' @aliases HmmGRanges
#' @rdname HmmGRanges-class
#' @export
setGeneric("HmmGRanges", function(states, feature_starts,
                                  feature_chrom, loglik)
           standardGeneric("HmmGRanges"))

setGeneric("isSnp", function(object) standardGeneric("isSnp"))

#' @export
setGeneric("loglikRatio", function(object) standardGeneric("loglikRatio"))

#' @export
setGeneric("tolerance", function(object) standardGeneric("tolerance"))

#' @export
setGeneric("loglik", function(object) standardGeneric("loglik"))
setGeneric("loglikLast", function(object) standardGeneric("loglikLast"))

#' @export
setGeneric("loglik<-", function(object,value) standardGeneric("loglik<-"))

#' @export
setGeneric("exceedsTolerance", function(object) standardGeneric("exceedsTolerance"))

#' @export
setGeneric("emissionParam", function(object) standardGeneric("emissionParam"))

#' @export
setGeneric("emissionParam<-", function(object,value) standardGeneric("emissionParam<-"))

#' @export
setGeneric("viterbi", function(object) standardGeneric("viterbi"))
#' @export
setGeneric("viterbi<-", function(object,value) standardGeneric("viterbi<-"))
#' @export
setGeneric("verbose", function(object) standardGeneric("verbose"))

setGeneric("CN_range", function(object) standardGeneric("CN_range"))
setGeneric("proportionOutlier", function(object) standardGeneric("proportionOutlier"))
setGeneric("temper", function(object) standardGeneric("temper"))

##setGeneric("scale")
setGeneric("gstudioPaths", function(object) standardGeneric("gstudioPaths"))
setGeneric("indexGenome", function(object) standardGeneric("indexGenome"))
setGeneric("read", function(object) standardGeneric("read"))
setGeneric("genotypes", function(object) standardGeneric("genotypes"))
setGeneric("SnpExperiment", function(object) standardGeneric("SnpExperiment"))
setGeneric("scaleBy", function(x, by) standardGeneric("scaleBy"))
setGeneric("scaleRead", function(x, params) standardGeneric("scaleRead"))
setGeneric("selectCols", function(object) standardGeneric("selectCols"))
setGeneric("datadir", function(object) standardGeneric("datadir"))
#' @export
setGeneric("fileName", function(object, label) standardGeneric("fileName"))
#' @export
setGeneric("filePaths", function(object) standardGeneric("filePaths"))

#' @export
setGeneric("lrrFile", function(object, label="lrr") standardGeneric("lrrFile"))
#' @export
setGeneric("bafFile", function(object, label="baf") standardGeneric("bafFile"))

#' @export
setGeneric("gtFile", function(object, label="gt") standardGeneric("gtFile"))

setGeneric("parseGStudio", function(object, param) standardGeneric("parseGStudio"))

#' @export
setGeneric("SnpDataFrame",
           function(x, row.names=NULL, check.names=TRUE,
                    isSnp)
           standardGeneric("SnpDataFrame"))

##setGeneric("SummarizedHMM", function(forward_backward,
##                                      loglik,
##                                      rowData=HmmGRanges())
##           standardGeneric("SummarizedHMM"))
