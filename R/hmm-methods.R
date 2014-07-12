assayList <- function(x) list(copyNumber(x)[, 1], baf(x)[, 1])

#' @export
updateHmmParams <- function(object, emission_param=EmissionParam(), transition_param=TransitionParam(),...){
  object <- NA_filter(object)
  assay_list <- assayList(object)
  chrom <- chromosome(object)
  emissions <- calculateEmission(assay_list, emission_param)
  transition_prob <- calculateTransitionProbability(object, transition_param)
  hmm_param <- HmmParam(emission=emissions,
                        emission_param=emission_param,
                        transition=transition_prob,
                        chrom=chrom)
  while(doUpdate(hmm_param)){
    hmm_param <- baumWelchUpdate(hmm_param, assay_list)
  }
  hmm_param
}

.hmm <- function(object,
                 emission_param=EmissionParam(),
                 transition_param,
                 ...){
  hmm_param <- updateHmmParams(object, emission_param=emission_param, transition_param=transition_param, ...)
  fit <- viterbi(hmm_param)
  hgr <- HmmGRanges(state(fit),
                    feature_starts=start(object),
                    feature_chrom=chromosome(object),
                    loglik=loglik(fit))
  hgr
}

#' @export
doUpdate <- function(param) {
  if(length(loglik(param)) == 0) return(TRUE)
  exceedsTolerance(param) && length(loglik(param)) < EMupdates(param)
}

#' @export
baumWelchUpdate <- function(param, assay_list){
  fit <- calculateViterbi(param)
  LL <- loglik(param)
  loglik(LL) <- loglik(fit)
  e_param <- updateParam(assay_list, emissionParam(param), fit)
  emit <- calculateEmission(assay_list, e_param)
  HmmParam(emission=emit,
           emission_param=e_param,
           transition=transition(param),
           chromosome=chromosome(param),
           loglik=LL,
           viterbi=fit,
           verbose=verbose(param))
}

setMethod(hmm2, "SnpArrayExperiment",
          function(object,
                   emission_param=EmissionParam(),
                   transition_param=TransitionParam(),
                   tolerance=2,
                   verbose=FALSE,
                   ...){
            ##transition_prob <- calculateTransitionProbability(object, transition_param)
            results <- foreach(j = seq_len(ncol(object))) %dopar% {
              obj <- object[, j]
              .hmm(obj,
                   emission_param=emission_param,
                   transition_param=transition_param,
                   tolerance=tolerance,
                   verbose=verbose, ...)
            }
            results <- setNames(GRangesList(results), colnames(object))
            results
          })


#' Retrieve genomic location of SNPs
#'
#' @export
setMethod("start", "oligoSnpSet", function(x) featureData(x)$position)

setMethod(hmm2, "oligoSnpSet",
          function(object, emission_param=EmissionParam(),
                   transition_param=TransitionParam(),
                   tolerance=2,
                   verbose=FALSE,
                   ...){
            object <- NA_filter(object)
            ## needs to be calculated for the object
            ##transition_prob <- calculateTransitionProbability(object, transition_param)
            .hmm(object,
                 emission_param=emission_param,
                 transition_param=transition_param,
                 tolerance=tolerance,
                 verbose=verbose, ...)
          })
