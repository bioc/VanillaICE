assayList <- function(x) list(copyNumber(x)[, 1], baf(x)[, 1])

#' Run the Baum-Welch algorithm to update HMM parameters
#'
#' This function is not intended to be called directly by the user. It
#' is exported in the package NAMESPACE for internal use by other BioC
#' packages.
#'
#' @param object a \code{\link{SnpArrayExperiment}} object
#' @param emission_param a \code{\link{EmissionParam}} object
#' @param transition_param a \code{\link{TransitionParam}} object
#' @export
updateHmmParams <- function(object, emission_param=EmissionParam(), transition_param=TransitionParam()){
  assay_list <- assayList(object)
  chrom <- chromosome(object)
  emissions <- calculateEmission(assay_list, emission_param)
  transition_prob <- calculateTransitionProbability(object, transition_param)
  hmm_param <- HmmParam(emission=emissions,
                        emission_param=emission_param,
                        transition=transition_prob,
                        chromosome=chrom)
  while(doUpdate(hmm_param)){
    hmm_param <- baumWelchUpdate(hmm_param, assay_list)
  }
  hmm_param
}

.hmm <- function(object,
                 emission_param=EmissionParam(),
                 transition_param,
                 ...){
  object <- NA_filter(object)
  hmm_param <- updateHmmParams(object, emission_param=emission_param, transition_param=transition_param, ...)
  fit <- viterbi(hmm_param)
  hgr <- HmmGRanges(state(fit),
                    feature_starts=start(object),
                    feature_chrom=chromosome(object),
                    loglik=loglikLast(loglik(hmm_param)),
                    emission_param=emissionParam(hmm_param))
  if(doPosterior(hmm_param)){
    probs <- posteriorProb(hmm_param, fit, chromosome(object))
    probs <- round(probs, 4)
    prCall <- mapply(function(state, index, prob) prob[index, state],
                     state=state(hgr),
                     index=seq_len(length(hgr)),
                     MoreArgs=list(prob=probs))
    hgr$prCall <- prCall
  }
  HMM(granges=hgr, param=hmm_param, posterior=probs)
}





posteriorProb <- function(hmm_param, fit, chrom){
  states <- paste0(chrom, ":", state(fit))
  emit <- emission(hmm_param)^(temper(hmm_param))
  lemit <- log(emit)
  ##
  ## compute the column sums of the log(emission) matrix for each
  ## segment
  ##
  ##  -- this is the log likeihood of the data (assuming independence
  ##  -- given the underling state) for each model (state)
  ##
  ## Find the indices for each segment
  ##
  diff_state <- states[-length(states)] != states[-1]
  segs <- cumsum(c(0, diff_state))
  index <- seq_len(nrow(lemit))
  indexlist <- split(index, segs)
  ##
  ## Calculate the column sum for each segment (see above)
  ##
  loglik <- foreach(i = indexlist) %do% colSums(lemit[i, , drop=FALSE])
  loglik <- do.call(rbind, loglik)
  ## For each segment, compute
  ##
  ## Pr(state | data, previous state) = Pr(data | state, previous state) * Pr(state|previous state)
  ##                                  = likihood * tau
  ##
  ## where tau is the transition probability and the likelihood is
  ## calculated as defined above
  ##
  NS <- numberStates(hmm_param)
  transition_index <- which(diff_state)
  previous_state <- as.integer(sapply(strsplit(states[transition_index], ":"), '[', 2))
  current_state <- as.integer(sapply(strsplit(states[transition_index+1], ":"), '[', 2))
  ##current_state <- states[transition_index+1]
  ##stopifnot(all(previous_state != current_state)) ##just to verify indexing
  tau <- transition(hmm_param)[transition_index]
  log_pr_different_state <- log((1-tau)/(NS-1))
  log_pr_same_state <- log(tau)
  ## multiply log likelihood + log transition probability
  probs <- matrix(NA, nrow(loglik), ncol(loglik))
  colnames(probs) <- colnames(loglik)
  ## not sure that this can be vectorized
  for(i in seq_len(length(tau))){
    x <- rep(log_pr_different_state[i], NS)
    x[previous_state[i]] <- log_pr_same_state[i]  ## sum(exp(x))=1
    tmp <- loglik[i+1, ] + x
    ## to avoid overflow
    tmp <- tmp-max(tmp)
    lik.tau <- exp(tmp)
    probs[i+1, ] <- lik.tau/sum(lik.tau)
  }
  probs[1, ] <- .firstSegment(loglik[1, ])
  ##if(any(is.nan(probs))) browser()
  ##probs[!is.finite(probs)] <- 1
  probs
}

.firstSegment <- function(loglik){
  tmp <- loglik + log(1/6)
  tmp <- tmp-max(tmp)
  lik.tau <- exp(tmp)
  lik.tau/sum(lik.tau)
}

#' Helper function to determine whether to update the HMM parameters via the Baum-Welch algorithm
#'
#' This function is not intended to be called directly by the user,
#' and is exported only for internal use by other BioC packages.
#'
#' @param param An object containing parameters for the HMM
#' @seealso \code{\link{HmmParam}}
#' @export
doUpdate <- function(param) {
  if(length(loglik(param)) == 0) return(TRUE)
  exceedsTolerance(param) && length(loglik(param)) < EMupdates(param)
}

#' Function for updating parameters for emission probabilities
#'
#' This function is not meant to be called directly by the user.  It
#' is exported in the package NAMESPACE for internal use by other BioC
#' packages.
#'
#' @param param  A container for the HMM parameters
#' @param assay_list list of log R ratios and B allele frequencies
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

#' @examples
#' data(snp_exp)
#' emission_param <- EmissionParam(temper=1/2)
#' fit <- hmm2(snp_exp, emission_param)
#' unlist(fit)
#' cnvSegs(fit)
#' ## There is too little data to infer cnv reliably in this trivial example.
#' ## To illustrate filtering options on the results, we select
#' ## CNVs for which
#' ## - the CNV call has a posterior probability of at least 0.5
#' ## - the number of features is 2 or more
#' ## - the HMM states are 1 (homozygous deletion) or 2 (hemizygous deletion)
#' fp <- FilterParam(probability=0.5, numberFeatures=2, state=c("1", "2"))
#' cnvSegs(fit, fp)
#' ## for parallelization
#' \dontrun{
#'    library(snow)
#'    library(doSNOW)
#'    cl <- makeCluster(2, type = "SOCK")
#'    registerDoSNOW(cl)
#'    fit <- hmm2(snp_exp, emission_param)
#' }
#' @aliases hmm2,SnpArrayExperiment-method
#' @rdname hmm2
setMethod(hmm2, "SnpArrayExperiment",
          function(object,
                   emission_param=EmissionParam(),
                   transition_param=TransitionParam(),
                   ...){
            ##transition_prob <- calculateTransitionProbability(object, transition_param)
            J <- ncol(object)
            j <- NULL
            results <- foreach(j= seq_len(ncol(object)), .packages="VanillaICE") %dopar% {
              hgr <- .hmm(object[, j],
                          emission_param=emission_param,
                          transition_param=transition_param,  ...)
              ##if(J > 1) emission(hgr) <- NULL
              emission(hgr) <- NULL
              seqinfo(hgr) <- seqinfo(object)
              hgr
            }
            names(results) <- colnames(object)
            ##unlist(HMMList(results))
            results <- HMMList(results)
            results
          })



#' Retrieve genomic location of SNPs
#' @param x a \code{oligoSnpSet} object
#' @aliases start,oligoSnpSet-method
#' @export
setMethod("start", "oligoSnpSet", function(x) featureData(x)$position)



#' @aliases hmm2,oligoSnpSet-method
#' @rdname hmm2
setMethod(hmm2, "oligoSnpSet",
          function(object, emission_param=EmissionParam(),
                   transition_param=TransitionParam(),
                   ...){
            cn <- copyNumber(object)[,1]
            if(is(cn, "integer")) stop("copy number in oligoSnpSet should be numeric and not integer")
            ## needs to be calculated for the object
            ##transition_prob <- calculateTransitionProbability(object, transition_param)
            J <- ncol(object)
            j <- NULL
            results <- foreach(j= seq_len(ncol(object)), .packages="VanillaICE") %dopar% {
              hgr <- .hmm(object[, j],
                          emission_param=emission_param,
                          transition_param=transition_param,  ...)
              if(J > 1) emission(hgr) <- NULL
              ##seqinfo(hgr) <- seqinfo(object)
              hgr
            }
            names(results) <- colnames(object)
##            .hmm(object,
##                 emission_param=emission_param,
            ##                 transition_param=transition_param, ...)
            results <- HMMList(results)
            results
          })
