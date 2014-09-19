setMethod(Viterbi, signature(state="missing"),
          function(state=integer(), loglik=numeric(),
                   forward_backward=new("matrix")){
            new("Viterbi", state=state, loglik=loglik,
                forward_backward=forward_backward)
          })

setMethod(Viterbi, signature(state="integer"),
          function(state, loglik=numeric(1),
                   forward_backward=matrix(0.0, length(state),6)){
            new("Viterbi", ##state=Rle(state),
                state=state,
                loglik=loglik,
                forward_backward=forward_backward)
          })

setMethod(Viterbi, signature(state="Rle"),
          function(state, loglik=numeric(1),
                   forward_backward=matrix(0.0, length(state),6)){
            new("Viterbi", state=as.integer(state), loglik=loglik,
                forward_backward=forward_backward)
          })

setMethod(forward_backward, "Viterbi", function(x) x@forward_backward)

#' Show method for objects of class \code{Viterbi}
#' @param object a \code{Viterbi} object
#' @aliases show,Viterbi-method
#' @rdname Viterbi-methods
setMethod(show, signature(object="Viterbi"),
          function(object){
            cat("class Viterbi\n")
            cat("state: vector of length", length(state(object)), "\n")
            cat("loglik:", round(loglik(object), 2), "\n")
            fv <- forward_backward(object)
            cat("forward_backward: ", nrow(fv), "x", ncol(fv), "matrix\n")
          })

## # Accessor for the Viterbi state path
## #

## #
## @param object an object of class \code{Viterbi}
## @return an integer vector



#' Accessor for the Viterbi state path
#'
#' The states are represented as integers: 1=homozygous deletion,
#' 2=hemizygous deletion, 3=diploid normal heterozygosity, 4=diploid
#' region of homozygosity, 5=single copy gain, 6=two or more copy gain.
#' @param object a \code{Viterbi} object
#' @aliases state,Viterbi-method
#' @rdname state-methods
#' @name state-methods
#' @export
setMethod(state, "Viterbi", function(object) object@state)


setValidity("Viterbi", function(object){
  msg <- TRUE
  ll <- loglik(object)
  if(length(ll) == 1){
    if(is.nan(ll)) msg <- "Invalid log lik.  Check that starting values are correctly specified, and that assay data is on the correct scale."
  }
  msg
})


## Viterbi algorithm
##
## Compute the most likely state-path via the Viterbi algorithm. This
## algorithm is exported for internal use by other BioC packages.
##
## @param param an object of class \code{HmmParam}
## @return an objec tof class \code{Viterbi}
## @export
calculateViterbi <- function(param=HmmParam()){
  NR <- nrow(param)
  NC <- ncol(param)
  emitv <- as.numeric(emission(param))
  ##nNA <- sum(is.na(emitv))/NC
  ##NR <- NR_total-nNA
  forvar <- matrix(0.0, NR*NC, 1)
  backvar <- matrix(0.0, NR*NC, 1)
  states <- integer(NR)
  scale_amount <- numeric(NR)
  emit <- as.matrix(emitv)
  res <- .C("viterbi2",
            emit,
            initial(param),
            transition(param),
            integer(NR),
            NC,
            NR,
            states,
            forvar,
            backvar,
            3L,
            scale_amount)
  forvar <- matrix(res[[8]], NR, NC)
  bakvar <- matrix(res[[9]], NR, NC)
  loglik <- -sum(log(res[[11]]))
  FV <- scaleForwardBackward(forvar, bakvar)
  colnames(FV) <- FV_columns()
  Viterbi(state=res[[7]],
          loglik=loglik,
          forward_backward=FV)
}

chromosomeRunLengths <- function(x){
  stateL <- split(x, names(x))
  tmp <- lapply(stateL, Rle)
  runval <- unlist(lapply(tmp, runValue))
  runlen <- unlist(lapply(tmp, runLength))
  Rle(runval, runlen)
}

setMethod(loglik, "Viterbi", function(object) object@loglik)
setMethod(statei, "Viterbi", function(object) as.integer(state(object)))
setMethod(statef, "Viterbi", function(object) factor(as.integer(state(object)), levels=c(1,2,3,4,5,6)))

scaleForwardBackward <- function(forward_var, backward_var){
  forward_backward <- forward_var * backward_var
  totals <- rowSums(forward_backward, na.rm=TRUE)
  forward_backward <- forward_backward/totals
  forward_backward
}
