
setMethod("loglik", "LogLik", function(object) object@loglik)
setMethod("loglikLast", "LogLik", function(object) object@loglik[length(object)])
setMethod("tolerance", "LogLik", function(object) object@tolerance)

setReplaceMethod("loglik", c("LogLik", "numeric"), function(object, value) {
  object@loglik <- c(loglik(object), value)
  object
})

#' @param x object of class \code{LogLik}
#' @export
#' @aliases length,LogLik-method
#' @rdname LogLik-class
setMethod("length", "LogLik", function(x) length(loglik(x)))

setMethod("loglikRatio", "LogLik", function(object) {
  N <- length(object)
  loglikLast(object)-loglik(object)[N-1]
})

#' @export
#' @param object a \code{LogLik} object
#' @rdname LogLik-class
#' @aliases show,LogLik-method
setMethod("show", "LogLik", function(object){
  cat("'LogLik' class\n")
  logliks <- paste(tail(round(loglik(object), 2)), collapse=", ")
  cat("  iterations:", length(object), "\n")
  cat("  log lik:", logliks, "\n")
  llr <- round(loglikRatio(object), 2)
  if(length(llr) > 0)
    if(llr < 0) message("LLR less than 0 do to constraints.")
  cat("  log lik ratio:", llr, "\n")
})

#' Constructor for LogLik class
#'
#' A container for the log likelihood of the Viterbi state path.
#' Stores the log likelihood from succesive updates of model
#' parameters.  When the difference between the log likelihoods at
#' iteration i and i-1 is below the tolerance, no additional updates
#' are performed.
#'
#' @param loglik length-one numeric vector for the log likelihood of the Viterbi state path
#' @param tolerance if the difference in the log-likelihood of the Viterbi state path after the Baum-Welch update is less than the specified tolerance, no additional Baum-Welch updates are required
#' @rdname LogLik
#' @seealso \linkS4class{LogLik}
#' @export
LogLik <- function(loglik=numeric(), tolerance=1L){
  new("LogLik", loglik=loglik, tolerance=tolerance)
}

setMethod("exceedsTolerance", "LogLik", function(object){
  if(length(object) == 1) return(TRUE)
  loglikRatio(object) > tolerance(object)
})
