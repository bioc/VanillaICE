setMethod("loglik", "LogLik", function(object) object@loglik)
setMethod("loglikLast", "LogLik", function(object) object@loglik[length(object)])
setMethod("tolerance", "LogLik", function(object) object@tolerance)

setReplaceMethod("loglik", c("LogLik", "numeric"), function(object, value) {
  object@loglik <- c(loglik(object), value)
  object
})

setMethod("length", "LogLik", function(x) length(loglik(x)))

setMethod("loglikRatio", "LogLik", function(object) {
  N <- length(object)
  loglikLast(object)-loglik(object)[N-1]
})

setMethod("show", "LogLik", function(object){
  cat("'LogLik' class\n")
  logliks <- paste(tail(round(loglik(object), 2)), collapse=", ")
  cat("  iterations:", length(object), "\n")
  cat("  log lik:", logliks, "\n")
  llr <- round(loglikRatio(object), 2)
  cat("  log lik ratio:", llr, "\n")
})

#' @export
LogLik <- function(loglik=numeric(), tolerance=1L){
  new("LogLik", loglik=loglik, tolerance=tolerance)
}

setMethod("exceedsTolerance", "LogLik", function(object){
  if(length(object) == 1) return(TRUE)
  loglikRatio(object) > tolerance(object)
})
