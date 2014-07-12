## EmissionParam: needs to know means and sds.  Number of states
## inferred by length of these vectors
##
## TransitionParam:  needs to know 'taup' and 'taumax'
##
## HmmParam: needs initial state, transition,
## and emission
##

setMethod(HmmParam, signature(emission="missing"),
          function(emission=matrix(0.0, 0, 0),
                   emission_param=EmissionParam(),
                   transition=rep(1-exp(-10), nrow(emission)),
                   chromosome=character(nrow(emission)),
                   loglik=LogLik(),
                   viterbi=Viterbi(),
                   verbose=FALSE){
            new("HmmParam",
                emission=emission,
                emission_param=emission_param,
                transition=transition,
                chromosome=chromosome,
                loglik=loglik,
                viterbi=viterbi,
                verbose=verbose)
          })

setMethod(HmmParam, signature(emission="matrix"),
          function(emission,
                   emission_param=EmissionParam(),
                   transition=rep(1-exp(-10), nrow(emission)),
                   chromosome=character(nrow(emission)),
                   loglik=LogLik(),
                   viterbi=Viterbi(),
                   verbose=FALSE){
            new("HmmParam",
                emission=emission,
                emission_param=emission_param,
                transition=transition,
                chromosome=chromosome,
                loglik=loglik,
                viterbi=viterbi,
                verbose=verbose)
          })



setMethod("emissionParam", "HmmParam", function(object) object@emission_param)

setMethod("verbose", "HmmParam", function(object) object@verbose)

setMethod("viterbi", "HmmParam", function(object) object@viterbi)
setReplaceMethod("viterbi", c("HmmParam", "Viterbi"), function(object, value){
  object@viterbi <- value
  object
})

setMethod("EMupdates", "HmmParam", function(object) EMupdates(emissionParam(object)))

setReplaceMethod("emissionParam", c("HmmParam", "EmissionParam"), function(object, value){
  object@emission_param <- value
  object
})

setMethod(show, signature(object="HmmParam"),
          function(object){
            cat(class(object), ":\n")
            emissions <- emission(object)
            cat("emission:", nrow(emissions), "rows,", ncol(emissions), "states \n")
            transitions <- transition(object)
            cat("transition:  min=", min(transitions), ", max=", max(transitions), "\n")
            init <- initial(object)
            cat("initial state probabilities: ", paste(round(init, 2), collapse=", "), "\n")
            cat("Log likelihood: \n")
            show(loglik(object))
          })

setMethod("loglik", "HmmParam", function(object) object@loglik)

updateLogLik <- function(object) loglik(viterbi(object))

setReplaceMethod("loglik", c("HmmParam", "LogLik"), function(object, value){
  object@loglik <- value
  object
})

setReplaceMethod("loglik", c("HmmParam", "numeric"), function(object, value){
  LL <- loglik(object)
  loglik(LL) <- value
  object@loglik <- LL
  object
})

setMethod("exceedsTolerance", "HmmParam", function(object){
  exceedsTolerance(loglik(object))
})

setMethod("emission", signature(object="HmmParam"), function(object) object@emission)

setMethod("chromosome", signature(object="HmmParam"), function(object) object@chromosome)

#' @export
setReplaceMethod("emission", signature(object="HmmParam"),
                 function(object,value){
                   object@emission <- value
                   object
                 })

setMethod("nrow", signature(x="HmmParam"), function(x){
  nrow(emission(x))
})

setMethod("ncol", signature(x="HmmParam"), function(x){
  ncol(emission(x))
})

setMethod("transition", signature(object="HmmParam"), function(object) {
  object@transition
})

setMethod("initial", signature(object="HmmParam"), function(object) {
  initial(emissionParam(object))
})

setValidity("HmmParam", function(object){
  emit <- emission(object)
  if(any(emit < 0, na.rm=TRUE)){
    msg <- "emissions must be non-negative"
  } else msg <- TRUE
  msg
})

setMethod("cn_means", "HmmParam", function(object) cn_means(emissionParam(object)))
setMethod("cn_sds", "HmmParam", function(object) cn_sds(emissionParam(object)))
setMethod("baf_means", "HmmParam", function(object) baf_means(emissionParam(object)))
setMethod("baf_sds", "HmmParam", function(object) baf_sds(emissionParam(object)))
