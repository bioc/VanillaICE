setMethod(show, "EmissionParam",
          function(object) {
            cat(class(object), ":\n")
            means <- round(cn_means(object), 2)
            cat("CN mean: ", paste(means, collapse=", "), "\n")
            sds <- paste(round(cn_sds(object), 2), collapse=", ")
            cat("CN sd: ", sds, "\n")
            cat("BAF mean: ", paste(round(baf_means(object), 2), collapse=", "), "\n")
            cat("BAF sd: ", paste(round(baf_sds(object), 2), collapse=", "), "\n")
            cat("max # of EM updates: ", EMupdates(object), "\n")
            init <- initial(object)
            cat("initial state probabilities: ", paste(round(init, 2), collapse=", "), "\n")
            cat("probability outlier:", proportionOutlier(object), "\n")
            cat("tempering scalar: ", temper(object), "\n")
            cat("See temper(), proportionOutlier(), cn_means(), cn_sds(),...\n")
          })

setMethod("initial", "EmissionParam", function(object) object@initial)

setMethod(EmissionParam, signature(cn_means="missing"),
          function(cn_means=CN_MEANS(),
                   cn_sds=CN_SDS(),
                   baf_means=BAF_MEANS(),
                   baf_sds=BAF_SDS(),
                   initial=rep(1/6, 6),
                   EMupdates=5L,
                   CN_range=c(-5, 3),
                   proportionOutlier=0.01,
                   temper=25){ ## flatten the peaks
            b_means <- setNames(baf_means, BAF_ALLELE_NAMES())
            b_sds <- setNames(baf_sds, BAF_ALLELE_NAMES())
            new("EmissionParam", cn_means=cn_means, cn_sds=cn_sds,
                baf_means=b_means, baf_sds=b_sds,
                initial=initial,
                EMupdates=EMupdates,
                CN_range=CN_range,
                proportionOutlier=proportionOutlier, temper=temper)
          })

setMethod(EmissionParam, signature(cn_means="numeric"),
          function(cn_means=cn_means,
                   cn_sds=CN_SDS(),
                   baf_means=BAF_MEANS(),
                   baf_sds=BAF_SDS(),
                   initial=rep(1/6, 6),
                   EMupdates=5L,
                   CN_range=c(-5,3),
                   proportionOutlier=0.01,
                   temper=25){
            b_means <- setNames(baf_means, BAF_ALLELE_NAMES())
            b_sds <- setNames(baf_sds, BAF_ALLELE_NAMES())
            new("EmissionParam", cn_means=cn_means, cn_sds=cn_sds,
                baf_means=b_means, baf_sds=b_sds,
                initial=initial,
                EMupdates=EMupdates,
                CN_range=CN_range,
                proportionOutlier=proportionOutlier,
                temper=temper)
          })

setValidity("EmissionParam", function(object){
  msg <- NULL
  if(length(cn_means(object)) != length(cn_sds(object)))
    msg <- c(msg, "Means and standard deviation vectors for copy number must be the same length\n")
  if(length(baf_means(object)) != length(baf_sds(object)))
    msg <- c(msg, "Means and standard deviation vectors for BAFs must be the same length\n")
  if(length(cn_means(object)) != 6)
    msg <- c(msg, "Copy number means must be a numeric vector of length 6\n")
  if(length(cn_sds(object)) != 6)
    msg <- c(msg, "Copy number sds must be a numeric vector of length 6\n")
  if(length(baf_means(object)) != 7)
    msg <- c(msg, "BAF means must be a numeric vector of length 7\n")
  if(length(baf_sds(object)) != 7)
    msg <- c(msg, "BAF sds must be a numeric vector of length 7\n")
})


setMethod(cn_means, "EmissionParam", function(object) object@cn_means)
setMethod(cn_sds, "EmissionParam", function(object) object@cn_sds)
setMethod(baf_means, "EmissionParam", function(object) object@baf_means)
setMethod(baf_sds, "EmissionParam", function(object) object@baf_sds)
setMethod(taup, "TransitionParam", function(object) object@taup)
setMethod(taumax, "TransitionParam", function(object) object@taumax)
setMethod(EMupdates, "EmissionParam", function(object) object@EMupdates)
setReplaceMethod("EMupdates", c("EmissionParam", "integer"), function(object, value) {
  object@EMupdates <- value
  object
})

setReplaceMethod("baf_sds", c("EmissionParam", "numeric"),
                 function(object, value){
                   object@baf_sds <- value
                   object
                 })

setReplaceMethod("cn_sds", c("EmissionParam", "numeric"),
                 function(object, value){
                   object@cn_sds <- value
                   object
                 })

setReplaceMethod("cn_means", c("EmissionParam", "numeric"),
                 function(object, value){
                   object@cn_means <- value
                   object
                 })

setReplaceMethod("baf_means", c("EmissionParam", "numeric"),
                 function(object, value){
                   object@baf_means <- value
                   object
                 })


setMethod(calculateEmission, signature(x="numeric"),
          function(x, param=EmissionParam()){
            means <- cn_means(param)
            sds <- cn_sds(param)
            limits <- CN_range(param)
            x <- threshold(x, limits)
            emit_cn <- mapply(dnorm, mean=as.list(means), sd=as.list(sds),
                              MoreArgs=list(x=x))
            p <- proportionOutlier(param)
            emit <- (1-p)*emit_cn + p*dunif(0, limits[1], limits[2])
            emit
          })

.hmm_states <- function() c("cn0", "cn1", "cn2", "cn2-loh", "cn3", "cn4")

.calculateBAFEmission <- function(x, param=EmissionParam()){
  means <- baf_means(param)
  sds <- baf_sds(param)
  ##
  ## copy number 1 (genotype A or B)
  ##
  pr <- mapply(dtrnorm, mean=as.list(means),
               sd=as.list(sds),
               MoreArgs=list(x=x))
  ##
  ## Temper the emission probabilities
  ##
  cn1 <- prC1(pr)
  cn2 <- prC2(pr)
  cn3 <- prC3(pr)
  cn4 <- prC4(pr)
  ## uniform for homozygous deletions
  x <- cbind(1, cn1, cn2, cn1, cn3, cn4)
  colnames(x) <- .hmm_states()
  x
}

EmissionView <- function(x, emit){
  colnames(emit) <- HMM_STATES()
  cbind(round(x,2), round(emit,2))
}


setMethod("CN_range", "EmissionParam", function(object) object@CN_range)
setMethod("proportionOutlier", "EmissionParam", function(object) object@proportionOutlier)

setMethod(calculateEmission, signature(x="list"),
          function(x, param=EmissionParam()){
            .calculateEmission(x, param)
          })

.calculateEmission <- function(x, param){
  means <- cn_means(param)
  sds <- cn_sds(param)
  limits <- CN_range(param)
  x[[1]] <- threshold(x[[1]], limits)
  emit_cn <- mapply(dnorm, mean=as.list(means), sd=as.list(sds),
                    MoreArgs=list(x=x[[1]]))
  emit_baf <- .calculateBAFEmission(x[[2]], param)
  ## must stay on probability scale for .viterbi2 (do not log)
  ##
  ## Temper the CN emission probabilities because of outliers
  scalar <- temper(param)
  p <- 1/scalar
  emit_cn <- (p)*emit_cn + (1-p)*dunif(0, limits[1], limits[2])
  ##emit_baf <- (1-1e-5)*emit_baf + 1e-5*1
  emit <- emit_cn * emit_baf
  colnames(emit) <- .hmm_states()
  return(emit)
}

setMethod(calculateEmission, signature(x="SummarizedExperiment"),
          function(x, param=EmissionParam()){
            x <- list(lrr(x)[,1], baf(x)[,1])
            .calculateEmission(x, param)
          })


setMethod("temper", "EmissionParam", function(object) object@temper)
