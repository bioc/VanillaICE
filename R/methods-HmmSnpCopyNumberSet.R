setMethod("initialize", "HmmSnpCopyNumberSet",
          function(.Object,
                   assayData = assayDataNew(
                     copyNumber=copyNumber,
                     cnConfidence=cnConfidence,
                     predictions=predictions),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   copyNumber=new("matrix"),
                   cnConfidence=new("matrix"),
                   predictions=new("matrix"),
                   scaleCopyNumber=numeric(),
                   log2=logical(),
                   locationCopyNumber=numeric(),
                   stateNames=character(),
                   initialStateProbability=NULL,
                   transitionProbability=NULL,
                   emissionProbability=NULL,
                   likelihood=NULL,
                   copyNumberIce=logical(),
                   distance=logical(),
                   ...) {
            assayData(.Object) <- assayData
            featureData(.Object) <- featureData
            phenoData(.Object) <- phenoData
            experimentData(.Object) <- experimentData
            annotation(.Object) <- annotation
            scaleCopyNumber(.Object) <- scaleCopyNumber
            locationCopyNumber(.Object) <- locationCopyNumber
            .Object@log2 <- log2
            predictions(.Object) <- predictions
            stateNames(.Object) <- stateNames
            initialStateProbability(.Object) <- initialStateProbability
            transitionProbability(.Object) <- transitionProbability
            emissionProbability(.Object) <- emissionProbability
            likelihood(.Object) <- likelihood
            copyNumberIce(.Object) <- copyNumberIce
            distance(.Object) <- distance
            .Object
          })

setValidity("HmmSnpCopyNumberSet", function(object){
  assayDataValidMembers(assayData(object), c("copyNumber", "cnConfidence", "predictions"))
})

##Copied from Biobase/R/methods-eSet.R
##methods used to subset assayData are not exported from Biobase
setMethod("[", "HmmSnpCopyNumberSet", function(x, i, j, ..., drop = FALSE) {
  if (missing(drop)) drop <- FALSE
  if (missing(i) && missing(j)) {
      if (length(list(...))!=0)
        stop("specify genes or samples to subset; use '",
             substitute(x), "$", names(list(...))[[1]],
             "' to access phenoData variables")
      return(x)
  }
  if(!missing(i)){
    if(!is.null(likelihood(x))){    
      if(nrow(likelihood(x)) == nrow(x))
        likelihood(x) <- likelihood(x)[i, ,drop=FALSE]
    }
  }
  if(!missing(i)){
    if(!is.null(emissionProbability(x))){
      if(nrow(emissionProbability(x)) == nrow(x))    
        emissionProbability(x) <- emissionProbability(x)[i, , drop=FALSE]
    }
  }      
  if (!missing(j))
    phenoData(x) <- phenoData(x)[j,, ..., drop = drop]
  if (!missing(i))
    featureData(x) <- featureData(x)[i,,..., drop=drop]
  ## assayData; implemented here to avoid function call
  if(!missing(i)){
    copyNumber(x) <- copyNumber(x)[i, , drop=FALSE]
    cnConfidence(x) <- cnConfidence(x)[i, , drop=FALSE]
    predictions(x) <- predictions(x)[i, , drop=FALSE]        
  }
  if(!missing(j)){
    copyNumber(x) <- copyNumber(x)[, j, drop=FALSE]
    cnConfidence(x) <- cnConfidence(x)[, j, drop=FALSE]
    predictions(x) <- predictions(x)[, j, drop=FALSE]        
  }  
  x
})

setMethod("bind", c("HmmSnpCopyNumberSet", "HmmSnpCopyNumberSet"),
          function(object1, object2, ...){
            fd <- rbind(fData(object1), fData(object2))
            fD <- new("AnnotatedDataFrame", data=fd,
                      varMetadata=varMetadata(featureData(object1)))
            new("HmmSnpCopyNumberSet",
                copyNumber=rbind(copyNumber(object1), copyNumber(object2)),
                cnConfidence=rbind(cnConfidence(object1), cnConfidence(object2)),
                featureData=fD,
                phenoData=phenoData(object1),
                annotation=annotation(object1),
                experimentData=experimentData(object1),
                predictions=rbind(predictions(object1), predictions(object2)),
                stateNames=stateNames(object1),
                chromosomeArm="both",
                distance=distance(object1),
                copyNumberIce=copyNumberIce(object1),
                transitionProbability=NULL,
                likelihood=NULL,
                log2=object1@log2,
                initialStateProbability=initialStateProbability(object1),
                locationCopyNumber=locationCopyNumber(object1),
                scaleCopyNumber=scaleCopyNumber(object1))
          })

setMethod("bind", c("NULL", "HmmSnpCopyNumberSet"), function(object1, object2, ...) object2)
setMethod("bind", c("HmmSnpCopyNumberSet", "NULL"), function(object1, object2, ...) object1)

##If there are more than one sample in object, it uses only the first.
setMethod("calculateEmissionProbability", "HmmSnpCopyNumberSet",
          function(object, ...){
            if(ncol(object) > 1) warning("calculateEmissionProbability uses the first sample in the object")
            object <- object[, 1]
            x <- as.vector(copyNumber(object))

            if(object@log2) x <- log2(x)

            locationCopyNumber(object) <- getCopyNumberLocation(object)
            B.cn <- locationCopyNumber(object)
            names(B.cn) <- stateNames(object)            

            if(!copyNumberIce(object)) {
              s <- scaleCopyNumber(object)
              emission.cn <- cbind(dnorm(x, B.cn[stateNames(object)[1]], s[1]),
                                   dnorm(x, B.cn[stateNames(object)[2]], s[2]),
                                   dnorm(x, B.cn[stateNames(object)[3]], s[3]))
            } else {
              s <- as.vector(cnConfidence(object))
              if(any(is.na(s))) stop("NA's in confidence scores")
              emission.cn <- cbind(dnorm(x, B.cn[stateNames(object)[1]], s),
                                   dnorm(x, B.cn[stateNames(object)[2]], s),
                                   dnorm(x, B.cn[stateNames(object)[3]], s))
            }
            emissionProbability(object) <- emission.cn
            rownames(emissionProbability(object)) <- featureNames(object)
            return(emissionProbability(object))
          })

setMethod("getCopyNumberLocation", "HmmSnpCopyNumberSet",
          function(object){
            if(is.null(locationCopyNumber(object))){
              isIce <- method(object) == "ICE"
              mu <- c(1, 2, 3)
              names(mu) <- stateNames(object)

              if(object@log2) mu <- log2(mu)
              if(annotation(object) == "Illumina500k"){
                print("Using ratios for the copy number location in the Illumina platform")
                mu <- c(log(1/2), 0, log(2))
              }
            } else { mu <- locationCopyNumber(object)}
            mu
          })

setMethod("getCopyNumberScale", "HmmSnpCopyNumberSet",
          function(object, I=1){
            if(is.null(scaleCopyNumber(object))){
              x <- copyNumber(object)[, I]
              x <- x[!is.na(x) & !is.nan(x)]
              if(object@log2)  x <- log2(x)
              x <- x[chromosome(object) != "X"]
              s <- diff(quantile(x, probs=c(0.16, (1-0.16)), na.rm=TRUE))/2
              s <- rep(s, numberStates(object))
              names(s) <- stateNames(object)
            } else {s <- scaleCopyNumber(object)}
            s
          })

setMethod("getInitialStateProbability", "HmmSnpCopyNumberSet", function(object){
  if(is.null(initialStateProbability(object))){
    isp <- c((1-0.99)/2, 0.99, (1-0.99)/2)
  } else isp <- initialStateProbability(object)
  isp
})

setMethod("getTransitionProbability", "HmmSnpCopyNumberSet",
          function(object){

            if(!is.null(transitionProbability(object))) return(transitionProbability(object))
            
            if(distance(object)){
              object2 <- object[order(position(object)), ]
              notTheSameOrder <- !(identical(featureNames(object2), featureNames(object)))
              if(notTheSameOrder) stop("SNPs must be ordered for distance calculations of transition probabilities")

              ##distance in units of 100 Mb ("~ 1cM") as in Beroukeim 2006 
              d <- (position(object)[2:nrow(object)] - position(object)[1:(nrow(object)-1)])/(100*10^6)
              theta <- 1-exp(-2*d)
              ##A is the probability of staying in the same state; it is 
              ##a vector of length T - 1
              A <- 1-theta
            } else {
              A <- matrix(0.01/3, ncol=3, nrow=3)
              diag(A) <- 1-0.01/3
              stateNames <- c("Deletion", "Normal", "Amplification")
              rownames(A) <- paste("current", stateNames)
              colnames(A) <- paste("next", stateNames)      
            }
            A
          })

##object may contain multiple chromosomes and multiple samples
setMethod("hmm", "HmmSnpCopyNumberSet",
          function(object,
                   verbose=TRUE,
                   arm=c("both", "p", "q")[1],
                   MISSING.CODE=NA,
                   keepEmission=FALSE,
                   keepLikelihood=FALSE,
                   keepTransitionProb=FALSE, ...){
            if(any(chromosome(object) == "X")){
              print("Chromosome X, XY, and M dropped from object")
              object <- object[chromosome(object) != "X" & chromosome(object) != "M" & chromosome(object) != "XY", ]
            }

            ##assumed to be the same for all samples/chromosomes
            data(chromosomeAnnotation, package="SNPchip", envir=environment())
            chromosomeArm(object) <- arm

            if(object@log2){
              if(min(copyNumber(object), na.rm=TRUE) < 0){
                warning("Negative values in copyNumber.  Not taking log")
                object@log2 <- FALSE
              }
            }
            locationCopyNumber(object) <- getCopyNumberLocation(object)
            scaleCopyNumber(object) <- getCopyNumberScale(object)            
            initialStateProbability(object) <- getInitialStateProbability(object)
            if(length(distance(object)) == 0) distance(object) <- TRUE
            objList <- split(object, as.character(chromosome(object)))

            ##Order each chromosome by physical position (important if genetic
            ##distance is used for transition probabilities)
            ##(chromosome-level ordering)
            myOrder <- function(object) object[order(position(object)), ]
            objList <- lapply(objList, myOrder)
            for(i in 1:length(objList)){
              if(verbose) print(paste("Fitting HMM to chromosome", names(objList)[i]))
              objList[[i]] <- hmmBySample(objList[[i]], 
                                          verbose=verbose,
                                          MISSING.CODE=MISSING.CODE,
                                          keepEmission=keepEmission,
                                          keepLikelihood=keepLikelihood,
                                          keepTransitionProb=keepTransitionProb)
            }
            if(length(objList) <= 1){
              object <- objList[[1]] ##else, return a list
            } else {
              object <- .unsplitHmm(objList, featureData(object), 
                                  copyNumber=do.call("rbind", lapply(object, copyNumber)),
                                  cnConfidence=do.call("rbind", lapply(object, cnConfidence)),
                                  predictions=do.call("rbind", lapply(object, predictions)),
                                  phenoData=phenoData(object[[1]]),
                                  annotation=annotation(object[[1]]),
                                  experimentData=experimentData(object[[1]]),
                                  initialStateProbability=intialStateProbability(object[[1]]),
                                  stateNames=stateNames(object[[1]]),
                                  copyNumberIce=copyNumberIce(object[[1]]),
                                  locationCopyNumber=locationCopyNumber(object[[1]]),
                                  distance=distance(object[[1]]))
            }
            object
          })

setMethod("locationCopyNumber", "HmmSnpCopyNumberSet", function(object) object@locationCopyNumber)
setReplaceMethod("locationCopyNumber", "HmmSnpCopyNumberSet",
                 function(object, value) {
                   object@locationCopyNumber <- value
                   object
                 })



setMethod("scaleCopyNumber", "HmmSnpCopyNumberSet", function(object) object@scaleCopyNumber)
setReplaceMethod("scaleCopyNumber", "HmmSnpCopyNumberSet",
                 function(object, value) {
                   object@scaleCopyNumber <- value
                   object
                 })

setMethod("show", "HmmSnpCopyNumberSet", function(object){
  cat(class( object ), " (storageMode: ", storageMode(object), ")\n", sep="")
  adim <- dim(object)
  if (length(adim)>1)
    cat("assayData:",
        if (length(adim)>1)
        paste(adim[[1]], "features,",
              adim[[2]], "samples") else NULL,
        "\n")
  cat("  element names:", paste(assayDataElementNames(object), collapse=", "), "\n")
  
  #######################################################
  ##General
  #######################################################
  cat("experimentData: use 'experimentData(object)'\n")
  pmids <- pubMedIds(object)
  if (length(pmids) > 0 && all(pmids != ""))
    cat("  pubMedIds:", paste(pmids, sep=", "), "\n")
  cat("phenoData\n")
  if(length(adim) > 1) show(phenoData(object))
  cat("featureData\n")      
  if(length(adim) > 1) show(featureData(object))
  cat("Annotation: ", annotation(object), "\n")
  cat("Hidden States: ", stateNames(object), "\n")
  cat("distance:", distance(object), "\n")
  cat("initialStateProbability: ", initialStateProbability(object), "\n")
  cat("chromosomeArm: ", chromosomeArm(object), "\n")


  #######################################################
  ##Class specific
  cat("locationCopyNumber: ", locationCopyNumber(object), "\n")
  cat("scaleCopyNumber: ", scaleCopyNumber(object), "\n")
  cat("log2: ", object@log2, "\n")
  cat("copyNumberIce: " , copyNumberIce(object), "\n")
})




