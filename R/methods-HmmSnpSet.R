setMethod("initialize", "HmmSnpSet",
          function(.Object,
                   assayData = assayDataNew(
                     calls=calls,
                     callsConfidence=callsConfidence,
                     copyNumber=copyNumber,
                     cnConfidence=cnConfidence,
                     predictions=predictions),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   calls = new("matrix"),
                   callsConfidence = new("matrix"),
                   copyNumber=new("matrix"),
                   cnConfidence=new("matrix"),
                   predictions=new("matrix"),
                   scaleCopyNumber=numeric(),
                   pCHOM=c(1-1e-3, 0.7, 1-1e-3, 0.7),
                   log2=logical(),
                   locationCopyNumber=numeric(),
                   stateNames=character(),
                   initialStateProbability=NULL,
                   transitionProbability=NULL,
                   emissionProbability=NULL,
                   likelihood=NULL,
                   copyNumberIce=FALSE,
                   callsIce=FALSE,
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
            pCHOM(.Object) <- pCHOM
            initialStateProbability(.Object) <- initialStateProbability
            transitionProbability(.Object) <- transitionProbability
            emissionProbability(.Object) <- emissionProbability
            likelihood(.Object) <- likelihood
            copyNumberIce(.Object) <- copyNumberIce
            callsIce(.Object) <- callsIce
            distance(.Object) <- distance
            .Object
          })

setValidity("HmmSnpSet", function(object){
  assayDataValidMembers(assayData(object), c("calls", "callsConfidence", "copyNumber", "cnConfidence", "predictions"))
})

##Copied from Biobase/R/methods-eSet.R
##methods used to subset assayData are not exported from Biobase
setMethod("[", "HmmSnpSet", function(x, i, j, ..., drop = FALSE) {
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
    calls(x) <- calls(x)[i, , drop=FALSE]
    callsConfidence(x) <- callsConfidence(x)[i, , drop=FALSE]
    predictions(x) <- predictions(x)[i, , drop=FALSE]    
  }
  if(!missing(j)){
    copyNumber(x) <- copyNumber(x)[, j, drop=FALSE]
    cnConfidence(x) <- cnConfidence(x)[, j, drop=FALSE]
    calls(x) <- calls(x)[, j, drop=FALSE]
    callsConfidence(x) <- callsConfidence(x)[, j,drop=FALSE]
    predictions(x) <- predictions(x)[, j, drop=FALSE]    
  }  
  x
})

setReplaceMethod("[", "HmmSnpSet", function(x, i, j,..., value){
  if(class(value) != "HmmSnpSet") stop("Value to replace must be a HmmTrioSet")
  ##should this be in the initialization of HmmTrioSet?
##  predictions(x) <- as.integer(matrix(NA, nrow(x), ncol(x)))
  ##only for one sample -- typically NULL
  if(!is.null(likelihood(x))){
    likelihood(x) <- matrix(NA, nrow(x), numberStates(value))
  }

  ##only for one sample -- typically NULL
  if(!is.null(emissionProbability(x))){
    emissionProbability(x) <- matrix(NA, nrow=nrow(x), ncol=numberStates(value))
  }
  
  distance(x) <- as.logical(NA)
  transitionProbability(x) <- as.numeric(rep(NA, nrow(x)))
  if(!missing(i) & missing(j)){
    predictions(x)[i, ] <- predictions(value)
  }
  if(!missing(i) & !missing(j)){
    predictions(x)[i, j] <- predictions(value)[, j]
  }  
  if(!missing(i)){
    likelihood(x)[i, ] <- likelihood(value)
    emissionProbability(x)[i, ] <- emissionProbability(value)
    distance(x) <- distance(value)
    ##difficult to subset transition probability, but I don't think we need to
    transitionProbability(x) <- transitionProbability(value)
  }
  x
})

setMethod("bind", c("HmmSnpSet", "HmmSnpSet"),
          function(object1, object2, ...){
            fd <- rbind(fData(object1), fData(object2))
            fD <- new("AnnotatedDataFrame", data=fd,
                      varMetadata=varMetadata(featureData(object1)))
            new("HmmSnpSet",
                copyNumber=rbind(copyNumber(object1), copyNumber(object2)),
                cnConfidence=rbind(cnConfidence(object1), cnConfidence(object2)),
                calls=rbind(calls(object1), calls(object2)),
                callsConfidence=rbind(callsConfidence(object1), callsConfidence(object2)),
                featureData=fD,
                phenoData=phenoData(object1),
                annotation=annotation(object1),
                experimentData=experimentData(object1),
                predictions=rbind(predictions(object1), predictions(object2)),
                stateNames=stateNames(object1),
                chromosomeArm="both",
##                method=method(object1),
                copyNumberIce=copyNumberIce(object1),
                callsIce=callsIce(object1),
                distance=distance(object1),
                initialStateProbability=NULL,
                transitionProbability=NULL,
                likelihood=NULL,
                log2=object1@log2,
                locationCopyNumber=locationCopyNumber(object1),
                scaleCopyNumber=scaleCopyNumber(object1),
                pCHOM=pCHOM(object1))
          })

setMethod("bind", c("NULL", "HmmSnpSet"), function(object1, object2, ...) object2)
setMethod("bind", c("HmmSnpSet", "NULL"), function(object1, object2, ...) object1)

setMethod("breakpoints", "HmmSnpSet",
           function(object, by="MB", decreasing=TRUE, reorder=TRUE,
                    removeNormal=TRUE, verbose=FALSE, ...){
             tmp <- callNextMethod(object,
                                   by=by,
                                   decreasing=decreasing,
                                   reorder=reorder,
                                   removeNormal=removeNormal,
                                   verbose=verbose)
             ##number of SNPs in region
             notNormal <- tmp[, "hiddenState"] != stateNames(object)[2]
             if(sum(notNormal) < 1) return(NULL)
             xBreaks <- tmp[notNormal, ]
             numberSnps <- function(breaks, object){
               I <- position(object) >= as.numeric(breaks["start"]) & position(object) <= as.numeric(breaks["last"])
               j <- match(breaks["sampleId"], sampleNames(object))
               nSnps <- sum(I)
               pHet <- sum(calls(object)[I, j] == 2, na.rm=TRUE)
               pNA <- sum(is.na(calls(object)[I, j]))
               x <- c(nSnps, pHet, pNA)
             }
             nSnps <- t(apply(xBreaks, 1, numberSnps, object))
             colnames(nSnps) <- c("numberSnps", "CHET", "is.na")
             tmp$numberSnps <- NA
             tmp$CHET <- NA
             tmp$is.na <- NA

             tmp[notNormal, "numberSnps"] <- nSnps[, "numberSnps"]
             tmp[notNormal, "CHET"] <- nSnps[, "CHET"]
             tmp[notNormal, "is.na"] <- nSnps[, "is.na"]
             tmp
           })

##If there are more than one sample in object, it uses only the first.
setMethod("calculateEmissionProbability", "HmmSnpSet",
          function(object, ...){
            if(ncol(object) > 1) warning("calculateEmissionProbability uses the first sample in the object")
            object <- object[, 1]
            x <- as.vector(copyNumber(object))

            if(object@log2) x <- log2(x)
              
            v <- as.vector(calls(object))
            v[v == 3] <- 1
            p <- as.vector(callsConfidence(object))

            locationCopyNumber(object) <- getCopyNumberLocation(object)
            B.cn <- locationCopyNumber(object)
            names(B.cn) <- stateNames(object)            

            if(!callsIce(object)) {
              emission.call <- getLohProbability(object)
              emission.call <- t(emission.call[, v])
            }  else {
              ##check with confidence scores are available.  If not,
              ##use the vanilla emission probability for calls.
              pCHOM <- pCHOM(object)
              emission.call <- callEmission(object,
                                            P.CHOM.LOH=pCHOM[1], #probability CHOM in a deletion
                                            P.CHOM.Normal=pCHOM[2], ##probability CHOM in normal region 
                                            SAMPLE=1)

              emission.call2 <- callEmission(object,
                                            P.CHOM.LOH=pCHOM[3],  #probability CHOM in copy-neutral region
                                            P.CHOM.Normal=pCHOM[4], ##probability CHOM in an amplification
                                            SAMPLE=1)
              
              ##There are 4 states and b.call has 2 columns.  Replicate the
              ##first 2 columns.  The joint emission probability will be a
              ##product of copy number emission probability and call emission
              ##probability
              emission.call <- cbind(emission.call, emission.call2)
              colnames(emission.call) <- stateNames(object)
            }
              
            if(!copyNumberIce(object)){
              s <- scaleCopyNumber(object)
              emission.cn <- cbind(dnorm(x, B.cn[stateNames(object)[1]], s[1]),
                                   dnorm(x, B.cn[stateNames(object)[2]], s[2]),
                                   dnorm(x, B.cn[stateNames(object)[3]], s[3]),
                                   dnorm(x, B.cn[stateNames(object)[4]], s[4]))
            } else{
              ##Write wrapper that does this according to class
              s <- as.vector(cnConfidence(object))
              if(any(is.na(s))) stop("NA's in confidence scores")
              emission.cn <- cbind(dnorm(x, B.cn[stateNames(object)[1]], s),
                                   dnorm(x, B.cn[stateNames(object)[2]], s),
                                   dnorm(x, B.cn[stateNames(object)[3]], s),
                                   dnorm(x, B.cn[stateNames(object)[4]], s))

            }
            emissionProbability(object) <- emission.call*emission.cn
            rownames(emissionProbability(object)) <- featureNames(object)
            return(emissionProbability(object))
          })

setMethod("getLohProbability", "HmmSnpSet",
          function(object){
            p <- pCHOM(object)
            B.call <- cbind(p, 1-p)
            rownames(B.call) <- stateNames(object)
            colnames(B.call) <- c("P(HOM|q)", "P(HET|q)")
            B.call
          })

setMethod("getCopyNumberLocation", "HmmSnpSet",
          function(object){
            if(is.null(locationCopyNumber(object))){
              isIce <- method(object) == "ICE"
              mu <- c(1, 2, 2, 3)
              names(mu) <- stateNames(object)
              if(object@log2) mu <- log2(mu)
              if(annotation(object) == "Illumina500k"){
                print("Using ratios for the copy number location in the Illumina platform")
                mu <- c(log(1/2), 0, 0, log(2))
                names(mu) <- stateNames(object)
              }
            } else {
              mu <- locationCopyNumber(object)
              names(mu) <- stateNames(object)
            }
            mu
        })


setMethod("getCopyNumberScale", "HmmSnpSet",
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


setMethod("getInitialStateProbability", "HmmSnpSet",
          function(object){
            if(is.null(initialStateProbability(object))){            
              isp <- c((1-0.99)/3, 0.99, (1-0.99)/3, (1-0.99)/3)
            } else isp <- initialStateProbability(object)
            isp
          })

setMethod("getTransitionProbability", "HmmSnpSet",
          function(object){
            if(!is.null(transitionProbability(object))) return(transitionProbability(object))
            if(distance(object)){
              object2 <- object[order(position(object)), ]
              notTheSameOrder <- !(identical(featureNames(object2), featureNames(object)))
              if(notTheSameOrder) stop("SNPs must be ordered for distance calculations of transition probabilities")
              ##distance in units of 100 Mb ("~ 1cM") as in Beroukeim 2006 
              d <- (position(object)[2:nrow(object)] - position(object)[1:(nrow(object)-1)])/(100*10^6)
              ##SNP-specific transition probabilities.  We will assume
              ##that the probability of transition from one state to
              ##any other state is the same. Only 1 number is needed
              ##for snp t, d_t
              theta <- 1-exp(-2*d)
              ##A is the probability of staying in the same state; it is 
              ##a vector of length T - 1
              A <- 1-theta
            } else {
              A <- matrix(c(0.98, 0.01, 0.005, 0.005,
                            (1-0.999)/3, 0.999, (1-0.999)/3, (1-0.999)/3,
                            0.05, 0.01, 0.98, 0.05,
                            0.05, 0.01, 0.05, 0.98), ncol=4, byrow=TRUE)
              rownames(A) <- paste("current", stateNames(object))
              colnames(A) <- paste("next", stateNames(object))
            }
            A
          })


##object may contain multiple chromosomes and multiple samples
setMethod("hmm", "HmmSnpSet",
          function(object,
                   verbose=TRUE,
                   log2=TRUE,
                   arm=c("both", "p", "q")[1],
                   MISSING.CODE=NA,
                   keepEmission=FALSE,
                   keepLikelihood=FALSE,
                   keepTransitionProb=FALSE,
                   SCALE=1, ...){
            object@log2 <- log2
            if(any(chromosome(object) == "X")){
              print("Chromosome X, XY, Y, and M dropped from object")
              object <- object[chromosome(object) != "Y" & chromosome(object) != "X" & chromosome(object) != "M" & chromosome(object) != "XY", ]
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

            ##using all autosomes (X chromosome has already been excluded)
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
                                          keepTransitionProb=keepTransitionProb,
                                          SCALE=SCALE)
            }
            if(length(objList) <= 1){
              object <- objList[[1]] ##else, return a list
            } else {
              object <- .unsplitHmm(objList, featureData(object),
                                    copyNumber=do.call("rbind", lapply(objList, copyNumber)),
                                    cnConfidence=do.call("rbind", lapply(objList, cnConfidence)),
                                    calls=do.call("rbind", lapply(objList, calls)),
                                    callsConfidence=do.call("rbind", lapply(objList, calls)),                       
                                    predictions=do.call("rbind", lapply(objList, predictions)),
                                    phenoData=phenoData(objList[[1]]),
                                    annotation=annotation(objList[[1]]),
                                    experimentData=experimentData(objList[[1]]),
                                    initialStateProbability=initialStateProbability(objList[[1]]),
                                    stateNames=stateNames(objList[[1]]),
                                    copyNumberIce=copyNumberIce(objList[[1]]),
                                    locationCopyNumber=locationCopyNumber(objList[[1]]),
                                    distance=distance(objList[[1]]),
                                    pCHOM=pCHOM(objList[[1]]),
                                    callsIce=callsIce(objList[[1]]))
            }
            object
          })


setMethod("locationCopyNumber", "HmmSnpSet", function(object) object@locationCopyNumber)
setReplaceMethod("locationCopyNumber", "HmmSnpSet",
                 function(object, value) {
                   object@locationCopyNumber <- value
                   object
                 })

setMethod("pCHOM", "HmmSnpSet", function(object) object@pCHOM)
setReplaceMethod("pCHOM", "HmmSnpSet",
                 function(object, value) {
                   if(length(value) != numberStates(object)) stop("replacement value should be a numerical vector of probabilities for each hidden state")                   
                   object@pCHOM <- value
                   object
                 })

setMethod("predictions", "HmmSnpSet", function(object) assayDataElement(object, "predictions"))
setReplaceMethod("predictions", signature(object="HmmSnpSet", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "predictions", value))

setMethod("scaleCopyNumber", "HmmSnpSet", function(object) object@scaleCopyNumber)
setReplaceMethod("scaleCopyNumber", "HmmSnpSet",
                 function(object, value) {
                   object@scaleCopyNumber <- value
                   object
                 })

setMethod("show", "HmmSnpSet", function(object){
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
  cat("pCHOM: ", pCHOM(object), "\n")
  cat("log2: ", object@log2, "\n")
  cat("copyNumberIce: ", copyNumberIce(object), "\n")
  cat("callsIce: ", callsIce(object), "\n")  
})

setMethod("ice2VanillaICE", "HmmSnpSet",
          function(object){
            attr(class(object), "package") <- "VanillaICE"
            object
          })





