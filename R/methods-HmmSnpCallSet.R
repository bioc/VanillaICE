setMethod("initialize", "HmmSnpCallSet",
          function(.Object,
                   assayData = assayDataNew(
                     calls=calls,
                     callsConfidence=callsConfidence,
                     predictions=predictions),
                   featureData = annotatedDataFrameFrom(assayData, byrow=TRUE),
                   phenoData = annotatedDataFrameFrom(assayData, byrow=FALSE),
                   calls=new("matrix"),
                   callsConfidence=new("matrix"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   predictions=new("matrix"),
                   stateNames=character(),
                   initialStateProbability=NULL,
                   transitionProbability=NULL,
                   emissionProbability=NULL,
                   likelihood=NULL,
                   callsIce=logical(),
                   distance=logical(),
                   pCHOM=c(1-1e-3, 0.7),
                   ...) {
            assayData(.Object) <- assayData
            featureData(.Object) <- featureData
            phenoData(.Object) <- phenoData
            experimentData(.Object) <- experimentData
            annotation(.Object) <- annotation
            predictions(.Object) <- predictions
            stateNames(.Object) <- stateNames
            initialStateProbability(.Object) <- initialStateProbability
            transitionProbability(.Object) <- transitionProbability
            emissionProbability(.Object) <- emissionProbability
            likelihood(.Object) <- likelihood
            callsIce(.Object) <- callsIce
            distance(.Object) <- distance
            pCHOM(.Object) <- pCHOM
            .Object
          })

setValidity("HmmSnpCallSet", function(object){
  assayDataValidMembers(assayData(object), c("calls", "callsConfidence", "predictions"))
})

##Copied from Biobase/R/methods-eSet.R
##methods used to subset assayData are not exported from Biobase
setMethod("[", "HmmSnpCallSet", function(x, i, j, ..., drop = FALSE) {
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
    calls(x) <- calls(x)[i, , drop=FALSE]
    callsConfidence(x) <- callsConfidence(x)[i, , drop=FALSE]
    predictions(x) <- predictions(x)[i, , drop=FALSE]            
  }
  if(!missing(j)){
    calls(x) <- calls(x)[, j, drop=FALSE]
    callsConfidence(x) <- callsConfidence(x)[, j,drop=FALSE]
    predictions(x) <- predictions(x)[, j, drop=FALSE]            
  }
  x
})

setMethod("bind", c("HmmSnpCallSet", "HmmSnpCallSet"),
          function(object1, object2, ...){
            fd <- rbind(fData(object1), fData(object2))
            fD <- new("AnnotatedDataFrame", data=fd,
                      varMetadata=varMetadata(featureData(object1)))
            new("HmmSnpCallSet",
                calls=rbind(calls(object1), calls(object2)),
                callsConfidence=rbind(callsConfidence(object1), callsConfidence(object2)),
                featureData=fD,
                phenoData=phenoData(object1),
                annotation=annotation(object1),
                experimentData=experimentData(object1),
                predictions=rbind(predictions(object1), predictions(object2)),
                stateNames=stateNames(object1),
                chromosomeArm="both",
                distance=distance(object1),
                transitionProbability=NULL,
                likelihood=NULL,
                initialStateProbability=initialStateProbability(object1),
                pCHOM=pCHOM(object1),
                callsIce=callsIce(object1))
          })

setMethod("bind", c("NULL", "HmmSnpCallSet"), function(object1, object2, ...) object2)
setMethod("bind", c("HmmSnpCallSet", "NULL"), function(object1, object2, ...) object1)


##If there are more than one sample in object, it uses only the first.
setMethod("calculateEmissionProbability", "HmmSnpCallSet",
          function(object, ...){
            if(ncol(object) > 1) warning("calculateEmissionProbability uses the first sample in the object")
            object <- object[, 1]
            v <- as.vector(calls(object))
            v[v == 3] <- 1
            p <- as.vector(callsConfidence(object))
            
            if(!callsIce(object)) {
              emission.call <- getLohProbability(object)
              emission.call <- t(emission.call[, v])
            } else {
              s <- as.vector(cnConfidence(object))
              if(any(is.na(s))) stop("NA's in confidence scores")
              pCHOM <- pCHOM(object)
              emission.call <- callEmission(object,
                                            P.CHOM.LOH=pCHOM[1],                                            
                                            P.CHOM.Normal=pCHOM[2],
                                            SAMPLE=1)
              colnames(emission.call) <- stateNames(object)
            }
            emissionProbability(object) <- emission.call
            rownames(emissionProbability(object)) <- featureNames(object)
            return(emissionProbability(object))
          })



setMethod("getInitialStateProbability", "HmmSnpCallSet", function(object){
  if(is.null(initialStateProbability(object))){
    isp <- c(1-0.99, 0.99)
  } else isp <- initialStateProbability(object)
  isp
})


setMethod("getLohProbability", "HmmSnpCallSet",
          function(object){
            p <- pCHOM(object)
            B.call <- cbind(p, 1-p)            
            rownames(B.call) <- stateNames(object)
            colnames(B.call) <- c("P(CHOM|q)", "P(CHET|q)")
            B.call
          })

setMethod("getTransitionProbability", "HmmSnpCallSet",
          function(object){
            ##check if NULL
            if(!is.null(transitionProbability(object))) return(transitionProbability(object))
            
            if(distance(object)){
              object2 <- object[order(position(object)), ]
              notTheSameOrder <- !(identical(featureNames(object2), featureNames(object)))
              if(notTheSameOrder) stop("SNPs must be ordered for distance calculations of transition probabilities")

              ##distance in units of 100 Mb ("~ 1cM") as in Beroukeim 2006 
              d <- (position(object)[2:nrow(object)] - position(object)[1:(nrow(object)-1)])/(100*10^6)
              theta <- 1-exp(-2*d)
              A <- 1-theta
            } else {
              A <- matrix(c(0.99, 0.01,
                            0.01, 0.99), ncol=2, byrow=TRUE)
              rownames(A) <- paste("current", stateNames(object))
              colnames(A) <- paste("next", stateNames(object))      
            }
            A
          })

##object may contain multiple chromosomes and multiple samples
setMethod("hmm", "HmmSnpCallSet",
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
            if(is.null(initialStateProbability(object))){
              initialStateProbability(object) <- getInitialStateProbability(object)
            }
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
              object <- .unsplitHmm(objList, featureData(objList),
                                  calls=do.call("rbind", lapply(objList, calls)),
                                  callsConfidence=do.call("rbind", lapply(objList, calls)),
                                  predictions=do.call("rbind", lapply(objList, predictions)),                       
                                  experimentData=experimentData(objList[[1]]),
                                  phenoData=phenoData(objList[[1]]),
                                  annotation=annotation(objList[[1]]),
                                  stateNames=stateNames(objList[[1]]),
                                  initialStateProbability=initialStateProbability(objList[[1]]),
                                  callsIce=callsIce(objList[[1]]),
                                  pCHOM=pCHOM(objList[[1]]),
                                  distance=distance(objList[[1]]))
            }
            object
          })

setMethod("pCHOM", "HmmSnpCallSet", function(object) object@pCHOM)
setReplaceMethod("pCHOM", "HmmSnpCallSet",
                 function(object, value) {
                   if(length(value) != numberStates(object)) stop("replacement value should be a numerical vector of probabilities for each hidden state")
                   object@pCHOM <- value
                   object
                 })

##setMethod("plotHmm", "HmmSnpCallSet",
##          function(object,
##                   col=c("lightblue", "red"),
##                   pch=c("o", "x"),
##                   cex=c(1, 1.4),
##                   col.states="black",
##                   lwd=2,
##                   lty=1,
##                   ylim=NULL,
##                   amount=0.1,
##                   bty="o",
##                   xaxt="s",
##                   xaxs="r",
##                   yaxt="s",
##                   type="s",
##                   log="y",
##                   outer=FALSE,
##                   cex.axis=1.2,
##                   side.xaxis=1,
##                   blackAndWhite=FALSE,
##                   col.predict=NULL,
##                   ybottom=ylim[1],
##                   ytop=1,
##                   add=FALSE,
##                   breaks,
##                   legend=FALSE,
##                   col.centromere="bisque",
##                   cytoband=TRUE,
##                   useLayout=TRUE,
##                   ...){
##            if(cytoband){
##              if(useLayout){
##                par(las=1, mar=c(0.5, 4, 1, 1), oma=c(4, 1, 2, 1))
##                layout(matrix(1:2, 2, 1), heights=c(1, 0.1))
##              }
##              xaxt <- "n"
##            }            
##            if(length(col) == 2) col <- col[calls(object)]
##            if(length(cex) == 2) cex <- cex[calls(object)]
##            if(length(pch) == 2) pch <- pch[calls(object)]
##            if(is.null(col.predict)){
##              col.predict <- c("yellow", "white")
##            }
##            data(chromosomeAnnotation, package="SNPchip", envir=environment())
##            centromere <- chromosomeAnnotation[as.character(unique(chromosome(object))), 1:2]            
##            if(!add){
##              calls(object)[calls(object) == 2] <- 5
##              calls(object)[calls(object) == 1] <- 4
##              plot(position(object), jitter(calls(object), amount=amount),
##                   ylab="",
##                   xlab="Mb",
##                   xaxt="n",
##                   yaxt="n",
##                   xaxs=xaxs,
##                   col=col,
##                   pch=pch,
##                   cex=cex,
##                   ylim=ylim,
##                   bty=bty,
##                   log=log)
##              if(xaxt == "s") axis(side=1, at=pretty(range(position(object))), labels=pretty(range(position(object)))/1e6, outer=outer, cex.axis=cex.axis)
##              if(yaxt == "s") axis(side=2, at=c(4, 5), labels=c("Hom", "Het"), cex.axis=cex.axis)
##              ##Draw centromere
##
##              rect(xleft=centromere[[1]], ybottom=ylim[1],
##                   xright=centromere[[2]], ytop=ylim[2],
##                   col=col.centromere,
##                   border=col.centromere)                            
##            }
##            if(missing(breaks)) {
##              print("Calculating breakpoints")
##              breaks <- breakpoints(object, removeNormal=FALSE)
##            }
##            if(!missing(breaks)){
##              rect(xleft=min(position(object)), ytop=ytop, ybottom=ybottom, xright=max(position(object)), col=col.predict[2])
##              for(i in 1:nrow(breaks)){
##                drawRect(start=breaks[i, "start"], last=breaks[i, "last"], object=object,
##                         ybottom=ybottom, ytop=ytop,
##                         col=col.predict)
##              }
##              rect(xleft=centromere[[1]], ybottom=ybottom,
##                   xright=centromere[[2]], ytop=ytop,
##                   col="white",
##                   border="white")
##              rect(xleft=min(position(object)), ytop=ytop, ybottom=ybottom, xright=max(position(object)))
##            }
##            if(legend){
##              legend("topleft", fill=col.predict, legend=stateNames(object), 
##                     bty="n", cex=1, ncol=2)
##            }
##            if(cytoband){
##              plotCytoband(object,
##                           xaxs="r",
##                           cytobandAxis=FALSE,
##                           xlim=range(position(object)))
##              at <- pretty(range(position(object)))
##              at[at == max(at)] <- max(position(object))
##              axis(1, at=at, labels=round(at/1e6, 0), cex.axis=1, outer=TRUE, xaxs="i")
##            }            
##          })

setMethod("predictions", "HmmSnpCallSet", function(object) object@predictions)
setReplaceMethod("predictions", "HmmSnpCallSet",
                 function(object, value) {
                   object@predictions <- value
                   object
                 })

setMethod("show", "HmmSnpCallSet", function(object){
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
  cat("pCHOM: ", pCHOM(object), "\n")
  cat("callsIce: ", callsIce(object), "\n")
})



