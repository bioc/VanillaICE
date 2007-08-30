setMethod("breakpoints", "hSet",
          function(object, by="MB", decreasing=TRUE, reorder=TRUE,
                   removeNormal=TRUE, verbose=FALSE, ...){

            
            tmp <- list()
            for(i in 1:ncol(object)){
              tmp[[i]] <- .breaks(object[, i])
              if(verbose){
                cat(sampleNames(object)[i], ", ")
              }
            }
            if(length(tmp) < 1) return()
            if(verbose) cat("\n")
            tmp <- do.call("rbind", tmp)
            if(removeNormal) {
              if(sum(tmp$hiddenState == stateNames(object)[2]) > 0){
                tmp <- tmp[tmp$hiddenState != stateNames(object)[2], , drop=FALSE]
              }
            }
            if(reorder){
              if(length(tmp) > 0){
                tmp <- tmp[order(tmp[, by], decreasing=decreasing), , drop=FALSE]
              }
            }
            ##Add information on cytobands
            data(cytoband, package="SNPchip", envir=environment())
            xBreaks <- tmp[tmp[, "hiddenState"] != stateNames(object)[2], , drop=FALSE]
            bands <- apply(xBreaks, 1, .locateCytoband, cytoband)
            tmp$cytoband <- NA
            tmp$cytoband[tmp[, "hiddenState"] != stateNames(object)[2]] <- bands
            tmp
          })


setMethod("getBreaks", "hSet",
          function(object, breaks){
            if(missing(breaks)) {
              if(length(unique(chromosome(object))) <= 1){
                breaks <- breakpoints(object,
                                      removeNormal=FALSE,
                                      by="start",
                                      decreasing=FALSE,
                                      verbose=FALSE)
              } else{
                obj <- split(object, chromosome(object))
                breaks <- lapply(obj, breakpoints, removeNormal=FALSE, decreasing=FALSE)
                breaks <- do.call("rbind", breaks)
              }
            }
            breaks
          })


##setMethod("breakpoints", "list",
##          function(object, by="MB", decreasing=TRUE, reorder=FALSE, verbose=FALSE, excludeNormal, ...){
##            tmp <- lapply(object, breakpoints, verbose=verbose, reorder=FALSE)
##            tmp <- do.call("rbind", tmp)
##            if(length(tmp) > 0){
##              tmp <- tmp[tmp$hiddenState != stateNames(object[[1]])[2], ]
##            }
##            if(reorder){
##              if(length(tmp) > 0){
##                tmp <- tmp[order(tmp[, by], decreasing=decreasing), ]
##              }
##            }
##            tmp
##          })

setMethod("chromosomeArm", "hSet", function(object) object@chromosomeArm)
setReplaceMethod("chromosomeArm", "hSet",
                 function(object, value) {
                   object@chromosomeArm <- value
                   object
                 })

setMethod("distance", "hSet", function(object) object@distance)
setReplaceMethod("distance", "hSet",
                 function(object, value) {
                   object@distance <- value
                   object
                 })

setMethod("emissionProbability", "hSet", function(object) object@emissionProbability)
setReplaceMethod("emissionProbability", "hSet",
                 function(object, value) {
                   object@emissionProbability <- value
                   object
                 })

setMethod("getPar", "hSet", function(object, ...){
  op <- switch(class(object),
               HmmSnpSet=new("ParHmmSnpSet", ...),
               HmmSnpCallSet=new("ParHmmSnpCallSet", ...),
               HmmSnpCopyNumberSet=new("ParHmmSnpCopyNumberSet", ...))

  object <- object[!is.na(chromosome(object)), ]
  chromosomeNames <- unique(chromosome(object))
  chromosomeNames <- chromosomeNames[order(as.numeric(chromosomeNames))]
  op$firstChromosome <- chromosomeNames[1]
  N <- length(chromosomeNames)
  S <- ncol(object)

  if(op$add.cytoband){
    ## plots: data, predictions, data, predictions, ... data, predictions, cytoband
    ## == S * 2 + 1
    op$heights <- c(rep(c(1, 0.2), S), 0.2)
    S <- S+1
  } else op$heights <- rep(c(1, 0.2), S)
  data(chromosomeAnnotation, package="SNPchip", envir=environment())
  w <- chromosomeAnnotation[chromosomeNames, "chromosomeSize"]
  op$widths <- w/min(w)

  R <- length(op$heights)
  C <- length(op$widths)
  nPlots <- R*C
  if(N*S > 1){
    ##plot all the data and cytobands, then plot predictions
    ##the data
    m <- matrix(1:(N*S), ncol=N, nrow=S, byrow=FALSE)
    ##Add rows for predictions
    if(op$add.cytoband) p <- matrix(1:(N*ncol(object)), ncol=N, nrow=ncol(object), byrow=FALSE)
    p <- p+max(m)
    ##combine the matrices
    mat <- matrix(NA, R, C)
    k <- 1; i <- 1
    while(i < R){
      mat[i, ] <- m[k, ]
      mat[i+1, ] <- p[k, ]
      k <- k+1
      i <- i+2      
    }
    mat[R, ] <- m[nrow(m), ]
    op$mat <- mat
  }
  op$ylim <- .calculateYlim(object=object, op=op)
##  if(op$use.chromosome.size){
##    op$xlim <- c(0, chromosomeSize(chromosomeNames))
##  } else {
##    op$xlim <- range(position(object))
##  }
 if(ncol(object) == 1) op$yaxt="s"
  op
})

setMethod("plotSnp", "hSet",
          function(object, op, breaks, ...){
            ##check if predictions is empty (a reason to keep this separate from assayData)
            ##define a class that extends assayData
            ## call it analyticData?
            old.par <- par(no.readonly=TRUE)
            on.exit(par(old.par))
            ##plot the data and cytoband as usual
            ##Note: need to set on.exit to FALSE
            if(missing(op)) op <- getPar(object, ...)
            callNextMethod(object=object, op=op, on.exit=FALSE, ...)
            plotPredictions(object, op, breaks)
          })





setMethod("initialStateProbability", "hSet", function(object) object@initialStateProbability)
setReplaceMethod("initialStateProbability", "hSet",
                 function(object, value) {
                   object@initialStateProbability <- value
                   object
                 })

setMethod("likelihood", "hSet", function(object) object@likelihood)
setReplaceMethod("likelihood", "hSet",
                 function(object, value) {
                   object@likelihood <- value
                   object
                 })

setMethod("numberStates", "hSet", function(object) length(stateNames(object)))

setMethod("stateNames", "hSet", function(object) object@stateNames)
setReplaceMethod("stateNames", "hSet",
                 function(object, value) {
                   object@stateNames <- value
                   object
                 })


setMethod("transitionProbability", "hSet", function(object) object@transitionProbability)
setReplaceMethod("transitionProbability", "hSet",
                 function(object, value) {
                   object@transitionProbability <- value
                   object
                 })


##Would be good to write this in C++
setMethod("viterbi", "hSet",
          function(object, SCALE=1, ...){
            N.STATES <- numberStates(object)
            delta <- psi <- matrix(NA, nrow(object), numberStates(object))
            rownames(delta) <- rownames(psi) <- featureNames(object)
            colnames(delta) <- colnames(psi) <- stateNames(object)
            delta[1, ] <- log(initialStateProbability(object) * emissionProbability(object)[1, ])
            psi[1, ] <- rep(0, numberStates(object))
            A <- transitionProbability(object)
            if(!distance(object)){
              stopifnot(is.matrix(A))
            }
            T <- nrow(object)
            ##t loops over SNP 2 - T
            for(t in 2:T){
              ##j loops over states 1 - N
              if(distance(object)){
                AA <- matrix((1-A[t-1])*1/(N.STATES-1), nr=N.STATES, nc=N.STATES)
                diag(AA) <- A[t-1]
                if(SCALE != 1){
                  AA[lower.tri(AA)] <- AA[upper.tri(AA)] <- AA[upper.tri(AA)]*SCALE
                  diag(AA) <- 1-sum(AA[1, 2:N.STATES])
                }
##                if((sum(rowSums(AA) - N.STATES) > 1e-3)) stop("transition probabilities must sum to 1")
              } else AA <- A                                          
              for(j in 1:N.STATES){
                ##Equation 105b
                delta[t, j] <- max(delta[t-1, ]+log(AA[, j]))+log(emissionProbability(object)[t, j])
                psi[t, j] <- order(delta[t-1, ]+log(AA[, j]), decreasing=TRUE)[1]
              }
            }
            likelihood(object) <- delta
            Pstar <- max(delta[nrow(delta), ])
            qhat <- rep(NA, nrow(delta))
            T <- nrow(delta)
            names(qhat) <- rownames(delta)  
            qhat[T] <- order(delta[T, ], decreasing=TRUE)[1]
            for(t in (T-1):1) qhat[t] <- psi[t+1, qhat[t+1]]
            qhat
          })

##hmmSample iterates through each sample in the object
##The object may contain only one chromosome, but multiple samples
##**note that this is general to all HmmClasses
setMethod("hmmBySample", "hSet",
          function(object,
                   verbose=FALSE,
                   MISSING.CODE=4,
                   keepEmission=FALSE,
                   keepLikelihood=FALSE,
                   keepTransitionProb=FALSE,
                   SCALE=1, ...){
            ##chromosome-level
            isParm <- chromosomeArm(object) == "p"
            isQarm <- chromosomeArm(object) == "q"
            isBoth <- chromosomeArm(object) == "both"
            if(isParm + isQarm + isBoth != 1) stop("Must specify 'p', 'q', or 'both' for chromosomeArm")
            ##sample-level
            ##Exclude missing calls
            data(chromosomeAnnotation, package="SNPchip", envir=environment())
            if(isParm | isBoth){
              parm <- object[position(object) < chromosomeAnnotation[as.character(chromosome(object)), 2], ]
              if(nrow(parm) < 1) {
                parm <- NULL
                if(isParm) stop("No SNPs on q arm.  Specify different arm.")
              } else {
                parm <- byArm(parm, verbose=verbose, MISSING.CODE=MISSING.CODE, arm="p", SCALE=SCALE)
              }
            }
            if(isQarm | isBoth){
              qarm <- object[position(object) > chromosomeAnnotation[as.character(chromosome(object)), 2], ]
              if(nrow(qarm) < 1) {
                q <- NULL
                if(isQarm) stop("No SNPs on q arm.  Specify different arm.")
              } else {
                qarm <- byArm(qarm, verbose=verbose, MISSING.CODE=MISSING.CODE, arm="q", SCALE=SCALE)
              }              
            }
            if(isParm) object <- parm
            if(isQarm) object <- qarm
            if(isBoth){
              ##combine q and p arms
              object <- bind(parm, qarm)
            }

            if(!keepEmission) emissionProbability(object) <- NULL
            if(keepEmission) warning("If fitting the HMM to multiple samples, only the emission probabilities in the last sample are saved")
            
            if(!keepLikelihood) likelihood(object) <- NULL
            if(keepLikelihood) warning("If fitting the HMM to multiple samples, only the likelihood from the last sample is saved")
            
            if(!keepTransitionProb) transitionProbability(object) <- NULL
            if(keepTransitionProb) warning("If fitting the HMM to multiple samples, only the transitionProb from the last sample is saved")
            object
          })

##return results from one arm
setMethod("byArm", "hSet",
          function(object, verbose, MISSING.CODE, arm="", SCALE=1, ...){
            for(j in 1:ncol(object)){
              excludeSample <- FALSE
              if(verbose) print(paste("   processing ", arm, " arm of sample", sampleNames(object)[j], " ..."))
              obj <- object[, j]
              if(class(obj) != "HmmSnpCopyNumberSet"){
                missing <- isMissing(obj, MISSING.CODE=MISSING.CODE)
                if(sum(!missing) < 5){
                  print("     **missing too many genotype calls...  skipping this sample**")
                  excludeSample <- TRUE
                } 
                obj <- obj[!missing, ]                
              } else missing <- rep(FALSE, nrow(obj))
              ##Calculate transition probabilities after subsetting
              if(!excludeSample){
                transitionProbability(obj) <- getTransitionProbability(obj)
                emissionProbability(obj) <- calculateEmissionProbability(obj)
                predictions(object)[!missing, j] <- viterbi(obj, SCALE=SCALE)
              } else {
                ##Exclude sample => return NA's for predictions
                predictions(object)[, j] <- rep(NA, nrow(object))
              }
            }
            object
          })




