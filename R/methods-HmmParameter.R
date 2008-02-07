setMethod("hmm", "HmmParameter",
          function(object, snpset, sn, ...){
            arm <- paste(chromosome(object), arm(object), sep="")
            beta <- data.frame(Beta(object))
            predictions <- apply(beta, 2, viterbi, pi=pi(object),
                                 tau=tau(object), states=states(object),
                                 arm=arm)
            rownames(predictions) <- featureNames(object)
            if(missing(sn)) sn <- sampleNames(snpset)
            colnames(predictions) <- sn

            ##  snpset <- snpset[!(is.na(chromosome(snpset))), ]
            snpset <- snpset[match(featureNames(object), featureNames(snpset)), ]
            breaks <- mapply(calculateBreakpoints, x=data.frame(predictions),
                             sampleNames=as.list(sn),
                             MoreArgs=list(position=position(snpset),
                               chromosome=chromosome(snpset),
                               states=states(object)),  ##for inheritance purposes
                             SIMPLIFY=FALSE)
            hmmOut <- new("HmmPredict",
                          states=states(object),
                          predictions=predictions,
                          breakpoints=breaks,
                          SnpClass=SnpClass(object),
                          featureData=featureData(object))
          })

setMethod("hmmOptions", "HmmParameter", function(object) object@hmmOptions)
setReplaceMethod("hmmOptions", c("HmmParameter", "HmmOptions"),
                 function(object, value){
                   object@hmmOptions <- value
                   object
                 })
setMethod("featureData", "HmmParameter", function(object) object@featureData)
setMethod("featureNames", "HmmParameter", function(object) featureNames(featureData(object)))
##setMethod("nrow", "HmmParameter", function(x) length(featureNames(x)))
setMethod("fData", "HmmParameter", function(object) pData(featureData(object)))
setMethod("sampleNames", "HmmParameter", function(object) colnames(Beta(object)))
setMethod("chromosome", "HmmParameter", function(object) featureData(object)$chromosome)
setMethod("arm", "HmmParameter", function(object) featureData(object)$arm)
setMethod("Beta", "HmmParameter", function(object) object@beta)
setMethod("show", "HmmParameter", function(object) str(object))
setMethod("SnpClass", "HmmParameter", function(object) SnpClass(hmmOptions(object)))
setMethod("states", "HmmParameter", function(object) states(hmmOptions(object)))
setReplaceMethod("states", c("HmmParameter", "character"),
                 function(object, value){
                   hmmOptions(object)@states <- value
                   object
                 })
setMethod("pi", "HmmParameter", function(object) object@pi)
setMethod("tau", "HmmParameter", function(object) object@tau)
setMethod("gt.ICE", "HmmParameter", function(object) gt.ICE(hmmOptions(object)))
setMethod("cn.ICE", "HmmParameter", function(object) cn.ICE(hmmOptions(object)))
setMethod("cn.location", "HmmParameter", function(object) cn.location(hmmOptions(object)))
setMethod("cn.robustSE", "HmmParameter", function(object) cn.robustSE(hmmOptions(object)))
setMethod("tau.scale", "HmmParameter", function(object) object@tau.scale)

setMethod("[", "HmmParameter",
          function(x, i, j, ..., drop = FALSE){
            if (missing(drop)) drop <- FALSE
            if (missing(i) && missing(j)) {
              if (length(list(...))!=0)
                stop("specify genes or samples to subset; use '",
                     substitute(x), "$", names(list(...))[[1]],
                     "' to access phenoData variables")
              return(x)
            }
            if(length(j) > 1 | missing(j)){
              stop("must specify 1 column (sample) to subset")
            }
            if(length(i) > 1 | missing(i)){
              stop("must specify 1 SNP to subset")
            }
            ##number of rows (SNPs)
            R <- length(featureNames(x))
            ##number of states 
            S <- length(pi(x))
            ##number of columns (samples)
            C <- ncol(Beta(x))
            if(!missing(i)){
              beta <- Beta(x)[, j]
              beta <- matrix(beta, nrow=R, ncol=S)
              beta <- beta[i, , drop=FALSE]
              colnames(beta) <- states(x)
              featureData <- featureData(x)[i, j]              
              rownames(beta) <- featureNames(featureData)

              transition <- list()
              states <- states(x)
              transition[[1]] <- matrix(tau(x)[(i-1), ], nrow=S, ncol=S)
              colnames(transition[[1]]) <- paste(states, "t+1", sep="_")
              rownames(transition[[1]]) <- paste(states, "t", sep="_")
              transition[[2]] <- matrix(tau(x)[(i-1), ], nrow=S, ncol=S)
              colnames(transition[[2]]) <- paste(states, "t+1", sep="_")
              rownames(transition[[2]]) <- paste(states, "t", sep="_")
              names(transition) <- rownames(tau(x))[(i-1):i]
            }
            results <- list()
            results[[1]] <- list(tau=transition, beta=beta)
            names(results) <- featureNames(featureData)
            results
          })

          
                   
