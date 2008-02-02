##One big matrix that could potentially span several chromosomes.
##featureData should be looked up in RSQLite database
setMethod("initialize", "HmmParameter",
          function(.Object,
                   snpset,
                   states=character(),
                   tau=matrix(),
                   tau.scale=matrix(),                   
                   pi=numeric(),
                   beta=matrix(),
                   notes=character(),
                   SCALE=2,##by default, twice as likely to return to normal state
                   ...){
            snpset <- snpset[!(is.na(chromosome(snpset))), ]
            chrom <- chromosome(snpset)
            chrom[chrom == "X"] <- 23
            chrom[chrom == "XY"] <- 24
            chrom[chrom == "M"] <- 25
            chrom <- as.numeric(chrom)            
            snpset <- snpset[order(chrom, position(snpset)), ]

##            browser()
            hmmOptions <- new("HmmOptions", snpset=snpset, states=states, ...)
            S <- length(hmmOptions@states)
            featureNames <- featureNames(snpset)

            data(chromosomeAnnotation, package="SNPchip", envir=environment())
            chrAnn <- as.matrix(chromosomeAnnotation)
            tmp <- position(snpset) <= chromosomeAnnotation[chromosome(snpset), "centromereStart"]
            chromosomeArm <- as.character(tmp)
            chromosomeArm[chromosomeArm == "TRUE"] <- "p"
            chromosomeArm[chromosomeArm == "FALSE"] <- "q"
            fD <- data.frame(cbind(chromosome(snpset), arm=chromosomeArm), row.names=featureNames,
                             stringsAsFactors=FALSE)
            colnames(fD) <- c("chromosome", "arm")
            varMetadata <- data.frame(labelDescription=c("chromosome", "arm"), row.names=colnames(fD))
            featureData <- new("AnnotatedDataFrame", data=fD, varMetadata=varMetadata)

            ###########################################################################            
            ##Initial State Probabilities
            ###########################################################################
            states <- hmmOptions@states
            if(length(pi) == 0){
              tmp <- rep(0, length(states))
              names(tmp) <- states
              tmp["N"] <- 0.999
              tmp[states != "N"] <- (1-0.999)/(length(states) - 1)
              pi <- log(tmp)
            }

            ###########################################################################            
            ##Emission Probabilities
            ###########################################################################
            ##Rows = number of SNPs * number of states, Columns = number of samples
            if(hmmOptions@SnpClass == "SnpCallSet" | hmmOptions@SnpClass == "oligoSnpSet"){
              x <- data.frame(calls(snpset), row.names=featureNames(snpset))
              confidence <- data.frame(callsConfidence(snpset), row.names=featureNames(snpset))
              print("Calculating emission probabilities for genotype calls.... ")
              beta.gt <- mapply(gtEmission, x=x, confidence=confidence,
                                MoreArgs=list(options=hmmOptions, featureNames=featureNames(snpset)), SIMPLIFY=TRUE)
            } else beta.gt <- 0

            if(hmmOptions@SnpClass == "SnpCopyNumberSet" | hmmOptions@SnpClass == "oligoSnpSet"){
              x <- data.frame(copyNumber(snpset), row.names=featureNames(snpset))
              confidence <- data.frame(cnConfidence(snpset), row.names=featureNames(snpset))
              cn.robustSE <- as.list(hmmOptions@cn.robustSE)
              print("Calculating emission probabilities for copy number estimates... ")
              beta.cn <- mapply(cnEmission, x=x, confidence=confidence,
                                robustSE=cn.robustSE,
                                MoreArgs=list(options=hmmOptions, featureNames=featureNames(snpset)), SIMPLIFY=TRUE)
            } else beta.cn <- 0

            ##Missing values in the emission probabilities are
            ##typically caused by missing values for genotype calls.
            ##Plugging in a constant such as zero for each of the
            ##hidden states should prevent the genotype call of the
            ##SNP from contributing to the likelihood, however the
            ##copy number estimate (if available) will still
            ##contribute
            beta.gt[is.na(beta.gt)] <- 0
            beta.cn[is.na(beta.cn)] <- 0
            beta <- beta.cn + beta.gt

            ###########################################################################            
            ##Transition Probabilities
            ###########################################################################                        
            ##Scaling transition probabilities between states
            if(missing(tau.scale) & S > 2){
              ##scalar parameters for the transition probability
              ##matrix.  By default, the probability of transitioning
              ##from an altered state back to a normal state is twice
              ##as likely as the probability of transitioning between
              ##two altered states              
              S <- length(states)
              tau.scale <- matrix(1, S, S)
              rownames(tau.scale) <- colnames(tau.scale) <- states
              x <- seq(1, (S-1), by=0.01)
              y <- rep(NA, length(x))
              for(i in 1:length(x)){
                y[i] <- ((S-1)-x[i])/(S-2) 
              }
              z <- abs(x/y - SCALE)
              z <- z == min(z)
              x <- x[z]
              y <- y[z]
              tau.scale[upper.tri(tau.scale)] <- y
              tau.scale[lower.tri(tau.scale)] <- y
              tau.scale[, "N"] <- x
              ##by default, do not tau.scale the probabilities of leaving the normal state
              tau.scale["N", ] <- 1
              diag(tau.scale) <- 1
            }
            if(S == 2) tau.scale <- matrix(1, nrow=S, ncol=S)
              
            ##We definine the transition probability matrix as a
            ##function of the distance between adjacent SNPs on the
            ##same chromosomal arm.  This requires breaking the object
            ##into a list by chromosome.  We then order the SNPs by
            ##physical position along the chromosome, calculating the
            ##distance between adjacent SNPs. If there are S states,
            ##the transition probability matrix from SNP t to SNP t+1
            ##has dimension S x S.  For each SNP t, we convert this
            ##matrix to a vector of length S^2.  We then return a
            ##matrix of the transition probabilities with R-1 rows and
            ##S^2 columns. For SNPs at the beginning or the end of a
            ##chromosome, the transition probability matrix will be
            ##NA's.  For SNPs at the beginning or end of a chromosomal
            ##arm, the transition probability matrix can be any number
            ##(it will not be used by the HMM).
            tau <- calculateTransitionProbability(object=snpset, options=hmmOptions, scale=tau.scale)
            if(any(is.na(tau))){
              warning("Missing values in transition probabilities")
              browser()
            }
            rownames(tau) <- paste(featureNames(snpset)[1:(nrow(snpset)-1)], featureNames(snpset)[2:nrow(snpset)], sep="->")
            .Object@hmmOptions <- hmmOptions
            .Object@tau <- tau
            .Object@tau.scale <- tau.scale
            .Object@pi <- pi
            .Object@beta <- beta
            .Object@featureData <- featureData
            .Object
          })

setMethod("hmmOptions", "HmmParameter", function(object) object@hmmOptions)
setMethod("featureData", "HmmParameter", function(object) object@featureData)
setMethod("featureNames", "HmmParameter", function(object) featureNames(featureData(object)))
##setMethod("nrow", "HmmParameter", function(x) length(featureNames(x)))
setMethod("pData", "HmmParameter", function(object) pData(featureData(object)))
setMethod("sampleNames", "HmmParameter", function(object) colnames(Beta(object)))
setMethod("chromosome", "HmmParameter", function(object) featureData(object)$chromosome)
setMethod("arm", "HmmParameter", function(object) featureData(object)$arm)
setMethod("Beta", "HmmParameter", function(object) object@beta)
setMethod("show", "HmmParameter", function(object) str(object))
setMethod("SnpClass", "HmmParameter", function(object) SnpClass(hmmOptions(object)))
setMethod("states", "HmmParameter", function(object) states(hmmOptions(object)))
setMethod("pi", "HmmParameter", function(object) object@pi)
setMethod("tau", "HmmParameter", function(object) object@tau)
setMethod("gt.ICE", "HmmParameter", function(object) gt.ICE(hmmOptions(object)))
setMethod("cn.ICE", "HmmParameter", function(object) cn.ICE(hmmOptions(object)))
setMethod("cn.location", "HmmParameter", function(object) cn.location(hmmOptions(object)))
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

          
                   
