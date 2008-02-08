##Default:  ordinary visualization of data with plotSnp for one sample / one chromosome
##Add a row for the prediction bar
##Add a row for the cytoband
setMethod("initialize", "ParHSet",
          function(.Object,
                   border.prediction=grey(0.4),
                   legend.prediction=TRUE,
                   legend.location.prediction="topleft", ...){
            .Object <- callNextMethod(.Object, ...)
            if(.Object$useLayout){
              if(.Object$add.cytoband){
                .Object$heights <- c(1, 0.1, 0.05)
                ##1 - plot the data
                ##2 - plot the cytoband
                ##3 - plot the predictions.  Plot the predictions
                ##between the data and the cytoband
                .Object$mat <- matrix(c(1, 3, 2), ncol=1)                
              } else {
                ##do not plot the cytoband
                .Object$heights <- c(1, 0.2)
                .Object$mat <- matrix(1:2, ncol=1)                                
              }
            }
            snpPar(.Object)$border.prediction <- border.prediction
            snpPar(.Object)$legend.prediction <- legend.prediction
            snpPar(.Object)$legend.location.prediction <- legend.location.prediction
            snpPar(.Object)$predictions.ylim <- c(-0.2, 1.8)
            .Object
          })

setMethod("initialize", "ParHmmSnpSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$col <- c("lightblue", "red", "lightblue", "green3")
            .Object$bg <- c("lightblue", "red", "lightblue", "green3")
            .Object$ylab <- "copy number"
            ##add col.predict as an argument
            require(RColorBrewer, quietly=TRUE) || stop("RColorBrewer package not available")
            col.predict <- brewer.pal(7, "BrBG")[c(7, 2)]
            snpPar(.Object)$col.predict <- c(col.predict[1], "white", "yellow", col.predict[2])
            .Object
          })

setMethod("initialize", "ParHmmSnpCopyNumberSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$ylab <- "copy number"
            require(RColorBrewer, quietly=TRUE) || stop("RColorBrewer package not available")
            col.predict <- brewer.pal(7, "BrBG")[c(7, 2)]
            snpPar(.Object)$col.predict <- c(col.predict[1], "white", col.predict[2])               
            .Object
          })

setMethod("initialize", "ParHmmSnpCallSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
##            if(.Object$add.cytoband){
##              .Object$heights <- c(1, 1, 0.8)
##              } else {
##                .Object$heights <- c(1, 1)
##              }
            .Object$outer.ylab <- FALSE
            .Object$ylim <- c(0, 1)
            .Object$ylab <- "genotype call"
            .Object$col <- c("lightblue", "red", "lightblue")
            snpPar(.Object)$col.predict <- c("black", "white")
            .Object
          })

##One big matrix that could potentially span several chromosomes.
##featureData should be looked up in RSQLite database
setMethod("initialize", "HmmParameter",
          function(.Object,
                   snpset,
                   states=character(),
                   tau=matrix(),
                   tau.scale=matrix(),                   
                   pi=numeric(),
                   notes=character(),
                   beta=matrix(),
                   SCALE=2,##by default, twice as likely to return to normal state
                   ...){
            snpset <- snpset[!(is.na(chromosome(snpset))), ]
            chrom <- chromosome2numeric(chromosome(snpset))
            snpset <- snpset[order(chrom, position(snpset)), ]

            hmmOptions <- new("HmmOptions", snpset=snpset, states=states, beta=beta, ...)
            S <- length(hmmOptions@states)
            featureNames <- featureNames(snpset)
            featureData <- addFeatureData(snpset)

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
            if(length(beta) == 1){
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
          }

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





