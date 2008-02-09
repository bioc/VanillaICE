.getCallEmission <- function(x,  ##genotype estimates from oligo
                             confidence, ##confidence scores from oligo
                             hapmapP,  ##obtained from hapmap
                             options){  ##HmmOptions object
  ##map observed call probablilities to a nonparametric estimate of
  ##the density using the reference p distributions when the truth is
  ##known
   gte <- x; rm(x)
  i11 <- hapmapP[, 1] == 3  ##called homozygous truth homozygous
  i12 <- hapmapP[, 1] == 4  ##called homozygous truth heterozygous
  i21 <- hapmapP[, 1] == 1  ##called het truth hom
  i22 <- hapmapP[, 1] == 2  ##called het truth het
  f11 <- density(hapmapP[i11, 2], from=1/1e3, to=1, n=1e3)
  f12 <- density(hapmapP[i12, 2], from=1/1e3, to=1, n=1e3)
  f21 <- density(hapmapP[i21, 2], from=1/1e3, to=1, n=1e3)
  f22 <- density(hapmapP[i22, 2], from=1/1e3, to=1, n=1e3)
  
  ###########################################################################
  ##distribution of observed call probabilities when the call was homozygous
  observedPcalledHom <- cut(as.vector(confidence[gte == 1]), breaks=f11$x, labels=FALSE) #called HOM, truth HOM w
  ##marginal distribution of observed call probabilities for heterozygous calls  
  observedPcalledHet <- cut(as.vector(confidence[gte == 2]), breaks=f21$x, labels=FALSE) #called HET, truth HOM

  ###########################################################################    
  ##P(phat | LOH, gte=HET)  -- wrong call
  ###########################################################################
  ##this is approximately P(phat | HOM, gte=HET)  
  pTruthIsLoh <- rep(NA, length(gte))  
  pTruthIsLoh[gte == 1] <- f11$y[observedPcalledHom]

  ###########################################################################    
  ##Calculate P(phat | LOH, gte=HOM) 
  ###########################################################################
  ##this is approximately P(phat | HOM, gte=HOM)  
  pTruthIsLoh[gte == 2] <- f21$y[observedPcalledHet]

  ###########################################################################
  ##Calculate P(phat | Normal, CHET)
  ###########################################################################  
  ##This is approx. P(phat | HET, CHET) * P(HET | CHET) + P(phat | HOM, CHET) * P(HOM|CHET)
  ##When the true state is normal, we need to estimate 4 things conditional on the call
  ## P(p | CHET, normal):
  ##chet1  P(p | HET, CHET)
  ##chet2  P(p | HOM, CHET)
  ##chet3  P(HET | CHET)
  ##chet4  P(HOM | CHET)
  pTruthIsNormal <- rep(NA, length(gte))
  chet1 <- f22$y[cut(confidence[gte == 2], breaks=f22$x, labels=FALSE)]
  chet2 <- f21$y[cut(confidence[gte == 2], breaks=f21$x, labels=FALSE)]
  gt.gte <- options@gt.gte
  pTruthIsNormal[gte == 2] <- chet1*(1-gt.gte[2]) + chet2*gt.gte[2]

  ###########################################################################
  ##Calculate P(phat | Normal, HOM)
  ###########################################################################    
  ##In a normal region truth can be HOM or HET
  ##P(phat | Normal, CHOM) is approx. P(phat | HET, CHOM) * P(HET|CHOM) + P(phat | HOM, CHOM) * P(HOM | CHOM)
  chom1 <- f12$y[cut(confidence[gte == 1], breaks=f12$x, labels=FALSE)]
  ##chom3 <- 1-0.9999  ##P(HET|CHOM)
  chom2 <- f11$y[cut(confidence[gte == 1], breaks=f11$x, labels=FALSE)]
  ##chom4 <- 0.9999    ##P(HOM|CHOM)
  ##probability that the tue state is normal when genotype call is homozygous
  pTruthIsNormal[gte == 1] <- chom1*gt.gte[1] + chom2*(1-gt.gte[1])

  ###########################################################################
  ##the stateProbability matrix is R x 2 (R = number of rows)
  stateProbability <- cbind(pTruthIsLoh, pTruthIsNormal)
  colnames(stateProbability) <- c("L", "N")
  rownames(stateProbability) <- c("P[p(CHOM)| q]", "P[p(CHET)|q]")[gte]
  ##stateProbability is R x 2 (R rows and 2 columns)
  gt.confidence.states <- options@gt.confidence.states
  stateProbability <- stateProbability[, gt.confidence.states]
  
  ###########################################################################
  ##Calculate P(CHET| q) and P(CHOM| q), where q is
  ##the true state.
  ###########################################################################  
  ##For q = LOH, we make the following simplifying
  ##assumptions:
  ##1.  P(CHET | q=LOH) \approx P(CHET | HOM) 
  ##2.  P(CHOM | q=LOH) \approx P(CHOM | HOM)
  ##For q = Normal,
  ##3.  P(CHET | q=Normal) \approx P(HET | q = Normal)
  ##    estimate the latter from the data, or use a fixed number.  e.g., 0.3
  ##4.  P(CHOM | q=Normal) \approx P(HOM | q = Normal)
  ##    estimate the latter quantity from the data, or use a fixed number 0.7
  ##P.CHET.LOH <- 1-P.CHOM.LOH
  ##P.CHET.LOH <- 1-gte.state["L"]
  ##P.CHET.Normal <- 1-P.CHOM.Normal
  ##P.CHET.Normal <- 1-gte.state["N"]
  gte.state <- options@gte.state
  B <- cbind(gte.state, 1-gte.state, NA)
  colnames(B) <- c("gte=HOM", "gte=HET", "gte=NA")
  rownames(B) <- options@states
  B <- t(B[, gte])
  emission.gt <- B * stateProbability
  rownames(emission.gt) <- names(gte)
  emission.gt
}

addFeatureData <- function(snpset){
  featureNames <- featureNames(snpset)
  data(chromosomeAnnotation, package="SNPchip", envir=environment())
  chrAnn <- as.matrix(chromosomeAnnotation)
  tmp <- position(snpset) <= chromosomeAnnotation[chromosome(snpset), "centromereStart"]
  chromosomeArm <- as.character(tmp)
  chromosomeArm[chromosomeArm == "TRUE"] <- "p"
  chromosomeArm[chromosomeArm == "FALSE"] <- "q"

  fD <- data.frame(cbind(chromosome(snpset), arm=chromosomeArm), row.names=featureNames,
                   stringsAsFactors=FALSE)
  fD$position <- position(snpset) 
  colnames(fD) <- c("chromosome", "arm", "position")
  varMetadata <- data.frame(labelDescription=c("chromosome", "chromosomal arm", "physical position"), row.names=colnames(fD))
  featureData <- new("AnnotatedDataFrame", data=fD, varMetadata=varMetadata)
  featureData
}

##.calculateYlim <- function(object, op){
##  if("copyNumber" %in% ls(assayData(object))){
##    ##only print this if there is more than one sample or more than 1 chromosome to plot
##    if(length(unique(chromosome(object))) > 1 || ncol(object) > 1){
##      print("one.ylim is FALSE. Calculating ylim based on the percentiles of the copy number distribution")
##    }
##    if(op$log == "y"){
##      ylim <- range(copyNumber(object), na.rm=TRUE)
##    } else{
##      ##use 1 and 99 quantiles to avoid outliers
##      ylim <- c(quantile(copyNumber(object), prob=0.001, na.rm=TRUE),
##                quantile(copyNumber(object), prob=0.999, na.rm=TRUE))
##    }
##  } else{
##    ##ylimits for genotypes??
##    y <- .getY(object)
##    ##jitter the genotype calls
##    y <- jitter(y, amount=0.05)
##    ylim <- range(y)
##  }
##  ylim
##}




##For each of the breaks in xBreak, find the cytoband(s)
.locateCytoband <- function(breaks, cytoband){
  cytoband <- cytoband[cytoband[, "chrom"] == breaks["chromosome"], ]
  ##include cytobands that begin before the break and end in the break or after the break
  bands <- cytoband[cytoband[, "chromEnd"] > as.numeric(breaks["start"]) & cytoband[, "chromStart"] < as.numeric(breaks["last"]), "name"]
  bands <- paste(bands, collapse=", ")
}

##setMethod("calculateBreakpoints", "HmmParameter",
##          function(object, x, position, chromosome, states, digits=2, sampleNames, ...){
calculateBreakpoints <- function(x, position, chromosome, states, digits=2, sampleNames){
            if(any(is.na(x))){
              Nmissing <- sum(is.na(x))
            } else Nmissing <- 0
            N <- length(x)
            if(length(unique(x)) == 1) return(NULL)
            d <- diff(x)
            if(sum(abs(d), na.rm=TRUE) < 1){
              print("No breakpoints in this region")
              return(NULL)
            }
            index <- c(1, (2:N)[d != 0 & !is.na(d)])
            index <- cbind(index[-length(index)], index[2:length(index)])
            index[, 2] <- index[, 2] - 1
            lastBreak <- index[nrow(index), 2]
            index <- rbind(index, c(lastBreak+1, N))  

            ##index is a matrix of indices in object where breaks occured.
            ##Replace by physical position.
            physical.positions <- matrix(position[as.vector(index)], nrow(index), ncol(index))
            physical.positions <- physical.positions/1e6
            colnames(physical.positions) <- c("start", "last")
            size <- (physical.positions[, "last"] - physical.positions[, "start"])
            N <- index[, 2] - index[, 1] + 1
            physical.positions <- cbind(physical.positions, size, N)
            colnames(physical.positions)[3:4] <- c("size", "N")
            predictedState <- function(i, x, chromosome){
              if(any(is.na(i))){
                stop("missing values")
              }
              pred <- x[(i[1]:i[2])]
              pred <- pred[!is.na(pred)]
              pred <- unique(pred)
              if(length(pred) > 1) {
                stop("predictions not unique")
              }
              chrom <- chromosome[i[1]:i[2]]
              chrom <- paste(unique(chrom), collapse=", ")
              c(pred, chrom)
            }
            ##The vector of predicted states (integer)
            ps <- t(apply(index, 1, predictedState, x, chromosome))
            states <- states[as.numeric(ps[, 1])]
            breaks <- data.frame(physical.positions)
            breaks$state <- states
            breaks$chrom <- ps[, 2]
            ##breaks$id <- rep(sampleNames(object), nrow(breaks))

            ##################################################
            ##Find positions of adjacent SNPs
            ##################################################
            adjacentSnps <- function(x, position, chromosome, digits){
              start <- as.numeric(x["start"])*1e6
              last <- as.numeric(x["last"])*1e6
              chrom <- x["chrom"]
              if(any(position < start & chromosome == chrom)){
                prevSnp <- max(position[position < start & chromosome == chrom])
              } else prevSnp <- NA
              if(any(position > last & chromosome == chrom)){                 
                nextSnp <- min(position[position > last & chromosome == chrom])
              } else nextSnp <- NA
              adj <- c(prevSnp, nextSnp)/1e6
              adj
            }
            adjacent <- t(apply(breaks, 1, adjacentSnps, position=position,
                                chromosome=chromosome, digits=digits))
            colnames(adjacent) <- c("prev", "next")
            breaks <- cbind(breaks, adjacent)
            breaks$id <- rep(sampleNames, nrow(breaks))
            colnames(breaks) <- c("start", "last", "size", "N", "state", "chr", "prev", "next", "id")

            ##Add columns for the name of the first SNP and the name of the last SNP
            colnames <- c("id", "chr", "state", "size", "N", "start", "last", "prev", "next")
            breaks <- breaks[, colnames]
            breaks
          }

calculateCnSE <- function(object,
                          referenceSet,
                          epsilon=0.1){
  if(missing(referenceSet)) stop("Require a referenceSet for estimating SNP-specific across sample variation in copy number estimates")
  require(genefilter) || stop("genefilter not available")
  is.autosome <- chromosome(object) %in% as.character(1:22) 
  object <- object[is.autosome, ]  
  referenceSet <- referenceSet[match(featureNames(object), featureNames(referenceSet)), ]
  robustSD <- function(X){
    diff(quantile(X, probs=c(0.16, (1-0.16)), na.rm=TRUE))/2
  }
  within.sd <- apply(copyNumber(object), 2, robustSD)
  across.sd <- rowSds(copyNumber(referenceSet), na.rm=TRUE)
  across.sd <- matrix(across.sd, nrow=nrow(object), ncol=ncol(object), byrow=FALSE)
  ##scale across.sd by the median sd of the sample
  median.across.sd <- median(across.sd)
  std.across.sd <- across.sd/median.across.sd
  SE <- within.sd*std.across.sd
  SE[SE == 0] <- epsilon
  SE
}

##Code this in C
viterbi <- function(beta, pi, tau, states, arm){
  S <- length(states)
  T <- nrow(tau) + 1
  beta <- matrix(beta, nr=T, nc=S)
  delta <- psi <- matrix(NA, T, S)
  colnames(delta) <- colnames(psi) <- states
  delta[1, ] <- pi + beta[1, ]
  psi[1, ] <- rep(0, S)

  for(t in 2:T){
    if(arm[t] != arm[t-1]){
      ##SNP t is on a different chromosome/chromosome arm
      delta[t, ] <- pi + beta[t, ]
      psi[t, ] <- rep(0, S)
    }
    AA <- matrix(tau[t-1, ], nr=S, nc=S)
    for(j in 1:S){
      ##Equation 105b
      delta[t, j] <- max(delta[t-1, ] + AA[, j]) + beta[t, j]
      psi[t, j] <- order(delta[t-1, ] + AA[, j], decreasing=TRUE)[1]
    }
  }
  Pstar <- max(delta[nrow(delta), ])
  qhat <- rep(NA, nrow(delta))
  T <- nrow(delta)
  names(qhat) <- rownames(delta)  
  qhat[T] <- order(delta[T, ], decreasing=TRUE)[1]
  for(t in (T-1):1) qhat[t] <- psi[t+1, qhat[t+1]]
  qhat    
}



