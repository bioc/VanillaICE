.breaks <- function(object){
  ##  object <- object[!is.na(predictions(object)), ]
  if(sum(is.na(predictions(object))) > 0){
    warning("NA values in predictions")
  }
  if(length(unique(chromosome(object))) > 1) stop("should only be one chromosome in the call to this function")
  object <- object[order(position(object)), ]
  
  pred <- predictions(object)##[!is.na(predictions(object))]
  N <- nrow(object)
  
  if(length(unique(pred)) == 1) return(NULL)
  d <- diff(pred)
  index <- c(1, (2:N)[d != 0])
  index <- cbind(index[-length(index)], index[2:length(index)])
  index[, 2] <- index[, 2] - 1
  lastBreak <- index[nrow(index), 2]
  index <- rbind(index, c(lastBreak+1, nrow(object)))  

  ##index is a matrix of indices in object where breaks occured.
  ##Replace by physical position.
  physical.positions <- matrix(position(object)[as.vector(index)], nrow(index), ncol(index))
  colnames(physical.positions) <- c("start", "last")
  size <- (physical.positions[, "last"] - physical.positions[, "start"])/1e6
  physical.positions <- cbind(physical.positions, size)
  colnames(physical.positions)[3] <- "MB"

  predictedState <- function(x, object)  {
    pred <- predictions(object)[(x[1]:x[2])]
    pred <- pred[!is.na(pred)]
    pred <- unique(pred)
    if(length(pred) > 1) {
      print("predictions not unique")
      browser()
    }
    pred
  }
  ##The vector of predicted states (integer)
  ps <- apply(index, 1, predictedState, object)
  states <- stateNames(object)[ps]

  breaks <- data.frame(physical.positions)
  breaks$hiddenState <- states
  breaks$chromosome <- rep(unique(chromosome(object)), nrow(breaks))
  breaks$id <- rep(sampleNames(object), nrow(breaks))

  ##################################################
  ##Find positions of adjacent SNPs
  ##################################################
  adjacentSnps <- function(x, object){
    if(sum(position(object) < as.numeric(x["start"])) > 0){
      prevSnp <- max(position(object)[position(object) < as.numeric(x["start"])])
    } else prevSnp <- NA
    if(sum(position(object) > as.numeric(x["last"])) > 0){                 
      nextSnp <- min(position(object)[position(object) > as.numeric(x["last"])])
    } else nextSnp <- NA
    adj <- c(prevSnp, nextSnp)
    adj
  }
  adjacent <- t(apply(breaks, 1, adjacentSnps, object))
  colnames(adjacent) <- c("previous_SNP", "next_SNP")
  breaks <- cbind(breaks, adjacent)  
  
  colnames <- c("id", "chromosome", "previous_SNP", "start", "last", "next_SNP", "MB", "hiddenState")
  breaks <- breaks[, colnames]
  colnames(breaks) <- colnames
  breaks
}



.startingValues <- function(pi=NULL, A=NULL, B=NULL, cn=NULL, call=NULL, N=4){
  if(is.null(pi)) pi <- c(0.02, 0.91, 0.05, 0.02)
  if(is.null(A)){
    A <- matrix(c(0.98, 0.02, 0, 0,
                  0.02, 0.94, 0.02, 0.02,
                  0, 0.01, 0.99, 0,
                  0, 0.01, 0, 0.99), nc = 4, byrow=TRUE)
  }
  if(is.null(cn)){
    ##use empirical estimates from the summary output
    cn <- matrix(c(-1, 1,
                   0, 1,
                   0, 1,
                   1, 1), nc=2, byrow=TRUE)
  }
  if(is.null(call)){
    call <- matrix(c(0.01, 0.99,
                     0.3, 0.7,
                     0.01, 0.99,
                     0.3, 0.7), nc=2, byrow=TRUE)
  }
  B <- list(cn=cn, call=call)
  lambda <- list(N=4, pi=pi, A=A, B=B)
  lambda
}


##input: object of class AnnotatedSnpSet
##output: a 2 X G matrix of
##        probabilities for each sample.  2 is the number of hidden
##        states for genotype calls (normal or LOH).  G is the number of SNPs.
getCallEmission <- function(object,
                            CALL.STATES=2,
                            SAMPLE=1,
                            P.HOM.CHOM=0.999,  ##P(HOM | CHOM)
                            P.HET.CHOM=1-0.999,##P(HET | CHOM)
                            P.HET.CHET=0.999,  ##P(HET | CHET)
                            P.HOM.CHET=1-0.999,##P(HOM | CHET)
                            P.CHOM.Normal,        ##P(CHOM | HOM)
                            P.CHOM.LOH,
                            hapmapP){
  ##map observed call probablilities to a nonparametric estimate of
  ##the density using the reference p distributions when the truth is
  ##known
  calls(object)[calls(object)[, SAMPLE] == 3, SAMPLE] <- 1

##  hapmapP <- getHapmapProbabilities(object)
  i11 <- hapmapP[, 1] == 3  ##called homozygous truth homozygous
  i12 <- hapmapP[, 1] == 4  ##called homozygous truth heterozygous
  i21 <- hapmapP[, 1] == 1  ##called het truth hom
  i22 <- hapmapP[, 1] == 2  ##called het truth het
  f11 <- density(hapmapP[i11, 2], from=1/1e3, to=1, n=1e3)
  f12 <- density(hapmapP[i12, 2], from=1/1e3, to=1, n=1e3)
  f21 <- density(hapmapP[i21, 2], from=1/1e3, to=1, n=1e3)
  f22 <- density(hapmapP[i22, 2], from=1/1e3, to=1, n=1e3)

  if(ncol(object) > 1) warning("Only using first sample")

##  warning("callsConfidence in oligo version > 1.1.1 contains the LLR.  The AssayData element pAcc contains the probability of a correct calls, and this is what is used in the annotation objects (at least for the 500k)...probably need to update the 100k platform")  
  observedP <- callsConfidence(object)[, SAMPLE]
##  observedP <- assayDataElement(object, "pAcc")[, SAMPLE]

  observedCalls <- calls(object)[, SAMPLE]

  ###########################################################################
  ##distribution of observed call probabilities when the call was homozygous
  observedPcalledHom <- cut(as.vector(observedP[observedCalls == 1]), breaks=f11$x, labels=FALSE) #called HOM, truth HOM w
  ##marginal distribution of observed call probabilities for heterozygous calls  
  observedPcalledHet <- cut(as.vector(observedP[observedCalls == 2]), breaks=f21$x, labels=FALSE) #called HET, truth HOM

  ##The P(phat | HOM, CHET) is approximately P(phat | LOH, CHET
  ##AND
  ## P(phat | HOM, CHOM) is approximately P(phat | LOH, CHOM)
  pTruthIsLoh <- rep(NA, nrow(object))  ##matrix(NA, nrow(observedP), ncol=S) 
  pTruthIsLoh[observedCalls == 1] <- f11$y[observedPcalledHom]
  pTruthIsLoh[observedCalls == 2] <- f21$y[observedPcalledHet]


  ##P(phat | Normal, CHET) is approx. P(phat | HET, CHET) * P(HET | CHET) + P(phat | HOM, CHET) * P(HOM|CHET)
  ##When the true state is normal, we need to estimate 4 things conditional on the call
  ## P(p | CHET, normal):
  ##chet1  P(p | HET, CHET)
  ##chet2  P(p | HOM, CHET)
  ##chet3  P(HET | CHET)
  ##chet4  P(HOM | CHET)
  pTruthIsNormal <- rep(NA, nrow(object))
  chet1 <- f22$y[cut(observedP[observedCalls == 2], breaks=f22$x, labels=FALSE)]
##  chet3 <- 0.999
  chet2 <- f21$y[cut(observedP[observedCalls == 2], breaks=f21$x, labels=FALSE)]
##  chet4 <- 1-0.999
  pTruthIsNormal[observedCalls == 2] <- chet1*P.HET.CHET + chet2*P.HOM.CHET
  

  ##In a normal region truth can be HOM or HET
  ##P(phat | Normal, CHOM) is approx. P(phat | HET, CHOM) * P(HET|CHOM) + P(phat | HOM, CHOM) * P(HOM | CHOM)
  chom1 <- f12$y[cut(observedP[observedCalls == 1], breaks=f12$x, labels=FALSE)]
##  chom3 <- 1-0.9999  ##P(HET|CHOM)
  chom2 <- f11$y[cut(observedP[observedCalls == 1], breaks=f11$x, labels=FALSE)]
##  chom4 <- 0.9999    ##P(HOM|CHOM)
  ##probability that the tue state is normal when genotype call is homozygous
  pTruthIsNormal[observedCalls == 1] <- chom1*P.HET.CHOM + chom2*P.HOM.CHOM
  stateProbability <- cbind(pTruthIsLoh, pTruthIsNormal)
  colnames(stateProbability) <- c("q=LOH", "q=Normal")
  rownames(stateProbability) <- c("P[p(CHOM)| q]", "P[p(CHET)|q]")[observedCalls]

  ##Next we need to calculate P(CHET| q) and P(CHOM| q), where q is
  ##the true state.  For q = LOH, we make the following simplifying
  ##assumptions:
  ##1.  P(CHET | q=LOH) \approx P(CHET | HOM) 
  ##2.  P(CHOM | q=LOH) \approx P(CHOM| HOM)
  ##For q = Normal,
  ##3.  P(CHET | q=Normal) \approx P(HET | q = Normal)
  ##    estimate the latter from the data, or use a fixed number.  e.g., 0.3
  ##4.  P(CHOM | q=Normal) \approx P(HOM | q = Normal)
  ##    estimate the latter quantity from the data, or use a fixed number 0.7
  P.CHET.LOH <- 1-P.CHOM.LOH
  P.CHET.Normal <- 1-P.CHOM.Normal
  B <- matrix(c(P.CHOM.LOH, P.CHET.LOH,
                P.CHOM.Normal, P.CHET.Normal), ncol=2, byrow=TRUE)
  colnames(B) <- c("CHOM", "CHET")
  rownames(B) <- c("LOH", "Norm")
  B <- t(B[, observedCalls])
  stateProbability <- B * stateProbability
  stateProbability
}


hapmapProbabilities <- function(x){
##  require("callsConfidence") || stop("callsConfidence package not available")
  get(switch(annotation(x),
         pd.mapping50k.hind240=data(hindPhat),
         pd.mapping50k.xba240=data(xbaPhat),
         pd.mapping250k.nsp=data(nspPhat), ##need to calculate phats for 500k
         pd.mapping250k.sty=data(styPhat)))
}

getHapmapProbabilities <- function(object){
  hapmapP <- list()
  ##################################################
  ##50k platforms  
  if(annotation(object) == "mapping100k"){
    hind <- object[enzyme(object) == "Hind", ]
    annotation(hind) <- "pd.mapping50k.hind240"
    hapmapP[[1]] <- hapmapProbabilities(hind)

    xba <- object[enzyme(object) == "Xba", ]
    annotation(xba) <- "pd.mapping50k.hind240"
    hapmapP[[2]] <- hapmapProbabilities(xba)
    names(hapmapP) <- c("hind", "xba")
  }
  if(annotation(object) == "pd.mapping50k.hind240")
    hapmapP[[1]] <- hapmapProbabilities(object)
  if(annotation(object) == "pd.mapping50k.xba240")
    hapmapP[[1]] <- hapmapProbabilities(object)

  ##################################################
  ##250k platforms
  if(annotation(object) == "mapping500k"){
    nsp <- object[enzyme(object) == "Nsp", ]
    annotation(hind) <- "pd.mapping250k.nsp"
    hapmapP[[1]] <- hapmapProbabilities(nsp)

    xba <- object[enzyme(object) == "Sty", ]
    annotation(xba) <- "pd.mapping250k.sty"
    hapmapP[[2]] <- hapmapProbabilities(xba)
    names(hapmapP) <- c("nsp", "sty")
  }
  if(annotation(object) == "pd.mapping250k.nsp")
    hapmapP[[1]] <- hapmapProbabilities(object)
  if(annotation(object) == "pd.mapping250k.sty")
    hapmapP[[1]] <- hapmapProbabilities(object)
  
  hapmapP
}

callEmission <- function(object, P.CHOM.Normal, P.CHOM.LOH, SAMPLE=1){

  if(!(annotation(object) == "mapping50k.xba240" | annotation(object) == "mapping50k.hind240" | annotation(object) == "mapping250k.nsp" | annotation(object) == "mapping250k.sty" | annotation(object) == "mapping100k" | annotation(object) == "mapping500k")){
    stop("Confidence estimates for genotype calls are only supported for the 50k and 250k Affymetrix SNP chips.\n  The slot annotation must contain one of the following identifiers: mapping50k.xba240, mapping50k.hind240, mapping100k,  mapping250k.nsp, mapping250k.sty, or mapping500k")
  }
  ##hapmapP is a list.
  hapmapP <- getHapmapProbabilities(object)
  emissionP <- matrix(NA, nr=nrow(object), nc=2)
  if(annotation(object) == "mapping100k"){
    if(!any("enzyme" %in% fvarLabels(object))) stop("enzyme (Xba or Hind) must be stored in featureData")
    emissionP[enzyme(object) == "Xba", ] <- getCallEmission(object=object[enzyme(object) == "Xba", SAMPLE], hapmapP=hapmapP$xba,
                      P.CHOM.Normal=P.CHOM.Normal,
                      P.CHOM.LOH=P.CHOM.LOH)
    emissionP[enzyme(object) == "Hind", ] <- getCallEmission(object=object[enzyme(object) == "Hind", SAMPLE], hapmapP=hapmapP$hind,
                      P.CHOM.Normal=P.CHOM.Normal,
                      P.CHOM.LOH=P.CHOM.LOH)
  }
  if(annotation(object) == "mapping500k"){
    if(!any("enzyme" %in% fvarLabels(object))) stop("enzyme (Nsp or Sty) must be stored in featureData")    
    emissionP[enzyme(object) == "Nsp", ] <- getCallEmission(object=object[enzyme(object) == "Nsp", SAMPLE], hapmapP=hapmapP$nsp,
                      P.CHOM.Normal=P.CHOM.Normal,
                      P.CHOM.LOH=P.CHOM.LOH)
    emissionP[enzyme(object) == "Sty", ] <- getCallEmission(object=object[enzyme(object) == "Sty", SAMPLE], hapmapP=hapmapP$sty,
                      P.CHOM.Normal=P.CHOM.Normal,
                      P.CHOM.LOH=P.CHOM.LOH)
  } 
  if(annotation(object) == "mapping50k.xba240" | annotation(object) == "mapping50k.hind240" | annotation(object) == "mapping250k.nsp" | annotation(object) == "mapping250k.sty"){
    emissionP <- getCallEmission(object=object[, SAMPLE], hapmapP=hapmapP[[1]],
                                 P.CHOM.Normal=P.CHOM.Normal,
                                 P.CHOM.LOH=P.CHOM.LOH)
  } 
  emissionP
}


.drawRect <- function(x, object, col){
  if(length(x) < 1) return()
  start <- as.numeric(x["start"])
  last <- as.numeric(x["last"])
  predict <- predictions(object)[position(object) >= start & position(object) <= last & chromosome(object) == x["chromosome"]]
  predict <- predict[!is.na(predict)]
  if(length(unique(predict)) > 1) {
    warning("predictions not unique")
    browser()
  }
  col <- col[unique(predict)]
  rect(xleft=start,
       ybottom=0,
       xright=last,
       ytop=1,
       col=col,
       border=col)
}

plotPredictions <- function(object, op, breaks, ...){
  if(missing(breaks)) breaks <- getBreaks(object)
  if(missing(op)){
    op <- switch(class(object),
                 HmmSnpSet=new("ParHmmSnpSet", ...),
                 HmmSnpCallSet=new("ParHmmSnpCallSet", ...),
                 HmmSnpCopyNumberSet=new("ParHmmSnpCopyNumberSet", ...))    
    op$use.layout <- FALSE
    def.par <- par(...)
    on.exit(def.par)
  }
  chr <- unique(chromosome(object))[order(as.numeric(unique(chromosome(object))))]
  for(i in chr){
    if(op$use.chromosome.size) op$xlim <- c(0, chromosomeSize(i)) else op$xlim <- range(position(object[chromosome(object) == i, ]))
    for(j in sampleNames(object)){
      tmp <- breaks[breaks[, "chromosome"] == i & breaks[, "id"] == j, , drop=FALSE]
      plot(0:1,
           type="n",
           xlim=op$xlim,
           ylim=c(-0.2, 1.5),
           yaxt="n",
           xaxt="n",
           yaxs="i",
           xaxs=op$xaxs)
      if(i == op$firstChromosome){
        mtext("HMM", side=2, at=0.5, cex=op$cex.axis, las=1)
      }
      apply(tmp, 1, .drawRect, object=object[chromosome(object) == i, match(j, sampleNames(object))], col=op$col.predict)      
      rect(xleft=centromere(i)[1],
           ybottom=0,
           xright=centromere(i)[2],
           ytop=1,
           col="white",
           border="white")
      rect(xleft=op$xlim[1],
           ytop=0,
           ybottom=1,
           xright=op$xlim[2],
           border=op$border.prediction)
      if(op$legend.prediction) legend(op$legend.location.prediction, fill=op$col.predict,
                                      legend=stateNames(object),
                                      ncol=length(stateNames(object)), bty="n",
                                      cex=op$cex.legend)
    }
  }  
}

##For each of the breaks in xBreak, find the cytoband(s)
.locateCytoband <- function(breaks, cytoband){
  cytoband <- cytoband[cytoband[, "chrom"] == breaks["chromosome"], ]
  ##include cytobands that begin before the break and end in the break or after the break
  bands <- cytoband[cytoband[, "chromEnd"] > as.numeric(breaks["start"]) & cytoband[, "chromStart"] < as.numeric(breaks["last"]), "name"]
  bands <- paste(bands, collapse=", ")
}

