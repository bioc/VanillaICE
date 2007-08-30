setClassUnion("NullOrNumeric", c("numeric", "NULL"))
setClassUnion("NullOrNumericOrMatrix", c("NULL", "numeric", "matrix"))
setClassUnion("NullOrMatrix", c("NULL", "matrix"))
setClass("hSet", contains="oligoSnpSet",
         representation(distance="logical",
                        initialStateProbability="NullOrNumeric",
                        transitionProbability="NullOrNumericOrMatrix",
                        emissionProbability="NullOrMatrix",
                        likelihood="NullOrMatrix",
                        chromosomeArm="character",
                        stateNames="character",
                        "VIRTUAL"),
         prototype=list(
           distance=TRUE,
           initialStateProbability=NULL,
           transitionProbability=NULL,
           emissionProbability=NULL,
           likelihood=NULL,
           chromosomeArm=character(),
           stateNames=character()))

setClass("HmmSnpCopyNumberSet", contains="hSet",
         representation(locationCopyNumber="numeric",
                        scaleCopyNumber="NullOrNumeric",
                        log2="logical",
                        predictions="matrix",
                        copyNumberIce="logical"),
         prototype=list(
           log2=logical(),
           initialStateProbability=c((1-0.99)/2, 0.99, (1-0.99)/2),
           locationCopyNumber=c(1, 2, 3),
           scaleCopyNumber=NULL,
           stateNames=c("deletion", "normal", "amplification"),
           copynumberIce=FALSE))


setClass("HmmSnpCallSet", contains="hSet", 
         representation(pCHOM="numeric",  #vector with length = number of states
                        predictions="matrix",
                        callsIce="logical"),
         prototype=list(
           initialStateProbability=c(1-0.99, 0.99),
           pCHOM=c(1-1e-3, 0.7),
           stateNames=c("loss", "retention"),
           callsIce=FALSE))

setClass("HmmSnpSet", contains=c("HmmSnpCopyNumberSet", "HmmSnpCallSet"),
         prototype=list(
           initialStateProbability=c((1-0.99)/3, 0.99, (1-0.99)/3, (1-0.99)/3),
           pCHOM=c(1-1e-3, 0.7, 1-1e-3, 0.7),
           stateNames=c("deletion", "normal", "LOH", "amplification"),
           locationCopyNumber=c(1, 2, 2, 3),
           copyNumberIce=FALSE,
           callsIce=FALSE))

##Classes for graphical parameters
setClass("ParHSet", contains="ParESet")
setClass("ParHmmSnpCallSet", contains="ParHSet")
setClass("ParHmmSnpCopyNumberSet", contains="ParHSet")
setClass("ParHmmSnpSet", contains="ParHSet")

setAs("oligoSnpSet", "HmmSnpSet", 
      function(from, to){
        object <- new("HmmSnpSet",
                      calls=calls(from),
                      callsConfidence=callsConfidence(from),
                      copyNumber=copyNumber(from),
                      cnConfidence=cnConfidence(from),
                      featureData=featureData(from),
                      phenoData=phenoData(from),
                      annotation=annotation(from),
                      experimentData=experimentData(from),
                      scaleCopyNumber=NULL,
                      predictions=matrix(NA, nrow=nrow(from), ncol=ncol(from)),
                      stateNames=c("deletion", "normal", "LOH", "amplification"),
                      copyNumberIce=FALSE,
                      callsIce=FALSE)
        object
      })


setAs("HmmSnpSet", "HmmSnpCopyNumberSet", 
      function(from, to){
        object <- new("HmmSnpCopyNumberSet",
                      copyNumber=copyNumber(from),
                      cnConfidence=cnConfidence(from),
                      featureData=featureData(from),
                      phenoData=phenoData(from),
                      annotation=annotation(from),
                      experimentData=experimentData(from),
                      locationCopyNumber=locationCopyNumber(from)[c(1, 2, 4)],
                      scaleCopyNumber=scaleCopyNumber(from)[c(1, 2, 4)],
                      predictions=matrix(NA, nrow=nrow(from), ncol=ncol(from)),
                      log2=from@log2,
                      copyNumberIce=copyNumberIce(from),
                      distance=distance(from),
                      stateNames=c("deletion", "normal", "amplification"))
        object
      })

setAs("HmmSnpSet", "HmmSnpCallSet", 
      function(from, to){
        object <- new("HmmSnpCallSet",
                      calls=calls(from),
                      callsConfidence=callsConfidence(from),
                      featureData=featureData(from),
                      phenoData=phenoData(from),
                      annotation=annotation(from),
                      experimentData=experimentData(from),
                      pCHOM=pCHOM(from)[1:2],
                      predictions=matrix(NA, nrow=nrow(from), ncol=ncol(from)),
                      callsIce=callsIce(from),
                      distance=distance(from),
                      stateNames=c("loss", "retention"))
        object
      })

setAs("oligoSnpSet", "HmmSnpCopyNumberSet", 
      function(from, to){
        object <- new("HmmSnpCopyNumberSet",
                      copyNumber=copyNumber(from),
                      cnConfidence=cnConfidence(from),
                      featureData=featureData(from),
                      phenoData=phenoData(from),
                      annotation=annotation(from),
                      experimentData=experimentData(from),
                      scaleCopyNumber=matrix(rep(NA, 3, ncol(from))),
                      predictions=matrix(NA, nrow=nrow(from), ncol=ncol(from)),
                      copyNumberIce=FALSE,
                      stateNames=c("deletion", "normal", "amplification"))
        object
      })

setAs("SnpCopyNumberSet", "HmmSnpCopyNumberSet", 
      function(from, to){
        object <- new("HmmSnpCopyNumberSet",
                      copyNumber=copyNumber(from),
                      cnConfidence=cnConfidence(from),
                      featureData=featureData(from),
                      phenoData=phenoData(from),
                      annotation=annotation(from),
                      experimentData=experimentData(from),
                      scaleCopyNumber=matrix(rep(NA, 3, ncol(from))),
                      predictions=matrix(NA, nrow=nrow(from), ncol=ncol(from)),
                      copyNumberIce=FALSE,
                      stateNames=c("deletion", "normal", "amplification"))
        object
      })

setAs("oligoSnpSet", "HmmSnpCallSet", 
      function(from, to){
        object <- new("HmmSnpCallSet",
                      calls=calls(from),
                      callsConfidence=callsConfidence(from),
                      featureData=featureData(from),
                      phenoData=phenoData(from),
                      annotation=annotation(from),
                      experimentData=experimentData(from),
                      predictions=matrix(NA, nrow=nrow(from), ncol=ncol(from)),
                      callsIce=FALSE,
                      stateNames=c("loss", "retention"))
        object
      })

setAs("SnpCallSet", "HmmSnpCallSet", 
      function(from, to){
        object <- new("HmmSnpCallSet",
                      calls=calls(from),
                      callsConfidence=callsConfidence(from),
                      featureData=featureData(from),
                      phenoData=phenoData(from),
                      annotation=annotation(from),
                      experimentData=experimentData(from),
                      predictions=matrix(NA, nrow=nrow(from), ncol=ncol(from)),
                      callsIce=FALSE,
                      stateNames=c("loss", "retention"))
        object
      })




                        
         
           
           
           
                        
                        
