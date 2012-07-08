##setGeneric("hmm", function(object, hmm.params, use.baf=FALSE, ...) standardGeneric("hmm"))
setGeneric("hmm", function(object, ...) standardGeneric("hmm"))
##setGeneric("hmm2", function(object, hmm.params, use.baf=FALSE, ...) standardGeneric("hmm2"))
setGeneric("cnEmission", function(object, stdev, ...) standardGeneric("cnEmission"))
setGeneric("gtEmission", function(object, hmm.params, gt.conf, is.snp, cdfName, ...) standardGeneric("gtEmission"))
setGeneric("bafEmission", function(object, ...) standardGeneric("bafEmission"))
setGeneric("sd", useAsDefault=function(x, na.rm=FALSE) stats::sd(x, na.rm))

setGeneric("nStates", function(object) standardGeneric("nStates"))
setGeneric("normalStateIndex", function(object) standardGeneric("normalStateIndex"))
setGeneric("emission", function(object) standardGeneric("emission"))
setGeneric("backwardVariable", function(object) standardGeneric("backwardVariable"))
setGeneric("forwardVariable", function(object) standardGeneric("forwardVariable"))
setGeneric("transitionProb", function(object) standardGeneric("transitionProb"))
setGeneric("viterbiStatePath", function(object) standardGeneric("viterbiStatePath"))
setGeneric("scaleFactor", function(object) standardGeneric("scaleFactor"))
setGeneric("fit", function(object) standardGeneric("fit"))
setGeneric("initialProb", function(object) standardGeneric("initialProb"))
setGeneric("normal2altered", function(object) standardGeneric("normal2altered"))
setGeneric("altered2normal", function(object) standardGeneric("altered2normal"))
setGeneric("altered2altered", function(object) standardGeneric("altered2altered"))
setGeneric("delta", function(object) standardGeneric("delta"))
setGeneric("pAA", function(object) standardGeneric("pAA"))

setGeneric("assayDataList", function(object) standardGeneric("assayDataList"))
setGeneric("featureDataList", function(object) standardGeneric("featureDataList"))



