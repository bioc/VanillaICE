setMethod("initialize", signature(.Object="Viterbi"),
	  function(.Object,
		   numberFeatures=0L,
		   numberStates=6L,
		   emission=matrix(0L, numberFeatures*numberStates, 1L),
		   initialProb=log(rep(1/numberStates, numberStates)),
		   arm=integer(numberFeatures),
		   transitionProb=matrix(numeric(max(c(numberFeatures-1, 0))), ncol=1L),
		   viterbiStatePath=integer(numberFeatures),
		   normalIndex=3L,
		   delta=matrix(0, nrow(emission), 1L),
		   normal2altered=1,
		   altered2normal=1,
		   altered2altered=1,
		   pAA=numeric(numberStates^2)){
		  callNextMethod(.Object,
				 numberFeatures=numberFeatures,
				 numberStates=numberStates,
				 emission=emission,
				 initialProb=initialProb,
				 arm=arm,
				 transitionProb=transitionProb,
				 viterbiStatePath=viterbiStatePath,
				 normalIndex=normalIndex,
				 delta=delta,
				 normal2altered=normal2altered,
				 altered2normal=altered2normal,
				 altered2altered=altered2altered,
				 pAA=pAA)
	  })

setMethod("delta", signature(object="Viterbi"),
	  function(object) object@delta)

setMethod("pAA", signature(object="Viterbi"),
	  function(object) object@pAA)

setValidity("Viterbi", function(object){
	ns <- nStates(object)
	##nf <- nrow(object)
	nf <- object@numberFeatures
	nf.ns <- ns*nf
	emit <- emission(object)
	if(ncol(emit) != 1){
		return("emission matrix should have a single column")
	}
	if(nrow(emit) != nf.ns){
		return("emission matrix should have number of features x number of states elements")
	}
	delta <- delta(object)
	if(!all(dim(emit) == dim(delta))){
		return("emission matrix and delta should have the same dimension")
	}
	statePath <- viterbiStatePath(object)
	if(length(statePath) != nf){
		return("viterbiStatePath should have the same length as the number of features")
	}
	if(ns > 0){
		ni <- normalIndex(object)
		if(ni <= 0 | ni > ns){
			return("the index for the normal state (normalIndex) \n must be an integer greater than zero and \n less than or equal to the number of states")
		}
		if(sum(exp(initialProb(object))) != 1)
			return("exp(initialProb) must sum to 1")
	}
	if(nf > 0){
		tp <- transitionProb(object)
		if(length(tp) != (nf-1)){
			return("transitionProb vector should be the number of features -1")
		}
	}
})

setMethod("show", signature(object="Viterbi"), function(object){
	cat("Number of features:", object@numberFeatures, "\n")
	cat("Number of states:", nStates(object), "\n")
})

setMethod("fit", signature(object="Viterbi"), function(object){
	fitViterbi(object)
})

fitViterbi <- function(object){
	tmp <- .C("viterbi",
		  emission(object),
		  initialProb(object),
		  transitionProb(object), ## log scale?
		  object@arm,
		  nStates(object),
		  object@numberFeatures,
		  viterbiStatePath(object),
		  delta(object),
		  object@normal2altered,
		  object@altered2normal,
		  object@altered2altered,
		  normalIndex(object),
		  pAA(object))
	new("Viterbi",
	    emission=emission(object),
	    initialProb=initialProb(object),
	    transitionProb=transitionProb(object),
	    arm=object@arm,
	    numberStates=nStates(object),
	    numberFeatures=object@numberFeatures,
	    viterbiStatePath=tmp[[7]],
	    delta=tmp[[8]],
	    normal2altered=object@normal2altered,
	    altered2normal=object@altered2normal,
	    altered2altered=object@altered2altered,
	    normalIndex=normalIndex(object),
	    pAA=pAA(object))
}
