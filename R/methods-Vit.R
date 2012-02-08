setMethod("nStates", signature(object="Vit"),
	  function(object){
		  object@numberStates
	  })

setMethod("initialProb", signature(object="Vit"),
	  function(object){
		  object@initialProb
	  })

setMethod("emission", signature(object="Vit"),
	  function(object){
		  object@emission
	  })

setMethod("viterbiStatePath", signature(object="Vit"),
	  function(object){
		  object@viterbiStatePath
	  })

setMethod("normalStateIndex", signature(object="Vit"),
	  function(object){
		  object@normalIndex
	  })

setMethod("normal2altered", signature(object="Vit"),
	  function(object){
		  object@normal2altered
	  })

setMethod("altered2normal", signature(object="Vit"),
	  function(object){
		  object@altered2normal
	  })

setMethod("altered2altered", signature(object="Vit"),
	  function(object){
		  object@altered2altered
	  })

setMethod("transitionProb", signature(object="Vit"),
	  function(object){
		  object@transitionProb
	  })

setMethod("nrow", signature(x="Vit"),
	  function(x){
		  x@numberFeatures
	  })
