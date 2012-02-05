setClass("HmmOptionList", contains="list")

setClass("Vit",
	 representation(emission="matrix",
			initialProb="numeric",
			transitionProb="matrix",
			arm="integer",
			numberStates="integer",
			numberFeatures="integer",
			viterbiStatePath="integer",
			normal2altered="numeric",
			altered2normal="numeric",
			altered2altered="numeric",
			normalIndex="integer",
			"VIRTUAL"))

setClass("Viterbi", contains="Vit",
	 representation(delta="matrix",
			pAA="numeric"))

setClass("Viterbi2", contains="Vit",
	 representation(forwardVariable="matrix",
			backwardVariable="matrix",
			scaleFactor="numeric"))






