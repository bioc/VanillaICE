setClassUnion("NULLorHmmPredict", c("NULL", "HmmPredict"))
setClass("HmmParameter",
	 representation(states="character",
			initialStateProbability="numeric",
			emission="array",
			genomicDistance="numeric",
			transitionScale="matrix"))
setClassUnion("NullOrSnpSet", c("NULL", "SnpLevelSet"))
setClassUnion("NullOrNumeric", c("NULL", "numeric"))
setClass("HmmOptions", ##contains="SnpLevelSet",
	 representation(states="character",
			snpset="SnpLevelSet",
			copyNumber.location="numeric",
			copyNumber.scale="NullOrNumeric",
			copyNumber.ICE="logical",
			calls.ICE="logical",
			probHomCall="numeric",
			term5="numeric"))







                        
         
           
           
           
                        
                        
