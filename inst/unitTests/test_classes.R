test_Viterbi2 <- function(){
	checkTrue(validObject(new("Viterbi")))
	checkTrue(validObject(new("Viterbi", numberFeatures=100L,
				  numberStates=2L,
				  normalIndex=1L)))
	checkTrue(validObject(new("Viterbi2")))
	checkTrue(validObject(new("Viterbi2", numberFeatures=100L,
				  numberStates=6L)))
}
