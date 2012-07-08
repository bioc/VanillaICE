test_bafconstraints <- function(){
	isNonDecreasing <- VanillaICE:::isNonDecreasing
	mus <- c(0, 1/4, 1/3, 1/2, 2/3, 3/4, 1)
	names(mus) <- c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B")
	checkTrue(isNonDecreasing(mus))

	makeMusBafNondecreasing <- VanillaICE:::makeMusBafNondecreasing
	mus[["AAAB"]] <- 0.35
	##trace(makeMusBafNondecreasing, browser)
	mus2 <- makeMusBafNondecreasing(mus)
	checkIdentical(mus2[["AAAB"]], 1/3)
	mus[["ABBB"]] <- 1/4
	mus2 <- makeMusBafNondecreasing(mus)
}
