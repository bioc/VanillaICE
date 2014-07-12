test_rowMAD <- function(){
  x <- matrix(rnorm(100), 10, 10)
  result1 <- apply(x, 1, mad)
  result2 <- rowMAD(x)
  checkIdentical(result1, result2)
}
