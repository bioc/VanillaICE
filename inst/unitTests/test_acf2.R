test_acf2 <- function(){
  set.seed(123)
  x <- rnorm(100)
  lag.max <- 10
  answer <- acf(x, lag.max=lag.max, plot=FALSE)[[1]][lag.max+1, , 1]
  checkEquals(acf2(x), answer)
}
