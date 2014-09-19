test_FilterParam <- function(){
  fp <- FilterParam()
  checkTrue(length(probability(FilterParam()))==1)
}
