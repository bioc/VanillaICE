test_FilterParam <- function(){
  fp <- FilterParam()
  checkTrue(length(probability(FilterParam()))==1)

  data(snp_exp)
  fit <- hmm2(snp_exp)
  checkTrue(is(unlist(fit), "HmmGRanges"))

  checkTrue(identical(length(fit), ncol(snp_exp)))
}
