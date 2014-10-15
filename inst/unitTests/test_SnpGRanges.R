test_SnpGRanges <- function(){
  g <- GRanges("chr1", IRanges(15L, 15L))
  sg <- SnpGRanges(g, isSnp=TRUE)
  checkTrue(validObject(sg))
  checkTrue(validObject(SnpGRanges()))
}
