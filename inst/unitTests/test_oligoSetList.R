test_oligoSetList <- function(){
  library(crlmm)
  foreach::registerDoSEQ()
  data(cnSetExample,package="crlmm")
  oligoList <- OligoSetList(cnSetExample)
  checkTrue(validObject(oligoList))
  checkTrue(validObject(oligoList[[1]]))
}
