test_oligoSetList <- function(){
  require(crlmm)
  foreach::registerDoSEQ()
  data(cnSetExample,package="crlmm")
  oligoList <- OligoSetList(cnSetExample)
  checkTrue(validObject(oligoList))
  checkTrue(validObject(oligoList[[1]]))
}
