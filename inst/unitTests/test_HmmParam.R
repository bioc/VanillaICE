test_EmissionParam <- function(){
  checkTrue(validObject(EmissionParam()))
  checkTrue(validObject(TransitionParam()))
  checkTrue(validObject(HmmParam()))
  param <- HmmParam()
  emission(param) <- NULL
}
