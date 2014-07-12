.notneedtest_SummarizedExperimentAccessors <- function(){
  library(VanillaICE)
  library(oligoClasses)
  library(GenomicRanges)
  ##
  ## only test if local
  ## Could probably get rid of this
  if(path.expand("~") == "/Users/rscharpf"){
    library(GenEpiData2014)
    ##load_all("~/Teaching/GenEpiData2014")
    data(cleft_exp, package="GenEpiData2014")
    checkIdentical(lrr(cleft_exp), assays(cleft_exp)[["cn"]])
    checkIdentical(baf(cleft_exp), assays(cleft_exp)[["baf"]])
    checkIdentical(sum(isSnp(cleft_exp)), 9041L)
    checkIdentical(start(cleft_exp), position(cleft_exp))
    checkIdentical(genomeBuild(cleft_exp), "hg18")
  }
}
