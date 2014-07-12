getOligoset <- function(){
  path <- tryCatch(system.file("extdata", package="VanillaICE", mustWork=TRUE), error=function(e) NULL)
  if(is.null(path)) path <- "~/Software/bridge/VanillaICE/inst/extdata"
  load(file.path(path, "oligosetForUnitTest.rda"))
  oligoset
}

getSE <- function() {
  se <- as(getOligoset(), "SnpArrayExperiment")
  copyNumber(se) <- copyNumber(se)/100
  baf(se) <- baf(se)/1000
  se
}

test_HmmGRanges <- function(){
  library(oligoClasses)
  library(IRanges)
  se <- getSE()
  states <- Rle(factor(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3)),
                as.integer(c(996, 102, 902, 50, 2467, 102, 76, 1822,
                             99, 900, 20, 160)))
  hgr <- HmmGRanges(states=states,
                    feature_starts=start(se),
                    feature_chrom=chromosome(se),
                    loglik=15.3)
  checkTrue(validObject(hgr))
}
