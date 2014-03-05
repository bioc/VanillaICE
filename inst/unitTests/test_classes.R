test_BeadStudioSetList <- function(){
  library(oligoClasses);library(VanillaICE);library(RUnit)
  ##validObject(new("BeadStudioSetList"))
  ##checkEquals(length(new("BeadStudioSetList")), 0L)
  path <- system.file("extdata", package="VanillaICE")
  fname <- list.files(path, pattern="LRRand", full.names=TRUE)
  ##trace(BeadStudioSetList, browser)
  obj1 <- VanillaICE:::BeadStudioSetList(fnames=fname,
                                         annotationPkg="genomewidesnp6Crlmm",
                                         genome="hg19")
  checkTrue(validObject(obj1))
  bSet <- obj1[[1]]
  checkTrue(validObject(bSet))

  library(oligoClasses)
  library(ff)
  ldPath(tempdir())
  foreach:::registerDoSEQ()
  obj2 <- VanillaICE:::BeadStudioSetList(fnames=fname,
                                         annotationPkg="genomewidesnp6Crlmm",
                                         genome="hg19")
  checkTrue(validObject(obj2))
  checkTrue(identical(as.numeric(baf(obj1)[[1]]),
                      as.numeric(baf(obj2)[[1]][,1, drop=FALSE])))

}

test_BeadStudioSet <- function(){
  path <- system.file("extdata", package="VanillaICE")
  fname <- file.path(path, "LRRandBAF.txt")
  bsSet <- BeadStudioSet(fname, annotationPkg="genomewidesnp6Crlmm", genome="hg19")
  checkTrue(validObject(bsSet))
}

test_oligoSetList <- function(){
  library(crlmm)
  data(cnSetExample,package="crlmm")
  oligoList <- OligoSetList(cnSetExample)
  checkTrue(validObject(oligoList))
  checkTrue(validObject(oligoList[[1]]))
}
