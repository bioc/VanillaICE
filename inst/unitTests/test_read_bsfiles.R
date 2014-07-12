test_read_bsfiles <- function(){
  library(data.table)
  path <- tryCatch(system.file("extdata", package="VanillaICE", mustWork=TRUE), error=function(e) NULL)
  if(is.null(path)) path <- "~/Software/VanillaICE/inst/extdata"
  file <- list.files(path, pattern="genome_studio", full.names=TRUE)

  tmp <- fread(file, nrows=0)
  keep <- c("SNP Name", "Allele1 - AB", "Allele2 - AB",
            "Log R Ratio", "B Allele Freq")
  select <- which(names(tmp)%in%keep)
  dat <- fread(file, select=select)
  checkTrue(all(dim(dat)==c(19L,3L)))

  ##
  ## Get feature annotation
  ##
  files <- rep(file, 5)
  pkg <- "human610quadv1bCrlmm"
  checkTrue(validObject(BeadStudioSet(files, header_info,
                                      "hg18", pkg)))
  checkTrue(validObject(BeadStudioSet(files[1], header_info,
                                      "hg18", pkg)))

  fname <- file.path(path, "LRRandBAF.txt")
  keep <- c("Name", "sample1.Log.R.Ratio", "sample1.B.Allele.Freq")
  labels <- setNames(bead_studio_variables()[3:4], keep[-1])
  header_info <- headerInfo(fname, skip=0, sep="\t",
                            keep=keep, labels=labels,
                            classes=c("character", rep("numeric", 2)))
  bsSet <- BeadStudioSet(fname, header_info,
                         annotationPkg="genomewidesnp6Crlmm",
                         genome="hg19")
  checkTrue(validObject(bsSet))

  ##
  ## BeadStudioSetList
  ##
  blist <- BeadStudioSetList(bsSet)
  checkTrue(validObject(blist))
  checkTrue(validObject(blist[[1]]))

  blist2 <- BeadStudioSetList(fname, header_info,
                              annotationPkg="genomewidesnp6Crlmm",
                              genome="hg19")
  checkTrue(validObject(blist2))
  checkTrue(validObject(blist2[[1]]))
}
