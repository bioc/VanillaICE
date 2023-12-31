simulateData <- function(answer,
                         r_means=c(-2, -0.5, 0, 0, 0.5, 1),
                         b_sd=0.05, r_sd=0.2){
  cn <- as.integer(answer)
  r <- rnorm(length(cn), mean=r_means[cn], r_sd)
  g <- sample(1:3, length(cn), replace=TRUE)
  g[cn==4] <- sample(c(1,3), sum(cn==4), replace=TRUE)
  bmeans <- baf_means(EmissionParam())
  b <- rep(NA, length(cn))
  rtrnorm <- VanillaICE:::rtrnorm
  b[g==1] <- rtrnorm(sum(g==1), 0, b_sd)
  b[g==3] <- rtrnorm(sum(g==3), 1, b_sd)
  b[g==2 & cn==3] <- rnorm(sum(g==2&cn==3), 0.5, b_sd)
  b[cn==5] <- rnorm(sum(cn==5), 1/3, b_sd)
  b[cn==2] <- rtrnorm(sum(cn==2), 0, b_sd)
  b[cn==1] <- runif(sum(cn==1))
  list(r, b, g)
}

getOligoset <- function(){
  path <- tryCatch(system.file("extdata", package="VanillaICE", mustWork=TRUE), error=function(e) NULL)
  if(is.null(path)) path <- "~/Software/bridge/VanillaICE/inst/extdata"
  load(file.path(path, "oligosetForUnitTest.rda"))
  oligoset
}

getSE <- function() {
  load(system.file("extdata", "snp_exp.rda", package="VanillaICE"))
  snp_exp
}

test_SnpArrayExperiment <- function(){
  library(VanillaICE)
  library(oligoClasses)
  checkTrue(validObject(VanillaICE:::SnpDataFrame()))
  checkTrue(validObject(VanillaICE:::SnpDataFrame(DataFrame())))
  sdf <- VanillaICE:::SnpDataFrame()
  checkException(VanillaICE:::SnpDataFrame(x=3))
  checkException(VanillaICE:::SnpDataFrame(DataFrame(x=3)))
  checkTrue(validObject(VanillaICE:::SnpDataFrame(DataFrame(x=3, isSnp=FALSE))))
  df <- DataFrame(x=3, isSnp=TRUE)
  checkTrue(validObject(sdf <- as(df, "SnpDataFrame")))
  sdf2 <- VanillaICE:::SnpDataFrame(DataFrame(x=3, isSnp=TRUE))
  checkIdentical(sdf2, sdf)
  df <- VanillaICE:::SnpDataFrame(DataFrame(x=3, isSnp=FALSE))
  checkIdentical(isSnp(df), FALSE)

  checkTrue(validObject(SnpGRanges()))
  checkTrue(validObject(SnpGRanges(GRanges())))
  x <- GRanges(Rle("chr3", 5), IRanges(1:5, width=1))
  checkException(SnpGRanges(x))
  checkTrue(validObject(SnpGRanges(x, isSnp=logical(5))))
  x$isSnp <- FALSE
  checkTrue(validObject(SnpGRanges(x)))

  y <- x <- matrix(10, 5, 2)
  colnames(x) <- colnames(y) <- letters[1:2]
  rowranges <- GRanges(Rle("chr1", 5), IRanges(1:5, width=1))
  checkException(validObject(SnpArrayExperiment(x, y)))
  checkException(validObject(SnpArrayExperiment(x, y, rowRanges=rowranges)))
  checkTrue(validObject(SnpArrayExperiment(x, y, rowRanges=rowranges,
                                           isSnp=rep(FALSE,5))))
  rowranges$isSnp <- rep(FALSE,5)
  checkTrue(validObject(SnpArrayExperiment(x, y, rowRanges=rowranges)))
  se <- SnpArrayExperiment(x, y, rowranges)
  checkIdentical(isSnp(se), rep(FALSE, 5))

  se <- getSE()
  answer <- Rle(as.integer(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3)),
                         as.integer(c(996, 102, 902, 50, 2467, 102, 76, 1822,
                                      99, 900, 20, 160)))
  r_means <- c(-2, -0.5, 0, 0, 0.5, 1)
  b_sd <- 0.05
  r_sd <- 0.2
  datlist <- simulateData(answer, r_means, b_sd, r_sd)
  r <- as.matrix(datlist[[1]])
  b <- as.matrix(datlist[[2]])
  g <- as.matrix(datlist[[3]])
  colnames(r) <- colnames(b) <- colnames(g) <- "a"
  validObject(SnpArrayExperiment(cn=r, baf=b,
                                 rowRanges=rowRanges(se)))
}


test_SnpArrayExperiment2 <- function(){
  library(BSgenome.Hsapiens.UCSC.hg18)
  sl <- seqlevels(BSgenome.Hsapiens.UCSC.hg18)
  require(data.table)
  extdir <- system.file("extdata", package="VanillaICE", mustWork=TRUE)
  features <- suppressWarnings(fread(file.path(extdir, "SNP_info.csv")))
  fgr <- GRanges(paste0("chr", features$Chr),
                 IRanges(features$Position, width=1),
                 isSnp=features[["Intensity Only"]]==0)
  fgr <- SnpGRanges(fgr)
  names(fgr) <- features[["Name"]]
  seqlevels(fgr) <- sl[sl %in% seqlevels(fgr)]
  seqinfo(fgr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(fgr),]
  fgr <- sort(fgr)
  files <- list.files(extdir, full.names=TRUE,
                      recursive=TRUE,
                      pattern="FinalReport")
  parsedDir <- tempdir()
##  views <- ArrayViews(rowRanges=fgr,
##                      sourcePaths=files,
##                      parsedPath=parsedDir)
##  dat <- fread(files[1], skip="[Data]")
##  select_columns <- match(c("SNP Name", "Allele1 - AB", "Allele2 - AB",
##                            "Log R Ratio", "B Allele Freq"), names(dat))
##  index_genome <- match(names(fgr), dat[["SNP Name"]])
##  scan_params <- CopyNumScanParams(index_genome=index_genome,
##                                 select=select_columns,
##                                 cnvar="Log R Ratio",
##                                 bafvar="B Allele Freq",
##                                 gtvar=c("Allele1 - AB", "Allele2 - AB"))
##  parseSourceFile(views, scan_params)
##  views2 <- views[, 4:5]
##  checkTrue(validObject(views2))
##  rr <- rowRanges(views2)
##  obj <- SnpGRanges(rr,
##                    isSnp=rep(TRUE, nrow(view2)))
##  checkTrue(validObject(obj))
##  checkTrue(validObject(SnpExperiment(views2)))
}
