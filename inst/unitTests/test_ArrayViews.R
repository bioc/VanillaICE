test_ArrayViews <- function(){
  require(BSgenome.Hsapiens.UCSC.hg18)
  require(data.table)
  require(foreach)
  require(oligoClasses)

  checkTrue(validObject(ArrayViews()))
  registerDoSEQ()
  extdir <- system.file("extdata", package="VanillaICE", mustWork=TRUE)

  features <- suppressWarnings(fread(file.path(extdir, "SNP_info.csv")))
  fgr <- GRanges(paste0("chr", features$Chr), IRanges(features$Position, width=1),
                 isSnp=features[["Intensity Only"]]==0)
  fgr <- SnpGRanges(fgr)
  names(fgr) <- features[["Name"]]
  seqlevels(fgr) <- seqlevels(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(BSgenome.Hsapiens.UCSC.hg18) %in% seqlevels(fgr)]
  seqinfo(fgr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(fgr),]
  fgr <- sort(fgr)
  checkTrue(validObject(fgr))

  files <- list.files(extdir, full.names=TRUE, recursive=TRUE, pattern="FinalReport")
  ids <- gsub(".rds", "", gsub("FinalReport", "", basename(files)))
  ## select a permanent location to store the parsed data

  dat <- fread(files[1], nrows=0)
  keep <- c("SNP Name", "Allele1 - AB", "Allele2 - AB", "Log R Ratio", "B Allele Freq")
  select <- which(names(dat)%in%keep)
  dat <- fread(files[1], select=select)
  index_genome <- match(names(fgr), dat[["SNP Name"]])

  scan_params <- CopyNumScanParams(index_genome=index_genome,
                                   select=as.integer(select))
  views <- ArrayViews(rowData=fgr,
                      sourcePaths=files,
                      sample_ids=ids,
                      parsedPath=tempdir())

  checkTrue(validObject(views))

  checkTrue(validObject(scan_params))

  checkTrue(identical(VanillaICE:::indexGenome(scan_params), index_genome[!is.na(index_genome)]))

  checkTrue(identical(VanillaICE:::selectCols(scan_params), select))

  checkTrue(identical(VanillaICE:::scale(scan_params), 1000))

  checkTrue(identical(VanillaICE:::cnvar(scan_params),  "Log R Ratio"))

  checkTrue(identical(VanillaICE:::bafvar(scan_params), "B Allele Freq"))

  checkTrue(identical(c("Allele1 - AB", "Allele2 - AB"), VanillaICE:::gtvar(scan_params)))

  if(FALSE){
    x <- unlist(x)
    checkTrue(is.null(x))

    r <- head(lrr(views))
    b <- head(baf(views))
    checkTrue(identical(dim(r), c(6L,6L)))
    checkTrue(identical(dim(b), c(6L,6L)))
    checkTrue(identical(rownames(b), rownames(views)[1:6]))

    ## Changing the colnames of the views object should not change the
    ## way that the parsed files are accessed (i.e., files are accessed
    ## by a name derived from the source files)
    colnames(views) <- letters[seq_len(ncol(views))]
    r2 <- head(lrr(views))
    colnames(r2) <- colnames(r) <- NULL
    checkIdentical(r, r2)

    ## Fit a 6-state HMM
    se <- SnpExperiment(views)
    if(FALSE){
      snp_exp <- se
      save(snp_exp, file="~/Software/bridge/VanillaICE/data/snp_exp.rda")
    }
    emission_param <- EmissionParam(temper=1/2, p_outlier=1/100)
    ## to few markers for fitting the HMM
    fit <- hmm2(se, emission_param=emission_param)
    checkTrue(validObject(fit))
  }
}

test_columnSubset <- function(){
  require(BSgenome.Hsapiens.UCSC.hg18)
  require(data.table)
  require(foreach)
  require(oligoClasses)
  registerDoSEQ()
  ##extdir <- system.file("extdata", package="VanillaICE", mustWork=TRUE)
  extdir <- "~/Software/bridge/VanillaICE/inst/extdata"
  features <- suppressWarnings(fread(file.path(extdir, "SNP_info.csv")))
  fgr <- GRanges(paste0("chr", features$Chr), IRanges(features$Position, width=1),
                 isSnp=features[["Intensity Only"]]==0)
  fgr <- SnpGRanges(fgr)
  names(fgr) <- features[["Name"]]
  seqlevels(fgr) <- seqlevels(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(BSgenome.Hsapiens.UCSC.hg18) %in% seqlevels(fgr)]
  seqinfo(fgr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(fgr),]
  fgr <- sort(fgr)

  files <- list.files(extdir, full.names=TRUE, recursive=TRUE, pattern="FinalReport")

  dat <- fread(files[1], nrows=0)
  keep <- c("SNP Name", "Allele1 - AB", "Allele2 - AB", "Log R Ratio", "B Allele Freq")
  select <- which(names(dat)%in%keep)
  dat <- fread(files[1], select=select)
  index_genome <- match(names(fgr), dat[["SNP Name"]])
  scan_params <- CopyNumScanParams(index_genome=index_genome,
                                   select=as.integer(select))
  parse_path <- tempdir()
  views <- ArrayViews(rowData=fgr,
                      sourcePaths=files,
                      parsedPath=parse_path)
  parseSourceFile(views, scan_params)

  sample_info <- read.csv(file.path(extdir, "sample_data.csv"), stringsAsFactors=FALSE)
  ind_id <- setNames(gsub(" ", "", sample_info$IndividualID), sample_info$File)
  colnames(views) <- ind_id[gsub(".txt", "", colnames(views))]
  views2 <- views[, c("22169_03", "22169_02", "22169_01")]
  r1 <- lrr(views)[, "22169_01"]
  r2 <- lrr(views2)[, "22169_01"]
  checkTrue(identical(r1, r2))
}
