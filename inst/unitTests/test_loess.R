.test_Loess <- function(){
  library(VanillaICE)
  library(foreach)
  registerDoSEQ()
  extdir <- system.file("extdata", package="VanillaICE", mustWork=TRUE)
  files <- list.files(extdir, pattern="FinalReport")
  features <- suppressWarnings(fread(file.path(extdir, "SNP_info.csv")))
  fgr <- GRanges(paste0("chr", features$Chr), IRanges(features$Position, width=1),
                 isSnp=features[["Intensity Only"]]==0)
  fgr <- SnpGRanges(fgr)
  names(fgr) <- features[["Name"]]
  seqlevels(fgr) <- seqlevels(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(BSgenome.Hsapiens.UCSC.hg18) %in% seqlevels(fgr)]
  seqinfo(fgr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(fgr),]
  fgr <- sort(fgr)
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
  l_files <- gsub(".txt", "_lrr.rds", basename(files))
  b_files <- gsub(".txt", "_baf.rds", basename(files))
  g_files <- gsub(".txt", "_gt.rds", basename(files))
  views <- ArrayViews(rowData=fgr,
                      sourcePaths=files,
                      sample_ids=ids,
                      parsedPath=tempdir(),
                      lrrFiles=l_files,
                      bafFiles=b_files,
                      gtFiles=g_files)
  parseSourceFile(views, scan_params)

  wavedir <- file.path(parsedPath(views), "loess")
  dir.create(wavedir)
  viewsLoess <- views
  lrrFile(viewsLoess) <- file.path(wavedir, basename(lrrFile(views)))

  load_all("~/Labs/Klein")
  for(k in seq_len(ncol(views))){
    cn <- Loess(views[, k])
    r <- scaleBy(cn, scale(scan_params))
    saveRDS(r, file=lrrFile(viewsLoess)[k])
  }

  ## without loess
  snp_exp <- SnpExperiment(views[, 4:5])
  param <- EmissionParam(temper=0.5)
  fit <- hmm2(snp_exp, param)
  filter_param <- FilterParam()
  g <- cnvFilter(fit, filter_param)

  snp_exp2 <- SnpExperiment(viewsLoess[, 4:5])
  param <- EmissionParam(temper=0.5)
  fit2 <- hmm2(snp_exp2, param)
  g2 <- cnvFilter(fit2, filter_param)
  checkTRUE(identical(state(g), state(g2)))
}
