genomestudio_data <- function(){
  dir <- "~/Ingo/ForExperimentDataPackage"
  extdir <- "~/Software/Illumina610SnpData/inst/extdata"
  gsfiles <- list.files(dir, pattern="FinalReport", full.names=TRUE)
  map_file <- file.path(dir, "Beaty_610Q_release_SNP Table.csv")
  features <- fread(map_file, nrows=10)
  features <- suppressWarnings(fread(map_file, select=c(3, 5:7)))
  features <- features[features$Chr %in% c(1:22, "X", "Y"), ]
  tmp <- as.data.frame(features)

  set.seed(123)
  index <- sample(seq_len(nrow(tmp)), 60e3)
  tmp <- tmp[index, ]
  rownames(tmp) <- NULL
  write.csv(tmp, file=file.path(extdir, "SNP_info.csv"),
            quote=FALSE, row.names=FALSE)

  files <- list.files(dir, full.names=TRUE, recursive=TRUE, pattern="FinalReport")
  outfiles <- gsub("Beaty_610Q_release_0409__", "", basename(files))

  for(j in seq_along(files)){
    dat <- fread(files[j], nrows=5)
    ##  looking at dat, we want to keep the following columns
    select <- c(1, 3, 12, 13, 20, 21)
    colnames(dat)[select]
    dat <- fread(files[j], select=select)
    dat <- as.data.frame(dat)
    dat <- dat[dat[[1]] %in% tmp$Name, ]
    write.csv(dat, file=file.path(extdir, outfiles[j]), quote=FALSE, row.names=FALSE)
  }
}
