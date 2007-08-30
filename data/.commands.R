library(ICE)
library(HmmDataSets)
data(na06993, package="HmmDataSets")

source("~/projects/Rpacks/ICE/data/.functions.R")
chromosome1 <- chr1Example(na06993)
hset <- as(chromosome1, "HmmSnpSet")
hset@log2 <- TRUE
distance(hset) <- TRUE
hset@locationCopyNumber <- c(log2(1), log2(2), log2(2), log2(3))
pCHOM(hset) <- c(0.999, 0.7, 0.999, 0.7)
copyNumberIce(hset) <- TRUE
callsIce(hset) <- TRUE
chromosome1 <- hmm(hset, verbose=TRUE)
save(chromosome1, file="~/projects/Rpacks/ICE/data/chromosome1.RData", compress=TRUE)

###########################################################################
##Illumina example
###########################################################################
library(ICE)
library(TrioData)
data(trioset)
class(trioset)

##colnames of assayData have .tar.gz extension....
#sampleNames(trioset) <- colnames(calls(trioset))
#stopifnot(validObject(trioset))

##update to oligoSnpSet defined in oligoClasses
##tmp <- updateObject(trioset)
##trioset <- tmp
##rm(tmp)
##save(trioset, file="~/projects/thomas/TrioData/data/trioset.RData", compress=TRUE)

hmmset <- as(trioset, "HmmSnpSet")
validObject(hmmset)
chromosome10 <- hmmset[chromosome(hmmset) == "10", grep("L07_2622", sampleNames(hmmset))]
##large deletion on q-arm
locationCopyNumber(chromosome10) <- log2(c(0.5, 1, 1, 1.5))
chromosome10 <- hmm(chromosome10, verbose=TRUE, log2=FALSE)
##save(chromosome10, file="~/projects/Rpacks/ICE/data/chromosome10.rda", compress=TRUE)

source("~/madman/Rpacks/SNPchip/R/AllClasses.R")
source("~/madman/Rpacks/SNPchip/R/methods-AnnotatedSnpCallSet.R")
