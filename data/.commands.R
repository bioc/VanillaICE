###########################################################################
##an oligoSnpSet example  
library(HmmDataSets)
data(na06993, package="HmmDataSets")
chr1 <- na06993[chromosome(na06993) == "1", ]

library(VanillaICE)
fvarLabels(chr1)[2] <- "chromosome"
fvarLabels(chr1)[3] <- "position"
pkgs <- c("pd.mapping50k.hind240", "pd.mapping50k.xba240")
annotation(chr1) <- paste(pkgs[order(pkgs)], collapse=",")

load("data/chromosome1.RData")
idx <- match(featureNames(chromosome1), featureNames(chr1))
chr1 <- chr1[idx, ]
identical(featureNames(chr1), featureNames(chromosome1))

featureData(chr1)$chromosome <- chromosome(chr1)
featureData(chr1)$position <- position(chr1)
identical(chromosome(chr1), fData(chr1)$chrom)

data(cephTriosAffymetrix100k)
ceph <- cephTriosAffymetrix100k; rm(cephTriosAffymetrix100k)
annotation(ceph) <- annotation(chr1)
source("~/madman/Rpacks/SNPchip/R/methods-oligoSnpSet.R")
ceph <- ceph[chromosome(ceph) == "1", ]
load("data/chr1.RData")
idx <- match(featureNames(chr1), featureNames(ceph))
ceph <- ceph[idx, ]
stopifnot(identical(featureNames(ceph), featureNames(chr1)))
copyNumber(chr1) <- log2(copyNumber(chr1))
copyNumber(ceph) <- log2(copyNumber(ceph))
source("R/functions.R")
cnSE <- calculateCnSE(object=chr1, referenceSet=ceph, epsilon=0.1)
cnConfidence(chr1) <- 1/cnSE
save(chr1, file="data/chr1.RData", compress=TRUE)

load("data/chr1.RData")
copyNumber <- cbind(copyNumber(chr1), copyNumber(chr1))
calls <- cbind(calls(chr1), calls(chr1))
cnConfidence <- cbind(cnConfidence(chr1), cnConfidence(chr1))
callsConfidence <- cbind(callsConfidence(chr1), callsConfidence(chr1))
colnames(copyNumber) <- colnames(calls) <- colnames(cnConfidence) <- colnames(callsConfidence) <- letters[1:2]
pd <- rbind(pData(chr1), pData(chr1))
rownames(pd) <- letters[1:2]
pD <- new("AnnotatedDataFrame", data=pd, varMetadata=varMetadata(chr1))
chr1.2 <- new("oligoSnpSet",
              copyNumber=copyNumber,
              calls=calls,
              cnConfidence=cnConfidence,
              callsConfidence=callsConfidence,
              featureData=featureData(chr1),
              phenoData=pD,
              annotation=annotation(chr1))
save(chr1.2, file="data/chr1.2.RData", compress=TRUE)


library("pd.mapping50k.xba240")
###########################################################################
##See the HmmDataSets package for chromosome1?
##Requires version 1.1.4 of VanillaICE to get the chromosome1 object
data(chromosome1)
phenoData <- phenoData(chromosome1)
save(phenoData, file="VanillaICE/inst/phenoData.RData")
featureData <- featureData(chromosome1)
save(featureData, file="VanillaICE/inst/featureData.RData", compress=TRUE)

copyNumber <- copyNumber(chromosome1)
write.csv(copyNumber, file="VanillaICE/inst/copyNumber.csv")

calls <- calls(chromosome1)
write.csv(calls, file="VanillaICE/inst/calls.csv")

cnConfidence <- 1/cnConfidence(chromosome1) ##confidence is inverse of standard error
write.csv(cnConfidence, file="VanillaICE/inst/cnConfidence.csv")
          
callsConfidence <- callsConfidence(chromosome1)
write.csv(callsConfidence, file="VanillaICE/inst/callsConfidence.csv")

chromosome1 <- new("oligoSnpSet",
                   copyNumber=copyNumber,
                   calls=calls,
                   cnConfidence=cnConfidence,
                   callsConfidence=callsConfidence)
featureData(chromosome1) <- featureData
phenoData(chromosome1) <- phenoData
annotation(chromosome1) <- "pd.mapping50k.hind240,pd.mapping50k.xba240"
validObject(chromosome1)
save(chromosome1, file="VanillaICE/data/chromosome1.RData", compress=TRUE)
###########################################################################
##a HmmParameter example


###########################################################################
##a HmmPredict example

