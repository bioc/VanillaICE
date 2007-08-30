chr1Example <- function(na06993){

  set.seed(15)
  chr1 <- na06993[chromosome(na06993) == "1", ]
  chr1 <- chr1[!is.na(chromosome(chr1)), ]
  chr1 <- chr1[order(position(chr1)), ]
  A <- 1001:1100
  B <- 6001:6100
  C <- 8001:8100
  D <- 2001:2200
  E1 <- 9001:9005  
  E2 <- 9101:9103


  calls(chr1) <- matrix(rbinom(nrow(chr1), 1, 0.3) + 1, nc=1)
  cnConfidence(chr1) <-  matrix(rgamma(nrow(chr1), shape=1, rate=2) + 0.3)
  
  ##################################################    
  ##Region A, LOH
  ##################################################      
  calls(chr1)[A] <- rbinom(100, 1, 0.001) + 1
  calls(chr1)[1050:1051] <- 2

  ##################################################  
  ##Region B, deletion
  ##################################################  
  epsilon <- diff(quantile(base::log(copyNumber(chr1), 2), probs=c(0.16, (1-0.16))))/2  
  calls(chr1)[B] <- rbinom(100, 1, 0.001) + 1
  copyNumber(chr1)[B] <- 2^(base::log(1, 2) + rnorm(100, 0, epsilon))
  calls(chr1)[6050] <- 2
  copyNumber(chr1)[6050:6051] <- 2^(base::log(2, 2) + rnorm(2, 0, epsilon/4))
  ##precise estimate simulated by a robust estimate of the copy number standard deviation
  cnConfidence(chr1)[6050:6051] <- epsilon/2


  ##################################################
  ##Region C,  Deletion
  ##################################################  
  calls(chr1)[C] <- 1
  copyNumber(chr1)[C] <- 2^(base::log(1, 2) + rnorm(100, 0, epsilon/2))
  copyNumber(chr1)[8050:8051] <- 2^(base::log(2, 2)) + rnorm(2, 0, epsilon/4)
  ##Noisy estimate
  cnConfidence(chr1)[8050:8051] <- epsilon*2

  ##################################################  
  ##Region D, Amplification
  ##################################################  
  copyNumber(chr1)[D] <- 2^(base::log(3, 2) + rnorm(200, 0, epsilon/2))
  copyNumber(chr1)[2100:2101] <- 2^(base::log(2, 2) + rnorm(2, 0, epsilon/4))
  ##precise estimate simulated by a robust estimate of the copy number standard deviation
  cnConfidence(chr1)[2100:2101] <- epsilon/2
  calls(chr1)[2101] <- 2
  
  ##To avoid simulating a breakpoint in a region of very close
  ##physical proximity (linkage disequilibrium)
##
##  position(chr1)[2102] <- position(chr1)[2103] - 1000
##  position(chr1)[2100:2101] <- position(chr1)[2100:2101] + 20000

  ##################################################  
  ##Region E, microdeletion, microamplification
  ##################################################
  copyNumber(chr1)[E1] <- 2^(base::log(1, 2) + rnorm(5, 0, epsilon/4))  
  copyNumber(chr1)[E2] <- 2^(base::log(3, 2) + rnorm(3, 0, epsilon/4))
  calls(chr1)[E1] <- 1
  cnConfidence(chr1)[E1] <- epsilon/2
  cnConfidence(chr1)[E2] <- epsilon/2  

  isXbaHom <- enzyme(chr1) == "Xba" & (calls(chr1) == 1 | calls(chr1) == 3)
  isXbaHet <- enzyme(chr1) == "Xba" & calls(chr1) == 2  
  isHindHom <- enzyme(chr1) == "Hind" & (calls(chr1) == 1 | calls(chr1) == 3)
  isHindHet <- enzyme(chr1) == "Hind" & calls(chr1) == 2

  ##We sample from confidence scores in the hapmap samples under the
  ##assumption that the genotype calls are correct
  require(callsConfidence) || stop("callsConfidence not available")
  data(xbaPhat)
  group <- xbaPhat[, 1] ##see referenceProbabilities.R in callsConfidence package
  
  callsConfidence(chr1)[isXbaHom] <- sample(xbaPhat[group == 3, 2], sum(isXbaHom), replace=FALSE)
  callsConfidence(chr1)[isXbaHet] <- sample(xbaPhat[group == 2, 2], sum(isXbaHet), replace=FALSE)

  data(hindPhat)
  group <- hindPhat[, 1] ##see referenceProbabilities.R in callsConfidence package  
  callsConfidence(chr1)[isHindHom] <- sample(hindPhat[group == 3, 2], sum(isHindHom), replace=FALSE)
  callsConfidence(chr1)[isHindHet] <- sample(hindPhat[group == 2, 2], sum(isHindHet), replace=FALSE)  
  ##save(chr1, file="~/projects/hmm/HmmDataSets/data/chr1.rda", compress=TRUE)
  chr1
}
