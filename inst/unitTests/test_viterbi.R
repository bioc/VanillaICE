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

lik <- function(r, b, g, answer, r_means, r_sd, b_sd){
  cn <- as.integer(answer)
  lik_cn <- dnorm(r, mean=r_means[cn], r_sd)
  lik_baf <- rep(NA, length(b))
  lik_baf[g==1] <- VanillaICE:::dtrnorm(b[g==1], 0, b_sd)
  lik_baf[g==3] <- VanillaICE:::dtrnorm(b[g==3], 1, b_sd)
  lik_baf[g==2 & cn==3] <- VanillaICE:::dtrnorm(b[g==2 & cn==3], 0.5, b_sd)
  lik_baf[cn==5] <- dnorm(b[cn==5], 1/3, b_sd)
  lik_baf[cn==2] <- VanillaICE:::dtrnorm(b[cn==2], 0, b_sd)
  lik_baf[cn==1] <- dunif(b[cn==1], 0, 1)
  loglik_true <- sum(log(lik_cn*lik_baf))
}


test_updating <- function(){
  ## update mean and standard deviation of BAFs
  library(GenomicRanges)
  oligoset <- getOligoset()
  set.seed(123)
  rtrnorm <- VanillaICE:::rtrnorm
  b <- c(rtrnorm(100, 0, 0.05),
         rtrnorm(100, 0.53, 0.05),
         rtrnorm(100, 1, 0.05))
  r <- rnorm(300, 0.1, 0.1)
  e_param <- EmissionParam()
  assay_list <- list(r, b)
  emissions <- calculateEmission(assay_list, e_param)
  t_param <- TransitionParam(taup=1e10, taumax=1)
  transition_prob <- calculateTransitionProbability(start(oligoset)[seq_along(b)], t_param)
  hmm_param <- HmmParam(emission=emissions,
                        emission_param=e_param,
                        transition=transition_prob,
                        chrom=rep("chr1", length(b)),
                        loglik=LogLik(tolerance=0.05))
  while(doUpdate(hmm_param)){
    hmm_param <- baumWelchUpdate(hmm_param, assay_list)
  }
  ## should be at least 3 updates
  checkTrue(length(loglik(hmm_param)) > 3)
  lim <- c(-3,2)
  x <- seq(lim[1], lim[2], by=0.001)
  means <- VanillaICE:::CN_PRIOR_MEANS()
  sds <- VanillaICE:::CN_PRIOR_SDS()
  plotDensity <- function(means, sds,  lim, col="black", lwd=1){
    x <- seq(lim[1], lim[2], length.out=1000)
    p <- mapply(dnorm, mean=means, sd=sds, MoreArgs=list(x=x))
    for(j in seq_len(ncol(p))){
      points(x, p[, j], pch=20, cex=0.2, lwd=lwd, col=col)
    }
  }
  if(FALSE){
    hist(r, breaks=100, col="grey", border="gray", xlim=lim)
    plotDensity(means, sds, lim=c(-3,2), col="black")
    plotDensity(cn_means(e_params4),
                cn_sds(e_params4), lim=c(-3,2), col="blue")
  }
  fit <- viterbi(hmm_param)
  checkTrue(mean(statei(fit) != as.integer(Rle(c(4,3,4), rep(100,3)))) < 0.005)
  expected_means <- setNames(VanillaICE:::BAF_PRIOR_MEANS(), names(baf_means(e_param)))
  expected_means[c("A", "AB", "B")] <- c(0, 0.53, 1)
  checkEquals(baf_means(emissionParam(hmm_param)), expected_means, tolerance=0.04)
  expected_means <- VanillaICE:::expandCnTwo(VanillaICE:::CN_PRIOR_MEANS())
  expected_means[3:4] <- 0.1
  checkEquals(cn_means(emissionParam(hmm_param)), expected_means, tolerance=0.02)
  ##
  ## b | genotype=A ~ gamma(1, 50)
  ## b | genotype=B ~ gamma(50, 1)
  ##
  ## mean: alpha/beta
  ## var: alpha/beta^2
  ##
  ## -- could we update this from the data
  ## sd = sqrt(alpha/beta^2)
  ## sd = sqrt(mean * 1/beta)
  ## sd2 = mean * 1/beta
  ## beta = mean/sd2
  ## alpha=beta*mean
  ## alpha=mean/sd2*mean
  ## alpha=1/sd2
  ##
  ##
  ## I think there is too much shrinkage of AB
  ##
  ##   for non-two states, the means should be at the prior means
  ##   for cn-two states, means should be near the true values
}

test_Viterbi <- function(){
  library(IRanges)
  library(VanillaICE)
  checkTrue(validObject(Viterbi()))
  checkTrue(validObject(Viterbi(state=integer())))
  checkTrue(validObject(Viterbi(state=2L)))
  checkTrue(validObject(state(Viterbi())))
}

test_baf_emission <- function(){
  set.seed(123)
  library(GenomicRanges)
  param <- EmissionParam()
  genotypes <- sample(1:3, 1000, replace=TRUE)
  mus <- c(0, 0.5, 1)[genotypes]
  sds <- c(0.01, 0.01, 0.01)[genotypes]
  b <- VanillaICE:::rtrnorm(1000, mus, sds)
  emit <- VanillaICE:::.calculateBAFEmission(b)

  head(cbind(b, emit))
  ## the emissions
  val <- apply(emit, 1, which.max)
  cn <- rep(NA, length(genotypes))
  cn[genotypes==1 | genotypes==3] <- 2L
  cn[genotypes==2] <- 3L
  ## data is super clean
  checkIdentical(mean(val!=cn), 0)

  hmm_param <- HmmParam(emission=emit)
  res <- calculateViterbi(hmm_param)
  checkIdentical(state(res), rep(3L,1000))

  transition_param <- TransitionParam(taup=1e12, taumax=1)
  emission_param <- EmissionParam(baf_sds=rep(0.02, 7))
  emit <- VanillaICE:::.calculateBAFEmission(b, emission_param)
  transition_pr <- calculateTransitionProbability(10e3, transition_param)
  hmm_param <- HmmParam(emission=emit,
                        transition=rep(transition_pr, nrow(emit)))
  res <- calculateViterbi(hmm_param)
  checkIdentical(state(res), rep(3L, 1000))
}

test_emission_for_list <- function(){
  set.seed(123)
  genotypes <- sample(1:3, 1000, replace=TRUE)
  mus <- c(0, 0.5, 1)[genotypes]
  sds <- c(0.01, 0.01, 0.01)[genotypes]
  b <- VanillaICE:::rtrnorm(1000, mus, sds)
  r <- rnorm(1000, 0, 0.2)

  ## without specifying parameters
  emissions <- calculateEmission(list(r,b))
  eview <- VanillaICE:::EmissionView(cbind(b, r), emissions)
  val <- apply(emissions, 1, which.max)
  checkTrue(mean(val != 3 & val != 4) < 0.2)
  index <- which(val != 3 & val != 4)
  head(eview[index, ])
  hmm_param <- HmmParam(emission=emissions)
  res <- calculateViterbi(hmm_param)
  checkIdentical(state(res), rep(3L, 1000))
  ##state(res)

  ## with parameters set
  param <- EmissionParam()
  sds <- baf_sds(param)
  sds[c(1,7)] <- c(0.02)
  sds[2:6] <- 0.05
  baf_sds(param) <- sds
  cn_sds(param) <- rep(mad(r), 6)
  emissions <- calculateEmission(list(r,b), param)
  eview <- VanillaICE:::EmissionView(cbind(b, r), emissions)
  val <- apply(emissions, 1, which.max)
  checkTrue(mean(val != 3 & val != 4) < 0.14)

  r <- rnorm(1000, 0, 0.1)
  emissions <- calculateEmission(list(r,b), param)
  eview <- VanillaICE:::EmissionView(cbind(b, r), emissions)
  val <- apply(emissions, 1, which.max)
  checkTrue(mean(val != 3 & val != 4) < 0.1)
  index <- which(val != 3 & val != 4)
  head(eview[index, ])
  hmm_param <- HmmParam(emission=emissions)
  fit <- calculateViterbi(hmm_param)
  checkIdentical(state(res), rep(3L, 1000))
}

test_summarized_exp <- function(){
  library(GenomicRanges)
  library(oligoClasses)
  genotypes <- sample(1:3, 1000, replace=TRUE)
  mus <- c(0, 0.5, 1)[genotypes]
  sds <- c(0.01, 0.01, 0.01)[genotypes]
  b <- as.matrix(VanillaICE:::rtrnorm(1000, mus, sds))
  r <- as.matrix(rnorm(1000, 0, 0.2))
  colnames(b) <- colnames(r) <- "a"
  rowdata <- GRanges(Rle("chr1", 1000),
                     IRanges(seq(1, 100e6, by=100e6/1000), width=25))
  se <- SummarizedExperiment(assays=SimpleList(cn=r,
                               baf=b),
                             rowData=rowdata)
  emit <- calculateEmission(se)
  param <- HmmParam(emission=emit)
  checkIdentical(state(calculateViterbi(param)),
                 rep(3L, 1000))
}

test_oligoset_comparison <- function(){
  ## SummarizedExperiment
  library(oligoClasses)
  library(VanillaICE)
  states <- as.integer(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3))
  nmarkers <- as.integer(c(996, 102, 902, 50, 2467, 102, 76, 1822,
                           99, 900, 20, 160))
  statepath <- rep(states, nmarkers)
  oligoset <- getOligoset()
  ## produces an error -- can't find replacement method for baf
  copyNumber(oligoset)[c(5, 6), ] <- NA
  fit <- hmm(oligoset, is.log=FALSE,
             TAUP=1e10,
             p.hom=1,
             cnStates=c(0.5, 1.5, 2, 2, 2.5, 3.2))
  checkEquals(state(fit[[1]]), states)

  r <- log2(copyNumber(oligoset)/100/2)
  r[is.na(r)] <- 0
  emissions <- calculateEmission(r[,1])
  hmm_param<- HmmParam(emission=emissions)
  ## close to right answer
  state(calculateViterbi(hmm_param))

  b <- baf(oligoset)[,1]/1000
  emissions <- calculateEmission(list(r[,1], b))
  hmm_param<- HmmParam(emission=emissions)
  state(calculateViterbi(hmm_param))
  answer <- IRanges::Rle(as.integer(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3)),
                         as.integer(c(996, 102, 902, 50, 2467, 102, 76, 1822,
                                      99, 900, 20, 160)))


  cn <- as.integer(answer)
  sds <- 0.2
  means <- c(-2, -0.5, 0, 0, 0.5, 1)
  r <- rnorm(length(cn), mean=means[cn], sds)
  g <- sample(1:3, length(cn), replace=TRUE)
  g[cn==4] <- sample(c(1,3), sum(cn==4), replace=TRUE)
  bmeans <- baf_means(EmissionParam())
  b <- rep(NA, length(cn))
  rtrnorm <- VanillaICE:::rtrnorm
  b[g==1] <- rtrnorm(sum(g==1), 0, 0.05)
  b[g==3] <- rtrnorm(sum(g==3), 1, 0.05)
  b[g==2 & cn==3] <- rnorm(sum(g==2&cn==3), 0.5, 0.05)
  b[cn==5] <- rnorm(sum(cn==5), 1/3, 0.05)
  b[cn==2] <- rtrnorm(sum(cn==2), 0, 0.02)

  e_param <- EmissionParam()
  emissions <- calculateEmission(list(r, b), e_param)
  ##Not picking up the LOH state -- had to do with transition probabilities
  t_param <- TransitionParam(taup=1e10, taumax=1)
  transition_prob <- calculateTransitionProbability(position(oligoset), t_param)
  hmm_param<- HmmParam(emission=emissions,
                       transition=transition_prob)
  state1 <- state(calculateViterbi(hmm_param))
  state2 <- IRanges::Rle(cn)
  checkTrue(mean(state1 != state2) < 0.005)
  eview <- VanillaICE:::EmissionView(cbind(r,b), emissions)
  head(eview[cn==4, ])
  colSums(eview[cn==4, ])
}



test_emission_update <- function(){
  library(oligoClasses)
  oligoset <- getOligoset()
  xx <- seq(max(position(oligoset))+10e3, by=50, length.out=50)
  positions <- c(position(oligoset), xx)
  answer <- IRanges::Rle(as.integer(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3, 1)),
                         as.integer(c(996, 102, 902, 50, 2467, 102, 76, 1822,
                                      99, 900, 20, 160, 50)))
  r_means <- c(-2, -0.5, 0, 0, 0.5, 1)
  b_sd <- 0.05
  r_sd <- 0.2
  set.seed(1)
  datlist <- simulateData(answer, r_means, b_sd, r_sd)
  r <- datlist[[1]]
  b <- datlist[[2]]
  g <- datlist[[3]]
  if(FALSE){
    par(mfrow=c(2,1), mar=c(0.5,4, 0.5, 3), las=1)
    plot(positions, r, pch=".", xaxt="n", xaxs="i")
    plot(positions, b, pch=".", xaxt="n", xaxs="i")
  }
  ## Does not take into account initial state probs and transition
  ## probs
  loglik_true <- lik(r, b, g, answer, r_means, r_sd, b_sd)
  ## Run HMM at true values  (over-parameterized because not all BAF and CN states occur)
  e_param_true <- EmissionParam(cn_means=r_means,
                                cn_sds=rep(r_sd, 6),
                                baf_means=c(0, 1/4, 1/3, 1/2, 2/3, 3/4, 1),
                                baf_sds=rep(b_sd, 7))
  emissions <- calculateEmission(list(r, b), e_param_true)
  t_param <- TransitionParam(taup=1e10, taumax=1)

  transition_prob <- calculateTransitionProbability(positions, t_param)
  hmm_param <- HmmParam(emission=emissions,
                        transition=transition_prob)
  ##trace(viterbi, browser)
  fit_true <- calculateViterbi(hmm_param)
  true_model_lik <- loglik(fit_true)

  e_param <- EmissionParam()
  emissions <- calculateEmission(list(r, b), e_param)
  ##Not picking up the LOH state -- had to do with transition probabilities
  t_param <- TransitionParam(taup=1e10, taumax=1)
  transition_prob <- calculateTransitionProbability(positions, t_param)
  hmm_param<- HmmParam(emission=emissions,
                       transition=transition_prob)
  fit <- calculateViterbi(hmm_param)
  loglik(fit)
  llr <- rep(NA, 5)
  for(i in seq_along(llr)){
    e_param <- updateParam(list(r, b), e_param, fit)
    emission(hmm_param) <- calculateEmission(list(r, b), e_param)
    fit2 <- calculateViterbi(hmm_param)
    llr[i] <- loglik(fit2)
    e_param <- e_param
  }
  ##checkTrue(all(diff(llr) > 0)) ## Not sure if this is guaranteed
  checkTrue(mean(statei(fit2) != as.integer(answer)) < 0.001)
  checkEquals(baf_sds(e_param), baf_sds(e_param_true), tolerance=0.1)
  checkEquals(cn_means(e_param), cn_means(e_param_true), tolerance=0.5)
}

test_GenomeAnnotatedDataFrameCoercision <- function(){
  library(Biobase)
  oligoset <- getOligoset()
  ## assumes probes are 25bp
  gr <- as(featureData(oligoset), "SnpGRanges")
  checkTrue(validObject(gr))
  se <- as(oligoset, "SnpArrayExperiment")
  checkTrue(validObject(se))
}


test_SnpArrayExperiment <- function(){
  library(oligoClasses)
  library(GenomicRanges)
  checkTrue(validObject(SnpDataFrame()))
  checkTrue(validObject(SnpDataFrame(DataFrame())))
  sdf <- SnpDataFrame()
  checkException(SnpDataFrame(x=3))
  checkException(SnpDataFrame(DataFrame(x=3)))
  checkTrue(validObject(SnpDataFrame(DataFrame(x=3, isSnp=FALSE))))
  df <- DataFrame(x=3, isSnp=TRUE)
  checkTrue(validObject(sdf <- as(df, "SnpDataFrame")))
  sdf2 <- SnpDataFrame(DataFrame(x=3, isSnp=TRUE))
  checkIdentical(sdf2, sdf)
  df <- SnpDataFrame(DataFrame(x=3, isSnp=FALSE))
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
  rowdata <- GRanges(Rle("chr1", 5), IRanges(1:5, width=1))
  checkException(validObject(SnpArrayExperiment(x, y)))
  checkException(validObject(SnpArrayExperiment(x, y, rowData=rowdata)))
  checkTrue(validObject(SnpArrayExperiment(x, y, rowData=rowdata, isSnp=rep(FALSE,5))))
  rowdata$isSnp <- rep(FALSE,5)
  checkTrue(validObject(SnpArrayExperiment(x, y, rowData=rowdata)))
  se <- SnpArrayExperiment(x, y, rowdata)
  checkIdentical(isSnp(se), rep(FALSE, 5))


  ##
  oligoset <- getOligoset()
  rowdata <- GRanges(Rle("chr1", nrow(oligoset)),
                     IRanges(position(oligoset), width=25L),
                     isSnp=isSnp(oligoset))
  answer <- IRanges::Rle(as.integer(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3)),
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
  ##cn_assays <- snpArrayAssays(cn=r, baf=b)
  validObject(SnpArrayExperiment(cn=r, baf=b, rowData=rowdata))
}

test_hmm2 <- function(){
  library(oligoClasses)
  library(GenomicRanges)
  oligoset <- getOligoset()
  fit1 <- hmm(oligoset,
              is.log=FALSE,
              TAUP=1e10,
              p.hom=1,
              cnStates=c(0.5, 1.5, 2, 2, 2.5, 3.2))
  ##trace(hmm2, signature="oligoSnpSet", browser)
  checkException(hmm2(oligoset))
  oligoset2 <- oligoset
  copyNumber(oligoset2) <- log2(copyNumber(oligoset2)/100/2)
  baf(oligoset2) <- baf(oligoset2)/1000
  ##trace(.hmm, browser)
  fit2 <- hmm2(oligoset2)
  sp1 <- rep(fit1[[1]]$state, fit1[[1]]$numberProbes)
  sp2 <- as.integer(Rle(statei(fit2), fit2$numberFeatures))
  checkTrue(mean(sp1==sp2) > 0.999)
  if(FALSE) plot(position(oligoset2), copyNumber(oligoset2), pch=".")
}

.test_streamlined_hmm <- function(){
  library(oligoClasses)
  oligoset <- getOligoset()
  se <- as(oligoset, "SnpArrayExperiment")
  assays(se)[["cn"]] <- lrr(se)/100
  assays(se)[["baf"]] <- baf(se)/1000
  hmm2(se)
  xx <- seq(max(position(oligoset))+10e3, by=50, length.out=50)
  positions <- c(position(oligoset), xx)
  rowdata <- GRanges(Rle("chr1", length(positions)),
                     IRanges(positions, width=25L))
  answer <- IRanges::Rle(as.integer(c(3, 4, 3, 5, 3, 2, 3, 3, 2, 3, 2, 3, 1)),
                         as.integer(c(996, 102, 902, 50, 2467, 102, 76, 1822,
                                      99, 900, 20, 160, 50)))
  r_means <- c(-2, -0.5, 0, 0, 0.5, 1)
  b_sd <- 0.05
  r_sd <- 0.2
  datlist <- simulateData(answer, r_means, b_sd, r_sd)
  r <- as.matrix(datlist[[1]])
  b <- as.matrix(datlist[[2]])
  g <- as.matrix(datlist[[3]])
  colnames(r) <- colnames(b) <- colnames(g) <- "a"
##  se <- CopyNumberExperiment(assays=SimpleList(cn=r,
##                               baf=b,
##                               gt=g),
##                             rowData=rowdata)
##  trace(hmm2, signature="SummarizedExperiment", browser)
  fit <- hmm2(se)
  e_param <- EmissionParam()
  emissions <- calculateEmission(list(r, b), e_param)
  ##Not picking up the LOH state -- had to do with transition probabilities
  t_param <- TransitionParam(taup=1e10, taumax=1)
  transition_prob <- calculateTransitionProbability(positions, t_param)
  hmm_param<- HmmParam(emission=emissions,
                       transition=transition_prob)
  fit <- calculateViterbi(hmm_param)

  x <- filters(b, r, starts)
  b <- x[[1]]
  r <- x[[2]]
  starts <- x[[3]]

  e_param <- EmissionParam()
  emissions <- calculateEmission(list(r, b), e_param)
  t_param <- TransitionParam(taup=1e10, taumax=1)
  transition_prob <- calculateTransitionProbability(, t_param)
  hmm_param<- HmmParam(emission=emissions,
                       transition=transition_prob)
  LL <- rep(NA, 10)
  delta <- 1
  i <- 1
  while(delta > 0.05){
    fit <- calculateViterbi(hmm_param)
    LL[i] <- loglik(fit)


    e_param <- updateParam(list(r, b), e_param, fit)
    emission(hmm_param) <- calculateEmission(list(r, b), e_param)
    if(i > 1) delta <- LL[i]-LL[i-1]
    if(i == 10) break()
    i <- i+1
  }
}

test_multiple_chromosomes <- function(){
  library(oligoClasses)
  se <- getSE()
  se <- sweepMode(se, MARGIN=2)
  checkEquals(colModes(copyNumber(se)), 0, tolerance=0.005)
  fit1 <- hmm2(se)

  rd <- rowData(se)
  index <- (length(rd)-100):length(rd)
  chrom <- chromosome(rd)
  chrom[index] <- "chr2"

  ##
  ## Force a distance between chromosomes start positions that will be
  ## negative, but for which the transition probability method for
  ## SnpArrayExperiment will put a weak transition probability
  ##
  start(rd)[index] <- start(rd)[seq_along(index)]
  end(rd)[index] <- end(rd)[seq_along(index)]
  rd <- GRanges(chrom, IRanges(start(rd), end(rd)), isSnp=rep(TRUE, length(rd)))
  rd2 <- SnpGRanges(rd)
  rowData(se) <- rd2
  fit2 <- hmm2(se)
  checkIdentical(state(fit2[[1]])[-length(state(fit2[[1]]))], state(fit1[[1]]))
}


test_cn_NAs <- function(){
  x <- matrix(rnorm(10*6), 10, 6)
  x[c(4,8), ] <- NA
  x <- as.numeric(x)
  NArows <- sum(is.na(x)/6)

  library(GenomicRanges)
  oligoset <- getOligoset()
  set.seed(123)
  rtrnorm <- VanillaICE:::rtrnorm
  b <- c(rtrnorm(100, 0, 0.05),
         rtrnorm(100, 0.53, 0.05),
         rtrnorm(100, 1, 0.05))
  r <- rnorm(300, 0.1, 0.1)
  r[c(10, 11, 100, 200)] <- NA

  starts <- start(oligoset)[seq_along(b)]
  starts[150] <- NA
  x <- NA_filter(list(b, r, starts))
  checkTrue(all(elementLengths(x) == 295L))

  ## test method on SnpArrayExperiment
  starts[150] <- starts[149]+1
  r <- setNames(as.matrix(r), "a")
  b <- setNames(as.matrix(b), "a")
  rowdata <- SnpGRanges(GRanges(Rle("chr1", 300),
                                IRanges(starts, width=1), isSnp=rep(TRUE, 300)))
  coldata <- DataFrame(row.names="a")
  se <- SnpArrayExperiment(cn=r, baf=b, rowData=rowdata, colData=coldata)
  checkTrue(validObject(se))
  se2 <- NA_filter(se)
  checkIdentical(nrow(se2), 296L)
}



test_pathologies <- function(){
  ##
  ## test shrinkage for big observations sequences
  ##
  se <- readRDS("~/Software/bridge/se_example.rds")
  se <- as(se, "SnpArrayExperiment")
  fit <- hmm2(se)[[1]]
  roi <- GRanges("chr17", IRanges(44165803, 44343902))
  result <- state(fit)[subjectHits(findOverlaps(roi, fit))]
  checkIdentical(as.integer(result), 5L)

  ##
  ## missing values
  ##

  ##
  ## user contract (misspecitified arguments)
  ##


  ##
  ## data not gaussian, extreme values
  ##

  ##
  ## germline mosaicism
  ##

}
