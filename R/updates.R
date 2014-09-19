prC1 <- function(pr) 1/2*pr[, "A"] + 1/2*pr[, "B"]
prC2 <- function(pr) 1/4*pr[, "A"] + 1/2*pr[, "AB"] +1/4*pr[, "B"]
prC3 <- function(pr) 1/8*pr[, "A"] + 3/8*pr[, "AAB"] + 3/8*pr[, "ABB"] + 1/8*pr[, "B"]
prC4 <- function(pr) 1/16*pr[, "A"] + 4/16*pr[, "AAAB"] +
  6/16*pr[, "AB"] + 4/16*pr[, "ABBB"] + 1/16*pr[, "B"]



isNonDecreasing <- function(x) any(diff(x) < 0)

shrink <- function(means, sds, state_freq,
                   prior_means,
                   prior_obs,
                   prior_sd){
  ## use the prior if any NaNs
  if(any(is.na(means))) means[is.na(means)] <- prior_means[is.na(means)]
  if(any(is.na(sds))) sds[is.na(sds)] <- prior_sd[is.na(sds)]
  ##
  ## posterior for means
  ##
  prior_prec <- prior_obs/prior_sd^2 ## n0 prior observations
  data_prec <- state_freq/sds^2
  posterior_means <- prior_prec/(prior_prec + data_prec) * prior_means + data_prec/(prior_prec + data_prec) * means
  ##
  ## posterior for precision
  ##    sigma2 ~ gamma(nu_n/2, nu_n*sigma^2_n/2)
  ##    n = # data observations
  ##    nu_0 = # prior observations
  ##    nu_n = nu_0 + n
  ##    kappa_n = kappa_0 + n
  ##    kappa_0 = # prior observations for mean mu0
  ##    mu0 = prior mean
  ##    sigma2_n = 1/nu_n * [nu_0*sigma^2_0 + (n-1)s^2 + kappa0*n/kappa_n(means - mu0)^2]
  ##
  ##  For states with 0 observations, the posterior should be the same
  ##  as the prior
  sigma2_0 <- prior_sd^2
  nu_0 <- kappa_0 <- prior_obs
  kappa_n <- nu_n <- prior_obs + state_freq
  nn <- pmax(state_freq, 0)
  sigma2_n <- 1/nu_n * (nu_0*sigma2_0 + pmax((nn-1),0)*sds^2 + kappa_0 * nn/kappa_n*(means - prior_means)^2)
  ##
  if(length(sigma2_n)==5){
    pooled <- sum(state_freq[-1]*sigma2_n[-1])/sum(state_freq[-1])
    sigma2_n[-1] <- pooled
  } else {
    pooled <- sum(state_freq[-c(1, 7)]*sigma2_n[-c(1,7)])/sum(state_freq[-c(1,7)])
    sigma2_n[-c(1,7)] <- pooled
  }
  posterior_sd <- sqrt(sigma2_n)
  list(means=posterior_means,
       sds=posterior_sd)
}

cn_alleles <- function(cn){
  switch(cn,
         one=c("A", "B"),
         two=c("A", "AB", "B"),
         three=c("A", "AAB", "ABB", "B"),
         four=c("A", "AAAB", "AB", "ABBB", "B"),
         NULL)
}

.update_baf <- function(b, phi, means, scalars, fv){
  phi <- t(t(phi)*scalars)
  row_total <- rowSums(phi, na.rm=TRUE)
  weights <- phi/row_total*fv
  total <- colSums(weights, na.rm=TRUE)
  checkWeight <- as.numeric(sum(weights[,1],na.rm=TRUE)/total[1])
  if(!all.equal(checkWeight, 1)) stop("Weights incorrect")
  sdFun <- function(b, w, mu, T){
    sqrt(sum(w*(b-mu)^2, na.rm=TRUE) / T)
  }
  meanFun <- function(b, w, T){
    sum(w*b, na.rm=TRUE) / T
  }
  mat <- matrix(NA, ncol(phi), 2)
  dimnames(mat) <- list(colnames(phi), c("mean", "sd"))
  for(j in seq_along(means)){
    mat[j, ] <- c(meanFun(b, weights[, j], total[j]),
                  sdFun(b, weights[, j], means[j],
                        total[j]))
  }
  return(mat)
}

updateBafParam <- function(b, param, model){
  means <- baf_means(param)
  sds <- baf_sds(param)
  FV <- forward_backward(model)
  phi <- mapply(dnorm, mean=as.list(means),
                sd=as.list(sds),
                MoreArgs=list(x=b))
  allele_prob <- scalars()
  mean_sd_list <- vector("list", 4)
  copynum <- c("one", "two", "three", "four")
  for(j in seq_along(copynum)){
    cn <- copynum[j]
    alleles <- cn_alleles(cn)
    mean_sd_list[[j]] <- .update_baf(b, phi[, alleles, drop=FALSE],
                                     means[alleles],
                                     allele_prob[[j]],
                                     FV[, cn])
  }
  tmp <- do.call(rbind, mean_sd_list)
  rownames(tmp) <- make.unique(rownames(tmp))
  mean_sd <- tmp[c("A.1", "AAAB", "AAB", "AB", "ABB", "ABBB", "B.1"), ]
  baf_means(param) <- setNames(mean_sd[, 1], names(means))
  baf_sds(param) <- setNames(mean_sd[, 2], names(means))
  ##
  ## shrink
  ##
  param <- shrinkBafParam(param, model)
  param
}

shrinkBafParam <- function(param, model){
  ## Rough estimate of the number of observations for estimating the
  ## mean
  sds <- baf_sds(param)
  means <- baf_means(param)
  nobs <- table(statef(model))
  nobs <- setNames(as.integer(nobs), names(nobs))
  ##
  ## For germline, there should be little shrinkage for the CN=2
  ## params.
  ##
  ##  The HMM does not allow us to directly estimate the number of
  ## observations for the different allele frequencies
  nAB <- max(1/2 * sum(nobs["3"]), 1)
  nB <- nA <- 1/2 * nAB + 1/2*sum(nobs[["4"]])
  ## moderate shrinkage
  nABB <- nAAB <- 3/8 * sum(nobs["5"])
  nABBB <- nAAAB <- 4/16 * sum(nobs["6"])
  ## note, this will not necessarily add up to the observed state
  ## frequencies
  state_freq <- setNames(c(nA, nAAAB, nAAB, nAB, nABB, nABBB, nB),
                         names(means))
  posteriors <- shrink(means=means,
                       sds=sds,
                       state_freq=state_freq,
                       prior_means=BAF_PRIOR_MEANS(),
                       prior_obs=max(1/10*sum(nobs), 25),
                       prior_sd=BAF_SDS())
  baf_means(param) <- posteriors[["means"]]
  baf_sds(param) <- posteriors[["sds"]]
  ##means_new <- makeMusBafNondecreasing(posterior_means)
  ##baf_means(param) <- means_new
  return(param)
}

updateCnParam <- function(cn, param, model){ ##fv, bv, j) {
  ##LIMIT_RANGE <- c(-4, 3)  ## threshold extreme log R ratios
  LIMIT_RANGE <- CN_range(param)
  nr <- length(cn)
  S <- 6L
  means <- cn_means(param)
  sds <- cn_sds(param)
  d <- dnorm(cn, mean=matrix(means, nr, S, byrow=TRUE),
             sd=matrix(sds, nr, S, byrow=TRUE))
  d1 <- d
  FV <- forward_backward(model)
  g1 <- FV * d1
  g2 <- FV * (1-d1)

  totalh <- apply(g1, 2, sum, na.rm=TRUE)
  means2 <- colSums(g1*cn, na.rm=TRUE)/totalh

  sds2 <- rep(NA, length(means))
  for(i in seq_along(sds2))  sds2[i] <- sqrt(sum(g1[, i]*(cn-means2[i])^2, na.rm=TRUE)/totalh[i])

  states <- table(statef(model))
  nobs <- setNames(as.integer(states), names(states))
  ##
  ## For copy number 2, states 3 and 4 are both observations of
  ## diploid copy number
  ##
  stats <- reduceCnTwo(nobs, means2, sds2)
  posteriors <- shrink(means=stats[["mean"]],
                       sds=stats[["sd"]],
                       state_freq=stats[["n"]],
                       prior_means=CN_PRIOR_MEANS(),
                       prior_sd=CN_PRIOR_SDS(),
                       prior_obs=max(1/4*sum(nobs), 50))
  ## a few constraints
  ## if(means[1] > -0.5) means[1] <- -0.5
  ## if(means[5] < 0.1) means[5] <- 0.1
  means <- setNames(posteriors[["means"]], names(means))
  cn_means(param) <- expandCnTwo(means)
  sds <- setNames(posteriors[["sds"]], names(means))
  cn_sds(param) <- expandCnTwo(sds)
  if(FALSE){
    r <- lrr(object)[chromosome(object) %in% paste0("chr", 1:22),1]
    hist(r, breaks=1000, col="gray", border="gray", freq=FALSE)
    mu <- modev(r)
    s <- cn_sds(param)[2]
    quants <- seq(0, 1, 0.001)
    x <- qnorm(quants, mu, s)
    probs <- dnorm(x, mu, s)
    points(x, probs, pch=20, cex=0.2)
    ac <- acf2(r)
  }
  param
}

updateParam <- function(assay_list, param, model){
  param <- updateBafParam(assay_list[[2]], param, model)
  updateCnParam(assay_list[[1]], param, model)
}



reduceCnTwo <- function(nobs, means, sds){
  n_diploid <- nobs[3:4]
  mean_diploid <- sum(n_diploid*means[3:4], na.rm=TRUE)/sum(n_diploid)
  sd_diploid <- sum(n_diploid*sds[3:4], na.rm=TRUE)/sum(n_diploid)
  means2 <- c(means[1:2], mean_diploid, means[5:6])
  sds2 <- c(sds[1:2], sd_diploid, sds[5:6])
  nobs <- c(nobs[1:2], sum(n_diploid), nobs[5:6])
  list(mean=means2, sd=sds2, n=nobs)
}

expandCnTwo <- function(x) x[c(1, 2, 3, 3, 4,5)]
