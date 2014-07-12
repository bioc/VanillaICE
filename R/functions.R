#' @export
threshold <- function(x, lim=c(-Inf,Inf), amount=0){
  notna <- !is.na(x)
  index1 <- which(x <= lim[1] & notna)
  index2 <- which(x >= lim[2] & notna)
  x1 <- runif(length(index1), lim[1], lim[1]+amount)
  x2 <- runif(length(index2), lim[2]-amount, lim[2])
  x[index1] <- x1
  x[index2] <- x2
  return(x)
}

setMethod(NA_index, "list", function(x){
  is_na <- sapply(x, is.na)
  which(rowSums(is_na) > 0)
})

setMethod(NA_filter, "list", function(x, i){
  if(missing(i)) i <- NA_index(x)
  if(length(i) > 0){
    x <- lapply(x, "[", -i)
  }
  x
})

setMethod(NA_filter, "numeric", function(x, i){
  if(missing(i)) i <- NA_index(x)
  if(length(i) > 0){
    x <- x[-i]
  }
  x
})

setMethod(NA_filter, "character", function(x, i){
  if(missing(i)) i <- NA_index(x)
  if(length(i) > 0){
    x <- x[-i]
  }
  x
})

NA_Rle <- function(object, i, nr){
  if(length(i) == 0) return(object)
  ## otherwise, fill in missing
  states <- rep(NA, nr)
  states[-i] <- statei(object)
  state(object) <- as(states, "Rle")
  object
}

autosomes <- function(x) unique(x[!(x %in% c("X", "Y") | x %in% c("chrX", "chrY"))])

paste1 <- function(x) paste(head(round(x,2), 3), collapse=", ")
HMM_STATES <- function() c("TwoCopyLoss",
                           "SingleCopyLoss",
                           "Diploid",
                           "CopyNeutralLOH",
                           "SingleCopyGain",
                           "TwoCopyGain")

CN_MEANS <- function() c(-2, -0.4, 0, 0, 0.4, 1)
CN_SDS <- function() rep(0.3, 6)
BAF_ALLELE_NAMES <- function() c("A", "AAAB", "AAB", "AB", "ABB", "ABBB", "B")
BAF_MEANS <- function() setNames(c(0, 0.1, 1/3, 0.5, 2/3, 0.9, 1),
                                 BAF_ALLELE_NAMES())
BAF_SDS <- function() setNames(rep(0.1, 7),
                               BAF_ALLELE_NAMES())

BAF_PRIOR_MEANS <- function() c(0, 1/4, 1/3, 1/2, 2/3, 3/4, 1)
CN_PRIOR_MEANS <- function() c(-2, -0.5, 0, 0.5, 1)
CN_PRIOR_SDS <- function() c(0.6, rep(0.3, 4))


scalars <- function()
  list(c(1/2, 1/2),
       c(1/4, 1/2, 1/4),
       c(1/8, 3/8, 3/8, 1/8),
       c(1/16, 4/16, 6/16, 4/16, 1/16))

FV_columns <- function() c("zero", "one", "two", "two_", "three", "four")


#' Rescale a numeric vector
#'
#' @param x numeric vector
#' @param l lower limit of rescaled \code{x}
#' @param u upper limit of rescaled \code{x}
#' @export
#' @keywords manip
rescale <- function(x, l, u){
  b <- 1/(u-l)
  a <- l*b
  (x+a)/b
}

getCrlmmReference <- function(x) paste(annotation(x), "Conf", sep="")

isLoaded <- function(dataset, environ=.vanillaIcePkgEnv)
	exists(dataset, envir=environ)

getVarInEnv <- function(dataset, environ=.vanillaIcePkgEnv){
	if (!isLoaded(dataset))
		stop("Variable ", dataset, " not found in .vanillaIcePkgEnv")
	environ[[dataset]]
}

loader <- function(theFile, envir, pkgname){
	theFile <- file.path(system.file(package=pkgname),
			     "extdata", theFile)
	if (!file.exists(theFile))
		stop("File ", theFile, "does not exist in ", pkgname)
	load(theFile, envir=envir)
}


dtrnorm <- function(x, mean, sd, lower=0, upper=1){
  dnorm(x, mean, sd)/(pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
}

trnorm <- function(x){
  function(mu, sigma, upper=1, lower=0){
    dnorm(x, mu, sigma)/(pnorm(upper, mu, sigma) - pnorm(lower, mu, sigma))
  }
}


rtrnorm <- function(n, mean, sd, lower=0, upper=1){
  x <- rnorm(n, mean, sd)
  x[x < lower] <- lower
  x[x > upper] <- upper
  x
}

setMethod(modev, "integer", function(x){
  tab <- table(x)
  as.integer(names(tab)[which.max(tab)])
})

setMethod(modev, "numeric", function(x){
  notna <- !is.na(x)
  xx <- x[notna]
  dens <- density(xx)
  x_mode <- dens$x[which.max(dens$y)]
  return(x_mode)
})

#' @export
rowModes <- function(x) apply(x, 1, modev)

#' @export
colModes <- function(x) apply(x, 2, modev)

#' Compute the median absolute deviation (MAD) for the rows of a matrix
#'
#' @param x matrix
#' @param ... additional arguments to rowMedians
#' @seealso \code{\link{mad}}
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' rowMedians(X)
#' @seealso \code{\link{mad}} \code{\link{rowMedians}}
#' @rdname robust-statistics
#' @export
rowMAD <- function(x, ...)  1.4826*rowMedians(abs(x-rowMedians(x, ...)), ...)

#' Calculate the column-wise median absolute deviation of a matrix
#'
#' @param x a numeric matrix
#' @param takeLog whether to first log2 transform the numeric matrix
#' @param ... additional arguments to \code{mad}
#' @seealso \code{\link{mad}}
#' @rdname robust-statistics
#' @return a matrix of the same dimension as the input. Within a column, the entries are identical
#' @examples
#' data(locusLevelData, package="oligoClasses")
#' sds <- robustSds(locusLevelData[["copynumber"]]/100,
#'                  takeLog=TRUE)
#' @export
robustSds <- function(x, takeLog=FALSE, ...){
  if(takeLog) x <- log2(x)
  sds <- apply(x, 2, "mad", ...)
  sds <- matrix(sds, nrow(x), ncol(x), byrow=TRUE)
  dimnames(sds) <- dimnames(x)
  return(sds)
}
