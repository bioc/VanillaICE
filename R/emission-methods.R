setMethod("gtEmission", c("numeric", "HmmOptions"),
          function(x, options, confidence, featureNames, ...){
            gtEmission(as.integer(x), options, confidence, featureNames, ...)
          })

setMethod("gtEmission", c("integer", "HmmOptions"),
          function(x, options, confidence, featureNames, ...){
            gte <- x; rm(x)
            fn <- featureNames
            if(options@SnpClass == "SnpCopyNumberSet") return()
            states <- options@states
            S <- length(states)
            
            ##gte: genotype estimates
            ##gte <- as.vector(gte)
            gte[gte == 3] <- 1
            gte[gte == 4 | is.na(gte)] <- 3 ##Make missing genotypes have value 3

            if(!options@gt.ICE) {
              gte.state <- options@gte.state
              gte.state <- t(cbind(gte.state, 1-gte.state, NA))
              gt.emission <- gte.state[gte, ]
            }  else {
              gte <- as.integer(gte)
              names(gte) <- fn
              gt.emission <- .gtEmission.ICE(x=gte,
                                             confidence=confidence,
                                             options=options)
            }
            log(gt.emission)
          })

setMethod(".gtEmission.ICE", c("integer", "HmmOptions"),
          function(x, options, confidence, ...){  ##, P.CHOM.Normal, P.CHOM.LOH, SAMPLE=1){
            gte <- x; rm(x)
            hapmapP <- list()
            pkgs <- strsplit(options@annotation, ",")[[1]]
            fn <- names(gte)

            .hapmapProbabilities <- function(annotationPackage){
              require("callsConfidence") || stop("callsConfidence package not available")
              get(switch(annotationPackage,
                         pd.mapping50k.hind240=data(hindPhat),
                         pd.mapping50k.xba240=data(xbaPhat),
                         pd.mapping250k.nsp=data(nspPhat), ##need to calculate phats for 500k
                         pd.mapping250k.sty=data(styPhat)))
            }            
  
            if(length(pkgs) > 1){
              hapmapP[[1]] <- .hapmapProbabilities(pkgs[1])
              hapmapP[[2]] <- .hapmapProbabilities(pkgs[2])
            } else {
              hapmapP <- .hapmapProbabilities(pkgs)
            }
            names(hapmapP) <- pkgs
            if(length(pkgs) > 1){
              require(RSQLite) || stop("RSQLite package not available")
              sql <- "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%'"
              object1 <- object2 <- options
              annotation(object1) <- pkgs[1]
              annotation(object2) <- pkgs[2]
              tmp1 <- dbGetQuery(db(object1), sql)
              tmp2 <- dbGetQuery(db(object2), sql)

              idx.pkg1 <- match(tmp1[["man_fsetid"]], fn)
              idx.pkg1 <- idx.pkg1[!is.na(idx.pkg1)]

              enzyme <- rep(NA, length(gte))
              names(enzyme) <- fn
              enzyme[idx.pkg1] <- pkgs[1]
              
              idx.pkg2 <- match(tmp2[["man_fsetid"]], fn)
              idx.pkg2 <- idx.pkg2[!is.na(idx.pkg2)]
              enzyme[idx.pkg2] <- pkgs[2]

              tmp1 <- .getCallEmission(x=gte[enzyme == pkgs[1]],
                                       confidence=confidence[enzyme == pkgs[1]],
                                       hapmapP=hapmapP[[1]],
                                       options=options)
              tmp2 <- .getCallEmission(x=gte[enzyme == pkgs[2]],
                                       confidence=confidence[enzyme == pkgs[2]],
                                       hapmapP=hapmapP[[2]],
                                       options=options)
              tmp <- rbind(tmp1, tmp2)
              idx <- match(fn, rownames(tmp))
              gt.emission <- tmp[idx, ]
              stopifnot(identical(rownames(gt.emission), fn))
            } else{
              gt.emission <- .getCallEmission(x=gte,
                                              confidence=confidence,
                                              hapmapP=hapmapP,
                                              options=options)
            } 
            gt.emission
          })

##If there is more than one sample in object, it uses only the first.
setMethod("cnEmission", c("numeric", "HmmOptions"),
          function(x, options, confidence, robustSE, ...){
            if(options@SnpClass == "SnpCallSet") return()            
            cne <- x; rm(x)
            fn <- names(cne)
            states <- options@states
            S <- length(states)
            cn.location <- options@cn.location
            cne <- matrix(cne, nrow=length(cne), ncol=S, byrow=FALSE)
            
            ##assume true copy number mean is the same for all samples
            cn.location <- matrix(cn.location, nrow=nrow(cne), ncol=S, byrow=TRUE)

            if(!options@cn.ICE){
              ##Use robust estimate of standard error-- sample-specific
              cn.SE <- matrix(robustSE, nrow=nrow(cn.location), ncol=ncol(cn.location), byrow=TRUE)
            } else{
              ##Use SNP-specific standard errors
              cn.SE <- 1/confidence
              cn.SE <- matrix(cn.SE, nrow=nrow(cne), ncol=ncol(cne), byrow=FALSE)
              colnames(cn.SE) <- states
              if(any(is.na(cn.SE))) stop("NA's in confidence scores.  Must  exclude these SNPs or plug in values for the confidence scores")
            }
            emission.cn <- dnorm(cne, cn.location, cn.SE)
            ##rownames(emission.cn) <- names(cne)
            ##colnames(emission.cn) <- states
            ##length should be R*S (R = number of SNPs, S = number of states)
            log(emission.cn)
          })
