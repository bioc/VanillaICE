setMethod("db", "HmmOptions",
          function(object) {
            require(annotation(object), character.only=TRUE) || stop(paste(annotation(object), "package not available"))
            get(annotation(object))@getdb()
        })

setMethod("annotation", "HmmOptions", function(object) annotation(object@snpset))
setReplaceMethod("annotation", "HmmOptions", function(object, value){
  object@annotation <- value
  object
})

##setMethod("cn.robustSE", "HmmOptions", function(object) object@cn.robustSE)
##setMethod("notes", "HmmOptions", function(object) object@notes)
setReplaceMethod("states", c("HmmOptions", "character"),
                 function(object, value){
                   object@states <- states
                   object
                 })

##setMethod("initialize", "HmmOptions",
##          function(.Object, snpset,
##                   states=character(),
##                   cn.location=numeric(),
##                   cn.SE=numeric(),
##                   cn.robustSE=numeric(),
##                   cn.ICE=FALSE,
##                   gt.gte=numeric(),
##                   gte.state=numeric(),
##                   gt.confidence.states=character(),
##                   gt.ICE=FALSE,
##                   annotation=character(),
##                   beta=matrix(),
##                   notes=""){
##            snpset <- snpset[!(is.na(chromosome(snpset))), ]
##            if(length(states) == 0){
##              states <- switch(class(snpset),
##                               SnpCallSet=c("L", "N"),
##                               SnpCopyNumberSet=c("D", "N", "A"),
##                               oligoSnpSet=c("D", "N", "L", "A"),
##                               stop("if hiddenStates not given,  snpset must be one of the following classes: SnpCallSet, SnpCopyNumberSet, oligoSnpSet"))
##            }
##            if(!("N" %in% states)) stop("one hidden state must have the name 'N' (normal)")          
##
##            if(length(beta) == 1){
##            if(length(cn.location) == 0 & class(snpset) != "SnpCallSet"){
##              cn.location <- switch(class(snpset),
##                                    SnpCopyNumberSet=c(1, 2, 3),
##                                    oligoSnpSet=c(1, 2, 2, 3),
##                                    stop("if cn.location not given, snpset must be one of the following classes: SnpCopyNumberSet, oligoSnpSet"))
##              if(length(cn.location) != length(states)) stop("Must specify cn.location and cn.SE")
##              names(cn.location) <- states
##              notes <- "Assume CNE have been transformed to log2 scale and that the CNE are approx. Gaussian"
##            }           
##            if(class(snpset) == "SnpCallSet"){
##              cn.ICE <- FALSE
##            }
##            if(class(snpset) == "SnpCopyNumberSet"){
##              gt.ICE <- FALSE
##            }
##            if(cn.ICE) {
##              ##print("use inverse of standard errors in the cnConfidence slot")
##              cn.robustSE <- "SNP-specific"
##              notes <- paste(notes, "Values in slot cnConfidence are assumed to be inverse standard errors", sep="; ")
##            } else{
##              if(class(snpset) != "SnpCallSet" && length(cn.robustSE) == 0){
##                autosomes <- chromosome(snpset) != "X" & chromosome(snpset) != "Y" & chromosome(snpset) != "XY"
##                snpset.autosomes <- snpset[autosomes, ]
##		f <- function(x){
##			diff(quantile(x, probs=c(0.16, 0.84), na.rm=TRUE))/2 ##+ quantile(x, probs=0.16)
##		}
##                cn.robustSE <- apply(copyNumber(snpset.autosomes), 2, f)
##                names(cn.robustSE) <- sampleNames(snpset)
##              }
##            }
##            if(gt.ICE){
##              require(callsConfidence) || warning("callsConfidence package is not available")
##              notes <- paste(notes, "***ICE for genotype calls is only developed for the 100k and 500k Affymetrix platforms.  We assume the chip has been pre-processed using oligo.  If not, set gt.ICE=FALSE***", sep="; ")              
##              gt.gte <- c(0.001, 0.001)
##              gt.confidence.states <- switch(class(snpset),
##                                             SnpCallSet=c("L", "N"),
##                                             oligoSnpSet=c("L", "N", "L", "N"),
##                                             stop("if phat.states not supplied, snpset must be a SnpCallSet or an oligoSnpSset"))
##              notes <- paste(notes, "gt.gte refers to the P(true gt | called gt) -- only probabilites for the two types of errors are provided; when calculating emission probabilities using ICE for genotype emission probabilities, we calculate P(confidence score | normal) and P(confidence score | loss.  gt.phat.state tells which states have normal AA/BB : AB and whih have loss (high AA/BB : AB). ", sep="; ")
##              names(gt.gte) <- c("P(gt=HET | gte=HOM)", "P(gt=HOM| gte=HET)")
##
##            } else{
##              gt.gte <- "not used if gt.ICE=FALSE"
##            }
##            if(class(snpset) == "SnpCallSet" | class(snpset) == "oligoSnpSet"){
####            if(class(snpset) != "TrioSet" & class(snpset) != "SnpCopyNumberSet"){
##              if(length(gte.state) == 0){
##                gte.state <- switch(class(snpset),
##                                    SnpCallSet=c(0.999, 0.75),
##                                    oligoSnpSet=c(0.999, 0.75, 0.999, 0.75),
##                                    stop("if gt.Prob not supplied, snpset must be one of the following classes: SnpCallSet, oligoSnpSet"))
##              }
##              names(gte.state) <- states
##              notes <- paste(notes, "gt.state is P(gte = HOM | state)", sep=", ")
##            }
##          }
##            .Object@SnpClass <- class(snpset)
##            .Object@states <- states
##            .Object@cn.location <- cn.location
##            .Object@gte.state <- gte.state
##            .Object@gt.gte <- gt.gte
##            .Object@cn.robustSE <- cn.robustSE
##            .Object@cn.ICE <- cn.ICE
##            .Object@gt.ICE <- gt.ICE
##            .Object@gt.confidence.states <- gt.confidence.states
##            .Object@annotation <- annotation(snpset)
##            .Object@notes <- notes
##            .Object
##          })

##setMethod("calculateDistance", "HmmOptions",
##          function(object){
##		  snpset <- object@snpset
##		  d <- diff(position(snpset))
##		  d
##          })
##setMethod("calculateDistance", "integer",
##	  function(object){
##		  (object[2:length(object)] - object[1:(length(object)-1)])##/(100*1e6)
##})
setMethod("physicalDistance", "HmmOptions", function(object) diff(position(snpset(object))))

##setMethod("snpdata", "HmmOptions", function(object) object@snpset)


setMethod("show", "HmmOptions", function(object){
	show(object@snpset)
	cat("\n hidden states: ", states(object))
	cat("\n copyNumber.location: ", object@copyNumber.location)
	cat("\n copyNumber.scale: ", object@copyNumber.scale)
	cat("\n copyNumber.ICE: ", object@copyNumber.ICE)
	cat("\n calls.ICE: ", object@calls.ICE)
	cat("\n probability of a Homozygous Call: ", object@probHomCall)
	cat("\n term5: ", object@term5, "\n")
})

setMethod("[", "HmmOptions",
          function(x, i, j, ..., drop = FALSE){
            if (missing(drop)) drop <- FALSE
            if (missing(i) && missing(j)) {
              if (length(list(...))!=0)
                stop("specify genes or samples to subset; use '",
                     substitute(x), "$", names(list(...))[[1]],
                     "' to access phenoData variables")
              return(x)
            }
            ##number of rows (SNPs)
            R <- nrow(x@snpset)
            ##number of states 
            S <- length(states(x))
            ##number of columns (samples)
            C <- ncol(x@snpset)
            if(!missing(i) & !missing(j)){
		    x@snpset <- x@snpset[i, j]
            }
            if(!missing(i) & missing(j)){
		    x@snpset <- x@snpset[i, ]
            }
	    if(missing(i) & !missing(j)){
		    x@snpset <- x@snpset[, j]
	    }
            x
          })

setMethod("calls.ICE", "HmmOptions", function(object) object@calls.ICE)
setMethod("copyNumber.ICE", "HmmOptions", function(object) object@copyNumber.ICE)
setReplaceMethod("copyNumber.ICE", c("HmmOptions", "logical"),
                 function(object, value){
			 object@copyNumber.ICE <- value
			 object
		 })
setReplaceMethod("calls.ICE", c("HmmOptions", "logical"),
                 function(object, value){
			 object@calls.ICE <- value
			 if(length(object@term5) != 2){
				 object@term5 <- c(0.999, 0.001)
				 print("using default values for term 5 slot. See HmmOptions help")
			 }
			 if(length(object@probHomCall) != 2){
				 object@probHomCall <- c(0.99, 0.7)
				 print("only probabilities for states LOH and Normal are allowed when calls.ICE is TRUE. Using defaults")				 
			 }
			 object
                 })

setMethod("probHomCall", "HmmOptions", function(object) object@probHomCall)
setReplaceMethod("probHomCall", signature(object="HmmOptions", value="numeric"),
		 function(object, value){
			 object@probHomCall <- value
			 object
		 })
setMethod("snpset", "HmmOptions", function(object) object@snpset)
setReplaceMethod("snpset", signature(object="HmmOptions", value="SnpLevelSet"),
		 function(object, value){
			 object@snpset <- value
			 object
		 })
setMethod("term5", "HmmOptions", function(object) object@term5)
setReplaceMethod("term5", signature(object="HmmOptions", value="numeric"),
		 function(object, value){
			 object@term5 <- value
			 object
		 })

##setMethod("SnpClass", "HmmOptions", function(object) object@SnpClass)
setMethod("states", "HmmOptions", function(object) object@states)
setMethod("copyNumber.location", "HmmOptions", function(object) object@copyNumber.location)
setReplaceMethod("copyNumber.location", signature(object="HmmOptions", value="numeric"),
		 function(object, value){
			 object@copyNumber.location <- value
			 object
		 })

setMethod("calculateDistance", "HmmOptions",
	  function(object){
		  calculateDistance(snpset(object))
	  })
##setReplaceMethod("states", "HmmOptions", function(object, value){
##  object@states <- value
##  object
##})
