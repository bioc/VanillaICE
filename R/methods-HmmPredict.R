setMethod("breakpoints", "HmmPredict", function(object) object@breakpoints)
setReplaceMethod("breakpoints", c("HmmPredict", "ANY"),
		 function(object, value){
			 object@breakpoints <- as.data.frame(value)
			 object
		 })
setMethod("states", "HmmPredict", function(object) object@states)
##setMethod("SnpClass", "HmmPredict", function(object) object@SnpClass)

setMethod("show", "HmmPredict", function(object){
	callNextMethod(object)
##	cat("\n SnpClass:", SnpClass(object), "\n")
	cat("hidden states: ", states(object), "\n")
	cat("breakpoints: \n")
	str(breakpoints(object))
})


setMethod("predictions", "HmmPredict", function(object) assayDataElement(object, "predictions"))
setReplaceMethod("predictions", signature(object="HmmPredict", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "HmmPredict", value))

setMethod("[", "HmmPredict",
          function(x, i, j, ..., drop = FALSE){
		  ##subset the features, but not the samples
		  x <- callNextMethod(x=x, i=i)
		  if(!missing(i) && !missing(j)){
			  ##What chromosomes remain?
			  chr <- unique(chromosome(x))
			  ##Keep the breakpoints in the above chromosomes
			  if(nrow(breakpoints(x)) > 0){
				  breakpoints(x) <- breakpoints(x)[breakpoints(x)$chr %in% chr, , drop=drop]
				  id <- sampleNames(x)[j]
				  breakpoints(x) <- breakpoints(x)[breakpoints(x)$sample %in% id, , drop=drop]
			  }
		  }            
		  if(missing(i) && !missing(j)){
			  id <- sampleNames(x)[j]
			  if(nrow(breakpoints(x)) > 0){
				  breakpoints(x) <- breakpoints(x)[breakpoints(x)$sample %in% id, , drop=drop]
			  }
		  }
		  ##Now subset the i as per use
		  object <- callNextMethod(x, j=j)		  
		  object
          })

setMethod("calculateBreakpoints", "HmmPredict",
	  function(object, ...){
		  regions <- list()
		  for(j in 1:ncol(object)){
			  regions[[j]] <- findBreaks(x=predictions(object)[, j],
						     states=states(object),
						     position=position(object),
						     chromosome=chromosome(object),
						     sample=sampleNames(object)[j])
		  }
		  regions <- do.call("rbind", regions)
		  return(regions)
	  })
#setMethod("findBreaks", "HmmPredict", function(object, ...){
#	  function(object, ...){
#		  regions <- list()
#		  for(j in 1:ncol(object)){
#			  regions[[j]] <- findBreaks(x=predictions(object)[, j],
#						     states=states(object),
#						     position=position(object),
#						     chromosome=chromosome(object),
#						     sample=sampleNames(object)[j])
#		  }
#		  regions <- do.call("rbind", regions)
#		  return(regions)
#	  })	


setMethod("summary", "HmmPredict",
	  function(object, by="chr", distance.unit=1000, digits=3, ...){
		  stop("summary method needs work")
		  states <- unique(calculateBreakpoints(object)[, "state"])
		  states <- states[states != "N"]		  
		  tab <- list()
		  for(i in 1:length(states)){
			  tab[[i]] <- makeTable(object,
						state=states[i],
						by=by,
						unit=distance.unit,
						digits=digits)
		  }
		  names(tab) <- states
		  tab
	  })
	  
##setMethod("selectSomeIndex",
##          signature(object="data.frame"),
##          function(object, maxToShow=5, byrow=TRUE, ...) {
##              len <-
##                if (byrow) dim(object)[[1]]
##                else dim(object)[[2]]
##              if (maxToShow < 3) maxToShow <- 3
##              if (len > maxToShow) {
##                  maxToShow <- maxToShow - 1
##                  bot <- ceiling(maxToShow/2)
##                  top <- len-(maxToShow-bot-1)
##                  list(1:bot, "...", top:len)
##              } else if (len >= 1) list(1:len, NULL, NULL)
##              else list(NULL, NULL, NULL)
##          })
                   
