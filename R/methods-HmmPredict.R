setMethod("breakpoints", "HmmPredict", function(object) object@breakpoints)
setReplaceMethod("breakpoints", c("HmmPredict", "ANY"),
		 function(object, value){
			 object@breakpoints <- as.data.frame(value)
			 object
		 })
setMethod("states", "HmmPredict", function(object) object@states)
setMethod("SnpClass", "HmmPredict", function(object) object@SnpClass)

setMethod("show", "HmmPredict", function(object){
	##	cat("Instance of class HmmPredict \n")	
	callNextMethod(object)
	cat("\n SnpClass:", SnpClass(object), "\n")
	cat("hidden states: ", states(object), "\n")
	cat("breakpoints: \n")
	str(breakpoints(object))
})


setMethod("predictions", "HmmPredict", function(object) assayDataElement(object, "predictions"))
setReplaceMethod("predictions", signature(object="HmmPredict", value="matrix"),
                 function(object, value) assayDataElementReplace(object, "HmmPredict", value))

setMethod("[", "HmmPredict",
          function(x, i, j, ..., drop = FALSE){
		  object <- callNextMethod()
		  if(!missing(i) && !missing(j)){
			  chr <- unique(chromosome(object))
			  f <- function(x, chr) x[x$chr %in% chr, ]
			  id <- sampleNames(object)[j]
			  breakpoints(object) <- breakpoints(object)[breakpoints(object)[, "id"] %in% id, drop=FALSE]
		  }            
		  if(missing(i) && !missing(j)){
			  id <- sampleNames(object)[j]
			  breakpoints(object) <- breakpoints(object)[breakpoints(object)[, "id"] %in% id]			  
		  }
		  object
          })

setMethod("calculateBreakpoints", "HmmPredict",
	  function(object, ...){
		  breaks <- list()
		  for(i in 1:ncol(object)){
			  x <- split(predictions(object)[, i], chromosome(object))
			  breaks[[i]] <- mapply(.calculatebreaks,
						x=x,
						chromosome=as.list(names(x)),
						MoreArgs=list(position=position(object),
						states=states(object),
						sampleNames=sampleNames(object)[i]),
						SIMPLIFY=FALSE)
			  breaks[[i]] <- do.call("rbind", breaks[[i]])
		  }
		  breaks <- do.call("rbind", breaks)
		  return(as.data.frame(breaks))
	  })

setMethod("summary", "HmmPredict",
	  function(object, by="chr", distance.unit=1000, digits=3, ...){
		  states <- unique(breakpoints(object)[, "state"])
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
                   
