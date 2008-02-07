setMethod("initialize", "HmmPredict",
          function(.Object,
                   states=character(),
                   predictions=matrix(),
                   breakpoints=list(),
                   SnpClass=character(),
                   featureData=new("AnnotatedDataFrame")){
            .Object@predictions <- predictions
            .Object@breakpoints <- breakpoints
            .Object@SnpClass <- SnpClass
            .Object@states <- states
            .Object@featureData <- featureData
            .Object
          })

setMethod("featureNames", "HmmPredict", function(object) rownames(predictions(object)))
setMethod("featureData", "HmmPredict", function(object) object@featureData)
setMethod("fData", "HmmPredict", function(object) pData(object@featureData))
setMethod("sampleNames", "HmmPredict", function(object) colnames(predictions(object)))
setMethod("predictions", "HmmPredict", function(object) object@predictions)
setMethod("position", "HmmPredict", function(object) fData(object)$position)
setMethod("chromosome", "HmmPredict", function(object) fData(object)$chromosome)
setMethod("breakpoints", "HmmPredict", function(object) object@breakpoints)
setMethod("states", "HmmPredict", function(object) object@states)
setMethod("SnpClass", "HmmPredict", function(object) object@SnpClass)

setMethod("show", "HmmPredict", function(object){
  cat("Instance of class HmmPredict \n")
  print("predictions:")
  str(predictions(object))
  cat("\n breakpoints: \n")
  breaks <- breakpoints(object)  
  i <- match(c("size", "start", "last"), colnames(breaks[[1]]))
  breaks <- lapply(breaks, function(object, i){ object[, i] <- round(object[, i], 3); object}, i=i)
  tmp <- lapply(breaks, selectSomeIndex)
  idx <- lapply(tmp, unlist)
  print(mapply(function(object, idx) object[idx, ], object=breaks, idx=idx, SIMPLIFY=FALSE))
  breaks
  cat("\n featureData \n")
  show(featureData(object))
})



setMethod("[", "HmmPredict",
          function(x, i, j, ..., drop = FALSE){
            if (missing(drop)) drop <- FALSE
            if (missing(i) && missing(j)) {
              if (length(list(...))!=0)
                stop("specify genes or samples to subset; use '",
                     substitute(x), "$", names(list(...))[[1]],
                     "' to access phenoData variables")
              return(x)
            }
            if(!missing(i) && missing(j)){
              x@predictions <- x@predictions[i, , drop=FALSE]
              x@featureData <- x@featureData[i, , drop=FALSE]
            }
            if(missing(i) && !missing(j)){
              x@predictions <- x@predictions[, j, drop=FALSE]
            }            
            if(!missing(i) && !missing(j)){
              x@predictions <- x@predictions[i, j, drop=FALSE]
              x@featureData <- x@featureData[i, , drop=FALSE]
            }
            x
          })

setMethod("selectSomeIndex",
          signature(object="data.frame"),
          function(object, maxToShow=5, byrow=TRUE, ...) {
              len <-
                if (byrow) dim(object)[[1]]
                else dim(object)[[2]]
              if (maxToShow < 3) maxToShow <- 3
              if (len > maxToShow) {
                  maxToShow <- maxToShow - 1
                  bot <- ceiling(maxToShow/2)
                  top <- len-(maxToShow-bot-1)
                  list(1:bot, "...", top:len)
              } else if (len >= 1) list(1:len, NULL, NULL)
              else list(NULL, NULL, NULL)
          })
                   
