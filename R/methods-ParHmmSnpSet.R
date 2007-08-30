setMethod("initialize", "ParHmmSnpSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$col <- c("lightblue", "red", "lightblue")
            .Object$bg <- c("lightblue", "red", "lightblue")            
            .Object$ylab <- "copy number"
            ##add col.predict as an argument
            require(RColorBrewer, quietly=TRUE) || stop("RColorBrewer package not available")
            col.predict <- brewer.pal(7, "BrBG")[c(7, 2)]
            snpPar(.Object)$col.predict <- c(col.predict[1], "white", "yellow", col.predict[2])
            .Object
          })
