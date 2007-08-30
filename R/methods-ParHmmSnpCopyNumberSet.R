setMethod("initialize", "ParHmmSnpCopyNumberSet",
          function(.Object, ...){
            .Object <- callNextMethod()
            .Object$ylab <- "copy number"
            require(RColorBrewer, quietly=TRUE) || stop("RColorBrewer package not available")
            col.predict <- brewer.pal(7, "BrBG")[c(7, 2)]
            snpPar(.Object)$col.predict <- c(col.predict[1], "white", col.predict[2])            
            .Object
          })
