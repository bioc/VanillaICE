setMethod("initialize", "ParHmmSnpCallSet",
          function(.Object, ...){
            .Object <- callNextMethod()
            if(.Object$add.cytoband){
              .Object$heights <- c(1, 1, 0.8)
              } else {
                .Object$heights <- c(1, 1)
              }
            .Object$outer.ylab <- FALSE
            .Object$ylab <- "genotype call"
            .Object$col <- c("lightblue", "red")
            snpPar(.Object)$col.predict <- c("yellow", "white")
            .Object
          })
