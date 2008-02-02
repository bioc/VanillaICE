##Default:  ordinary visualization of data with plotSnp for one sample / one chromosome
##Add a row for the prediction bar
##Add a row for the cytoband
setMethod("initialize", "ParHSet",
          function(.Object,
                   border.prediction=grey(0.4),
                   legend.prediction=TRUE,
                   legend.location.prediction="topleft", ...){
            .Object <- callNextMethod(.Object, ...)
            if(.Object$useLayout){
              if(.Object$add.cytoband){
                .Object$heights <- c(1, 0.1, 0.05)
                ##1 - plot the data
                ##2 - plot the cytoband
                ##3 - plot the predictions.  Plot the predictions
                ##between the data and the cytoband
                .Object$mat <- matrix(c(1, 3, 2), ncol=1)                
              } else {
                ##do not plot the cytoband
                .Object$heights <- c(1, 0.2)
                .Object$mat <- matrix(1:2, ncol=1)                                
              }
            }
            snpPar(.Object)$border.prediction <- border.prediction
            snpPar(.Object)$legend.prediction <- legend.prediction
            snpPar(.Object)$legend.location.prediction <- legend.location.prediction
            snpPar(.Object)$predictions.ylim <- c(-0.2, 1.8)
            .Object
          })

setMethod("initialize", "ParHmmSnpSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$col <- c("lightblue", "red", "lightblue", "green3")
            .Object$bg <- c("lightblue", "red", "lightblue", "green3")
            .Object$ylab <- "copy number"
            ##add col.predict as an argument
            require(RColorBrewer, quietly=TRUE) || stop("RColorBrewer package not available")
            col.predict <- brewer.pal(7, "BrBG")[c(7, 2)]
            snpPar(.Object)$col.predict <- c(col.predict[1], "white", "yellow", col.predict[2])
            .Object
          })

setMethod("initialize", "ParHmmSnpCopyNumberSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$ylab <- "copy number"
            require(RColorBrewer, quietly=TRUE) || stop("RColorBrewer package not available")
            col.predict <- brewer.pal(7, "BrBG")[c(7, 2)]
            snpPar(.Object)$col.predict <- c(col.predict[1], "white", col.predict[2])               
            .Object
          })

setMethod("initialize", "ParHmmSnpCallSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            if(.Object$add.cytoband){
              .Object$heights <- c(1, 1, 0.8)
              } else {
                .Object$heights <- c(1, 1)
              }
            .Object$outer.ylab <- FALSE
            .Object$ylab <- "genotype call"
            .Object$col <- c("lightblue", "red", "lightblue")
            snpPar(.Object)$col.predict <- c("black", "white")
            .Object
          })

setMethod("getPar", "HmmPredict",
          function(object, snpset, ...){
            op <- switch(SnpClass(object),
                         oligoSnpSet=new("ParHmmSnpSet", ...),
                         SnpCallSet=new("ParHmmSnpCallSet", ...),
                         SnpCopyNumberSet=new("ParHmmSnpCopyNumberSet", ...),
                         stop("Object must be oligoSnpSet, SnpCallSet, or SnpCopyNumberSet"))
            snpset <- snpset[!is.na(chromosome(snpset)), ]
            chromosomeNames <- unique(chromosome(snpset))
            chromosomeNames <- chromosomeNames[order(as.numeric(chromosomeNames))]
            op$firstChromosome <- chromosomeNames[1]
            N <- length(chromosomeNames)
            S <- ncol(snpset)

            if(op$add.cytoband){
              ## plots: data, predictions, data, predictions, ... data, predictions, cytoband
              ## == S * 2 + 1
              op$heights <- c(rep(c(1, 0.2), S), 0.1)
              S <- S+1
            } else op$heights <- rep(c(1, 0.2), S)
            data(chromosomeAnnotation, package="SNPchip", envir=environment())
            w <- chromosomeAnnotation[chromosomeNames, "chromosomeSize"]
            op$widths <- w/min(w)

            R <- length(op$heights)
            C <- length(op$widths)
            nPlots <- R*C
            if(N*S > 1){
              ##plot all the data and cytobands, then plot predictions
              ##the data
              m <- matrix(1:(N*S), ncol=N, nrow=S, byrow=FALSE)
              ##Add rows for predictions
              if(op$add.cytoband) p <- matrix(1:(N*ncol(snpset)), ncol=N, nrow=ncol(snpset), byrow=FALSE)
              p <- p+max(m)
              ##combine the matrices
              mat <- matrix(NA, R, C)
              k <- 1; i <- 1
              while(i < R){
                mat[i, ] <- m[k, ]
                mat[i+1, ] <- p[k, ]
                k <- k+1
                i <- i+2      
              }
              mat[R, ] <- m[nrow(m), ]
              op$mat <- mat
            }
            op$ylim <- .calculateYlim(object=snpset, op=op)
            ##if plotting multiple chromosomes, xlim is not set
            ##Might be better to make a matrix, with rownames as the chromosomes
            if(op$use.chromosome.size){
              op$xlim <- matrix(NA, nrow=length(chromosomeNames), ncol=2)    
              op$xlim[, 1] <- rep(0, nrow(op$xlim))
              op$xlim[, 2] <- chromosomeSize(chromosomeNames)
              rownames(op$xlim) <- chromosomeNames    
            } else{
              objList <- split(snpset, chromosome(snpset))
              objList <- objList[chromosomeNames]
              op$xlim <- t(sapply(objList, function(snpset) range(position(snpset))))
            }  
            if(ncol(snpset) == 1) op$yaxt="s"
            op
          })



