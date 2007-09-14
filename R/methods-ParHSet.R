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

setMethod("plotSnp", c("ParHSet", "ANY"), ##"hSet",
          function(object, snpset, breaks, ...){
            old.par <- par(no.readonly=TRUE)
            on.exit(par(old.par))
            callNextMethod()            
            plotPredictions(object=snpset, op=object, breaks=breaks)
          })




