setMethod("plotSnp", c("ParHSet", "SnpLevelSet"), ##"hSet",
          function(object, snpset, hmmPredict, ...){
            old.par <- par(no.readonly=TRUE)
            on.exit(par(old.par))
	    if(!is.null(object$abline.v)){
		    if(length(sampleNames(snpset)) == 1){
			    v <- breakpoints(hmmPredict)[, c("start", "last")]*1e6
			    v <- v[v != min(v) & v != max(v)]
			    object$abline.v <- v
		    }
	    }
            callNextMethod(object, snpset)
            predictions <- data.frame(predictions(hmmPredict))
            breaks <- breakpoints(hmmPredict)
            ##plots all chromosomes in a sample
	    breaks <- split(breaks, breaks[, "id"])
	    idx <- match(sampleNames(snpset), names(breaks))
	    breaks <- breaks[idx]
            mapply(plotPredictions,
                   predictions=predictions,
                   breakpoints=breaks,
                   MoreArgs=list(op=object,
                     position=position(snpset),
                     chromosome=chromosome(snpset),
                     states=states(hmmPredict)),
                   SIMPLIFY=FALSE)
	    return()
          })

setMethod("plotSnp", c("ParHmmSnpSet", "SnpLevelSet"), ##"hSet",
          function(object, snpset, hmmPredict, ...){
            idx <- match(featureNames(hmmPredict), featureNames(snpset))
            idx <- idx[!is.na(idx)]
            snpset <- snpset[idx, ]
            idx <- match(featureNames(snpset), featureNames(hmmPredict))
            hmmPredict <- hmmPredict[idx, ]
            if(!identical(featureNames(snpset), featureNames(hmmPredict)))
              stop("featureNames in snpset and hmmPredict arguments are not identical")
            callNextMethod(object, snpset, hmmPredict, hmmParams, ...)
	    return()
          })

setMethod("plotSnp", c("ParHmmSnpCopyNumberSet", "SnpLevelSet"), ##"hSet",
          function(object, snpset, hmmPredict, ...){
            idx <- match(featureNames(hmmPredict), featureNames(snpset))
            idx <- idx[!is.na(idx)]
            snpset <- snpset[idx, ]
            idx <- match(featureNames(snpset), featureNames(hmmPredict))
            hmmPredict <- hmmPredict[idx, ]
            if(!identical(featureNames(snpset), featureNames(hmmPredict)))
              stop("featureNames in snpset and hmmPredict arguments are not identical")
            callNextMethod(object, snpset, hmmPredict, hmmParams, ...)
          })

setMethod("plotSnp", c("ParHmmSnpCallSet", "SnpLevelSet"), ##"hSet",
          function(object, snpset, hmmPredict, ...){
            idx <- match(featureNames(hmmPredict), featureNames(snpset))
            idx <- idx[!is.na(idx)]
            snpset <- snpset[idx, ]
            idx <- match(featureNames(snpset), featureNames(hmmPredict))
            hmmPredict <- hmmPredict[idx, ]
            if(!identical(featureNames(snpset), featureNames(hmmPredict)))
              stop("featureNames in snpset and hmmPredict arguments are not identical")
            callNextMethod(object, snpset, hmmPredict, hmmParams, ...)            
          })

plotPredictions <- function(predictions,
                            op,
                            breakpoints,
                            position,
                            chromosome,
                            states, ...){
	chr <- unique(chromosome)
	for(i in chr){
	  tmp <- breakpoints[breakpoints[, "chr"] == i, , drop=FALSE]
	  plot(0:1,
	       type="n",
	       xlim=op$xlim[i, ],
	       ylim=op$predictions.ylim,
	       yaxt="n",
	       xaxt="n",
	       xaxs=op$xaxs,
	       yaxs="i",
	       xaxs=op$xaxs)
	  if(i == op$firstChromosome){
		  mtext("HMM", side=2, at=0.5, cex=op$cex.axis, las=1)
	  }
	  .drawRect <- function(x, object, position, op){
		  col <- op$col.predict
		  if(length(x) < 1) return()
		  start <- max(as.numeric(x["start"]) * 1e6, op$xlim[1])
		  last <- min(as.numeric(x["last"]) * 1e6, op$xlim[2])
		  predict <- predictions[position >= start & position <= last & chromosome == x["chr"]]
		  predict <- predict[!is.na(predict)]
		  if(length(unique(predict)) > 1) {
			  browser()
			  stop("predictions not unique")
		  }
		  col <- col[unique(predict)]
		  rect(xleft=start,
		       ybottom=0,
		       xright=last,
		       ytop=1,
		       col=col,
		       border=col)
	  }
	  if(!is.null(tmp)){
		  ##best to draw the biggest regions first and the smallest regions last.
		  tmp <- tmp[order(tmp[, "size"], decreasing=TRUE), ]
		  apply(tmp, 1, .drawRect, position=position, op=op)
	  }
	  ##white rectangle at centromere
	  rect(xleft=centromere(i)[1],
	       ybottom=0,
	       xright=centromere(i)[2],
	       ytop=1,
	       col="white",
	       border="white")
	  rect(xleft=op$xlim[i, 1],
	       ytop=1,
	       ybottom=0,
	       xright=op$xlim[i, 2],
	       border=op$border.prediction)
	  if(op$legend.prediction){
		  legend(op$legend.location.prediction,
			 fill=op$col.predict[unique(predictions)],
			 legend=states[unique(predictions)],
			 ncol=length(states[unique(predictions)]),
			 bty="n",
			 cex=op$cex.legend)
	  }
  }
}
setMethod("show", "ParHSet", function(object) str(snpPar(object)))




##setMethod("plotSnp", c("ParHmmSnpCopyNumberSet", "HmmSnpCopyNumberSet"), ##"hSet",
##          function(object, snpset, breaks, ...){
##            callNextMethod()
##          })
##

