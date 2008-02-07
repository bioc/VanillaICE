setMethod("getPar", c("ParHSet", "SnpLevelSet"),
          function(object, snpset, ...){
            object <- callNextMethod(object, snpset, ...)
            
            ##only need to change the layout
            N <- length(unique(chromosome(snpset)))
            S <- ncol(snpset)
            if(object$add.cytoband){
              ## plots: data, predictions, data, predictions, ... data, predictions, cytoband
              ## == S * 2 + 1
              object$heights <- c(rep(c(1, 0.2), S), 0.1)
              S <- S+1
            } else object$heights <- rep(c(1, 0.2), S)
            R <- length(object$heights)
            C <- length(object$widths)
            nPlots <- R*C
            if(N*S > 1){
              ##plot all the data and cytobands, then plot predictions
              ##the data
              m <- matrix(1:(N*S), ncol=N, nrow=S, byrow=FALSE)
              ##Add rows for predictions
              if(object$add.cytoband) p <- matrix(1:(N*ncol(snpset)), ncol=N, nrow=ncol(snpset), byrow=FALSE)
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
              object$mat <- mat
            }
            if(ncol(snpset) == 1) object$yaxt="s"
            object
          })
