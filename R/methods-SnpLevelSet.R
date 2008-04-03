setMethod("calculateDistance", "SnpLevelSet",
          function(object){
		  d <- (position(object)[2:nrow(object)] - position(object)[1:(nrow(object)-1)])##/(100*1e6)
          })


setMethod("calculateTransitionProbability", c("SnpLevelSet", "HmmOptions"),
          function(object, options, scale, ...){
            ##distance in units of 100 Mb ("~ 1cM") as in Beroukeim 2006
            states <- options@states
            S <- length(states)
            d <- (position(object)[2:nrow(object)] - position(object)[1:(nrow(object)-1)])/(100*1e6)
            

            ##SNP-specific transition probabilities.  We will assume that the
            ##probability of remaining in the same state is the same for all
            ##states.
##            exp(-2*d)  #probability of staying in the same state. (this is 1-theta)
          
            theta <- 1-exp(-2*d)

            tmp <- matrix(1:(S^2), nrow=S, ncol=S)
            rownames(tmp) <- colnames(tmp) <- states
            columns <- diag(tmp)

            tNames <- list()
            for(i in 1:S){
              tNames[[i]] <- paste(states[i], states, sep="->")
            }
            tNames <- do.call("rbind", tNames)
            tNames <- as.vector(tNames)
            
            ##The probability of staying in the same state is 1 - theta
            A <- array(NA, c(nrow(object)-1, S, S))
            
            A <- matrix(NA, nrow=(nrow(object) - 1), ncol=S^2)
            colnames(A) <- tNames
            A[, columns] <- 1-theta ## should recycle appropriately
            epsilon <- theta/(S-1) ##or 1-(1-theta)
            A[, -columns] <- epsilon  ##The rowSums should be S (the number of states)
  
            ##Scale the probability of transitioning to a state that
            ##is not normal
            tNames <- matrix(tNames, S, S)
            t.vector <- as.vector(tNames)
            ##Scale the transition probabilities by taking the product
            ##of the scale vector and each row in the transition
            ##probability matrix.  Easy to do this by transposing the
            ##A matrix (columns in the transposed matrix of length
            ##S^2), then multiplying by the scale vector will
            ##automatically multiply each column by scale.  Transpose
            ##the product to get the desired result
            scale <- as.vector(scale)                          
            tau <- t(t(A) * scale)

            ##may not add up to 1 for snps on different chromosomes -- distance may be a negative number
            idx <- tau <= 0
            idx <- rowSums(idx) > 0
            if(any(idx)) {
              tau[idx, columns] <- 0.9999
              tau[idx, -columns] <- (1-0.9999)/(S - 1)
            }
            log(tau)
          })





