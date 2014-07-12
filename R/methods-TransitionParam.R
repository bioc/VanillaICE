setMethod(TransitionParam, signature(taup="missing"),
          function(taup=1e10, taumax=1-5e-8){
            new("TransitionParam", taup=taup, taumax=taumax)
          })

setMethod(TransitionParam, signature(taup="numeric"),
          function(taup=1e10, taumax=1-5e-8){
            new("TransitionParam", taup=taup, taumax=taumax)
          })

setMethod(show, "TransitionParam",
          function(object) {
            cat(class(object), ":\n")
            cat("TAUP: ", taup(object), "\n")
            cat("max tau: ", taumax(object), "\n")
          })

setMethod("calculateTransitionProbability", signature(x="numeric"),
          function(x, param=TransitionParam()){
            if(length(x) == 1){
              ## assume a distance between markers is provided
              d <- x
            } else d <- diff(x)
            TAUP <- taup(param)
            TAUMAX <- taumax(param)
            probability <- exp(-2*d/TAUP)
            if(TAUMAX < 1){
              probability[probability > TAUMAX] <- TAUMAX
            }
            ##if(min(probability) < 0) stop("Invalid transition probabilities. Check that x values are sorted")
            ##if(max(probability) > 1) stop("Invalid transition probabilities. Check that x values are sorted")
            probability
          })
