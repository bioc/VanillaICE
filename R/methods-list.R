setMethod(".unsplitHmm", c("list", "AnnotatedDataFrame"),
          function(from, annotatedDataFrame, ...){
            object <- switch(class(from[[1]]),
                             HmmSnpSet=new("HmmSnpSet", ...),
                             HmmSnpCallSet=new("HmmSnpCallSet", ...),
                             HmmSnpSet=new("HmmSnpSet", ...),
                             stop("Must be a class defined in VanillaICE"))
            featureData(object) <- annotatedDataFrame[match(featureNames(object), featureNames(annotatedDataFrame)), , ]
            stopifnot(identical(rownames(copyNumber(object)), rownames(fData(object))))            
            object
          })
