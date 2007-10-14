setMethod(".unsplitHmm", c("list", "AnnotatedDataFrame"),
          function(from, annotatedDataFrame, ...){
            object <- switch(class(from[[1]]),
                             HmmSnpSet=new("HmmSnpSet", ...),
                             HmmSnpCallSet=new("HmmSnpCallSet", ...),
                             HmmSnpCopyNumberSet=new("HmmSnpCopyNumberSet", ...),
                             stop("Must be a class defined in VanillaICE"))
            featureData(object) <- annotatedDataFrame[match(featureNames(object), featureNames(annotatedDataFrame)), , ]
            stopifnot(validObject(object))
            object
          })
