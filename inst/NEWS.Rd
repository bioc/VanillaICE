\name{NEWS}
\title{News for Package 'VanillaICE'}

\section{Changes in version 1.18.0}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Using NEWS.Rd

      \item oligoSnpSet, CopyNumberSet, and BeadStudioSet (defined in
      oligoClasses) require that all assayData elements are integer
      matrices. In particular, copy number should be scaled by 100 and
      coerced to integers and B allele frequencies should be scaled by
      1000 and coerced to integers.  The function integerMatrix(x, 100)
      in the oligoClasses package will coerce a matrix of doubles, 'x',
      to a matrix of integers scaled by a factor of 100.  The hmm
      methods for these classes expect that the assay data elements are
      appropriately scaled integer matrices.

    }
  }
}