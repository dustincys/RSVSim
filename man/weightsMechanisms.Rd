\name{weightsMechanisms}
\alias{weightsMechanisms}

\docType{data}

\title{
  Weights for SV formation mechanisms
}

\description{
  The default weights of all supported SV formation mechanisms for all SVs.
  It is loaded automatically when the repeat bias feature has been turned on (\code{repeatBias=TRUE} in function \code{simulateSV}).
  (But this dataset is intended for internal use; no need to handle it manually.)
}

\usage{
  data("weightsMechanisms")
}

\format{
    Formal class 'data.frame'
}


\examples{
  data("weightsMechanisms")
  weightsMechanisms
}

\keyword{datasets}
