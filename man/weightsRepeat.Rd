\name{weightsRepeat}
\alias{weightsRepeat}

\docType{data}

\title{
  Weights for repeat region bias
}

\description{
  The default weights of all supported repeat regions for all supported SV formation mechanisms.
  It is loaded automatically when the repeat bias feature has been turned on (\code{repeatBias=TRUE} in function \code{simulateSV}).
  (But this dataset is intended for internal use; no need to handle it manually.)
}

\usage{
  data(weightsRepeats)
}

\format{
    Formal class 'data.frame'
}


\examples{
  data(weightsRepeats)
  weightsRepeats
}

\keyword{datasets}
