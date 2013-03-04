\name{segDups}
\alias{segDups}

\docType{data}

\title{
  Segmental duplications
}

\description{
  A list with coordinates of segmental duplications detected by RepeatMasker.
  It is loaded automatically when the repeat bias feature has been turned on (\code{repeatBias=TRUE} in function \code{simulateSV}).
  (But this dataset is intended for internal use; no need to handle it manually.)
}

\usage{
  data(segmentalDuplications)
}

\format{
    Formal class 'GRanges'
}

\source{
  Downloaded from the UCSC Browser's RepeatMasker track.
  Meyer et al., The UCSC Genome Browser database: extensions and updates 2013, 2013, Nucleic Acids Res, 41(Database issue), 64-69.
  Smit et al., RepeatMasker Open-3.0., 1996-2010, <http://www.repeatmasker.org>.
}

\examples{
  data(segmentalDuplications)
  segDups
}

\keyword{datasets}
