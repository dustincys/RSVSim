\name{estimateSVSizes}

\alias{estimateSVSizes}
\alias{estimateSVSizes,numeric,missing,ANY,ANY,character-method}
\alias{estimateSVSizes,numeric,missing,ANY,ANY,missing-method}
\alias{estimateSVSizes,numeric,numeric,ANY,ANY,missing-method}


\title{Draw random structural variation sizes from a beta distribution}

\description{\pkg{RSVSim} can implement structural variations (SVs) of specific sizes. \command{estimateSVSizes} draws random values for SV sizes from a beta distribution. It can fit the distribution according to given SV sizes or default values.}

\usage{estimateSVSizes(n, svSizes, minSize, maxSize, default, hist=TRUE)}

\arguments{
 \item{n}{The number of SVs to simulate}
 \item{svSizes}{A numeric vector with SV sizes to calculate the parameters for the beta distribution}
 \item{minSize}{The minimum returned SV size}
 \item{maxSize}{The maximum returned SV size}
 \item{hist}{TRUE or FALSE to show a histogram of the SV sizes (or not)}
 \item{default}{Setting this to "deletions", "insertions", "inversions" or "tandemDuplications" loads default shape parameters for the beta distribution for these SV types (see details below)}
}

\details{
  \itemize{
    \item \code{minSize} and \code{maxSize} are optional and taken from the given set of \code{svSizes} if omitted
    \item The default shape parameters for deletions, insertions, inversions and tandem duplications were estimated from sequencing studies in the Database of Genomic Variants release 2012-03-29. In total, 1.129 deletions, 490 insertions, 202 inversions and 145 tandem duplications between 500bp and 10kb were used to estimate the shape of the distribution. 
  }
}

\value{A numeric vector with beta distributed values between \code{minSize} and \code{maxSize}.}

\references{Database of Genomic Variants: \url{http://dgvbeta.tcag.ca/dgv/app/home?ref=NCBI36/hg18}}

\note{
  \itemize{
    \item It is intended to run this function separately for every SV type and then provide the output to the respective size argument of \code{\link{simulateSV}}
    \item It is recommended to use a \code{minSize} and \code{maxSize} that is consistent with the minimum/maximum values in \code{svSizes}.
    \item When using default shape parameters for the beta distribution, it works best to simulate SVs that do not differ too much in size (aroung 500bp-10kb).
  }
}

\author{Christoph Bartenhagen}

\examples{
  ## estimate sizes for 20 SVs from a given set of values:
  svSizes = c(10,20,30,40,60,80,100,150,200,250,300,400,500,750,1000)
  estimateSVSizes(n=20, svSizes=svSizes, hist=FALSE) 

  ## when using the default shape parameters for deletions:
  estimateSVSizes(n=20, minSize=500, maxSize=10000, hist=FALSE, default="deletions")
}
