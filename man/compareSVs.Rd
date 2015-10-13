\name{compareSV}

\alias{compareSV}
\alias{compareSV,character,character-method}
\alias{compareSV,character,data.frame-method}
\alias{compareSV,data.frame,data.frame-method}

\title{Compare the simulation with a set of SVs}

\description{A typical use case of structural variant (SV) simulation with \code{\link{simulateSV}} is the evaluation of SV detection algorithms. The function \code{compareSV} looks for breakpoint overlaps between the output of the simulation (ground truth) and the output of an SV detection program and computes the sensitivity and precision. There is currently no common standard format for SVs. Because the main information about SVs is their position in the genome and, sometimes, the breakpoint sequence (depending on the SV detection algorithm), \code{compareSV} expects the SV detections in a simple BED- or BEDPE format.}

\usage{compareSV(querySVs, simSVs, tol=200)}

\arguments{
  \item{querySVs}{The set detected of SVs. Either a filename for a table of SVs on disk or a \code{data.frame}. The table has to be in BED- or BEDPE-format (see details below). The \code{querySVs} may only contain SVs of one type and not a mixture of several kinds of SVs.}
  \item{simSVs}{The set of simulated SVs as returned by the function \code{\link{simulateSV}}. It can be either a filename for a table of SVs on disk or a \code{data.frame}.}
  \item{tol}{The tolerance in bp, when comparing the (approximate) positions in \code{querySVs} with the exact, simulated breakpoints \code{simSVs}. Positions in \code{querySVs} may only differ by at most \code{tol} base pairs from the simulated breakpoints.}
}

\details{
  An overlap is defined as the overlap between the breakpoints / breakpoint regions up to the given tolerance in bp. Overlap does not mean the whole affected region between the start and end of the SV. \cr
  The comparison has to be done for each type of SV separately.
  It is required to use the returned tables from \code{\link{simulateSV}} for the argument \code{simSVs}. For \code{querySVs}, the tables have to be in BED- or BEDPE-format, a simple tab-separated table with genome coordinates and an optional column with the breakpoint sequence (bpSeq):
  \itemize{
    \item Deletions: (1) BED with columns chr, start, end (, bpSeq) for exact breakpoints or (2) BEDPE with columns chr, start1, end1, chr, start2, end2 (, bpSeq) for approximate breakpoint regions
    \item Insertions: BEDPE with columns chrA, startA, endA, chrB, startB, endB (, bpSeq). Typically, a complete insertion is reported in two rows, one for the breakpoint at the 5' and one at the 3' end
    \item Inversions: (1) BED with columns chr, start, end (, bpSeq1, bpSeq2) or (2) BEDPE with columns chr, start1, end1, chr, start2, end2 (, bpSeq1, bpSeq2). Inversions, the larger ones, typically have two breakpoint sequences (one for the 5' and on for the 3' end)
    \item Tandem duplications: (1) BED with columns chr, start, end (, bpSeq) or (2) BEDPE with columns chr, start1, end1, chr, start2, end2 (, bpSeq)
    \item Translocations: BEDPE with columns chrA, startA, endA, chrB, startB, endB (, bpSeq1, bpSeq2)
  }
  The BEDPE-format is required for insertions and translocations, since they involve two regions on the genome. For other SVs, the BEDPE allows to specify approximate regions for each breakpoint. For example a deletion "chr:start-end" can be given as BED-file with columns chr, start, end or (a little redundant) as BEDPE-file with columns chr, start, start, chr, end, end or (with some tolerance) like chr, start-tol, start+tol, chr, end-tol, end+tol. The tolerance can also be regulated by the function argument \code{tol}.
}

\value{
The table of simulated SVs, as given in the function argument \code{querySVs}, but with additional columns for the overlapping region in \code{querySVs} and the percentage overlap between the breakpoint sequences (if they were provided as a column in \code{querySVs}). \cr
Furthermore, the function prints the sensitivity and precision in the R console.
}

\references{More information about the BED-format can be found on the BEDTools homepage: \url{http://code.google.com/p/bedtools}}

\author{Christoph Bartenhagen}

\seealso{\code{\link{simulateSV}}}

\examples{

  ## Toy example: Artificial genome with two chromosomes
  genome = DNAStringSet(c("AAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTT", "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCC"))
  names(genome) = c("chr1","chr2")

  #######################
  ## Example 1: Deletions
  ## Simulation of 5 deletions of 5bp each
  sim = simulateSV(output=NA, genome=genome, dels=5, sizeDels=5, bpSeqSize=10, seed=246)
  simSVs = metadata(sim)$deletions

  ## An SV detection in BED format may look like this:
  ## Four of five deletions were detected; two with exact and two with an approximate breakpoint
  ## Two additional deletions were detected, which were not part of the simulation
  ## The column with the breakpoint sequence is optional, the column names not important (BED-files have no header)
  querySVs = data.frame(
    chr=c("chr1","chr1","chr1","chr2","chr2","chr2"), 
    start=c(4,12,20,10,21,34), 
    end=c(8,16,28,14,31,38), 
    bpSeq=c("AAAAAAAAAA", "AAAAAAAAAT", "AAAATTTTTT", "GGGGGGGGGG", "GGGGGGGCCC", "CCCCCCCCCC")
  )
  
  ## Compare the SVs with 0bp tolerance:
  ## Only the two exact detections have an overlap
  simSVs_overlap1 = compareSV(querySVs, simSVs, tol=0)
  simSVs_overlap1

  ## Increasing the breakpoint tolerance to +/- 3bp :
  ## Now, the overlap also includes the more imprecise detections
  ## And the sensitivity and precision are better
  ## Note that for deletion2, the breakpoint sequence matches only by 50%
  simSVs_overlap2 = compareSV(querySVs, simSVs, tol=3)
  simSVs_overlap2

  ############################
  ## Example 2: Translocations
  ## Simulation of 2 translocations (only one of them is balanced):
  sim = simulateSV(output=NA, genome=genome, trans=2, percBalancedTrans=0.5, bpSeqSize=10, seed=246)
  simSVs = metadata(sim)$translocations
  
  ## Detected translocations have to be given in BEDPE-format (i.e. at least six columns with chr,start,end for breakpoints on both chromosomes)
  ## In this example, the breakpoints were approximated up to 1 or 2bp
  ## Optional breakpoint sequences are missing
  querySVs = data.frame(
    chr=c("chr2", "chr1", "chr2"), 
    start1=c(25,3,9), 
    end1=c(29,7,12), 
    chr2=c("chr1","chr2","chr1"),
    start2=c(22,10,3),
    end2=c(25,13,4)
  )

  simSVs_overlap = compareSV(querySVs, simSVs, tol=0)
  simSVs_overlap

}
