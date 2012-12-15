\name{simulateSV}

\alias{simulateSV}
\alias{simulateSV,ANY}
\alias{simulateSV,ANY-method}

\title{Structural Variant Simulation}

\description{A tool for simulating deletions, insertions, inversions, tandem duplications and translocations in any genome available as FASTA-file or BSgenome data package. Structural variations (SVs) are placed within the given genome, or only a subset of it, in a random, non-overlapping manner or at given genomic coordinates.}

\usage{simulateSV(output=".", genome, chrs, dels=0, ins=0, invs=0, dups=0, trans=0, size, sizeDels=10, sizeIns=10, sizeInvs=10, sizeDups=10, regionsDels, regionsIns, regionsInvs, regionsDups, regionsTrans, maxDups=10, percCopiedIns=0, percBalancedTrans=1, transInsert=0, bpSeqSize=100, random=TRUE, seed)}

\arguments{
  \item{output}{Output directory for the rearranged genome and SV lists; turn this off by passing \code{NA} (default: current directory)}
  \item{genome}{The genome as \code{DNAStringSet} or as filename pointing to a FASTA-file containing the genome sequence}
  \item{chrs}{Restrict simulation to certain chromosomes only (default: all chromosomes available)}
  \item{dels}{Number of deletions}
  \item{ins}{Number of insertions}
  \item{invs}{Number of inversions}
  \item{dups}{Number of tandem duplications}
  \item{trans}{Number of translocations}
  \item{size}{Size of SVs in bp (a single numeric value); a quick way to set a size, which is applied to all simulated SVs}
  \item{sizeDels}{Size of deletions: Either a single number for all deletions or a vector with a length for every single deletion}
  \item{sizeIns}{Size of insertions: Either a single number for all insertions or a vector with a length for every single insertion}
  \item{sizeInvs}{Size of inversions: Either a single number for all inversions or a vector with a length for every single inversion}
  \item{sizeDups}{Size of tandem duplications: Either a single number for all tandem duplications or a vector with a length for every single tandem duplication}
  \item{regionsDels}{\code{GRanges} object with regions within the genome where to place the deletions}
  \item{regionsIns}{\code{GRanges} object with regions within the genome where to place the Insertions}
  \item{regionsInvs}{\code{GRanges} object with regions within the genome where to place the inversions}
  \item{regionsDups}{\code{GRanges} object with regions within the genome where to place the tandem duplications}
  \item{regionsTrans}{\code{GRanges} object with regions within the genome where to place the translocations}
  \item{maxDups}{Maximum number of repepeats for tandem duplications}
  \item{percCopiedIns}{Percentage of copy-and-paste-like insertions (default: 0, i.e. no inserted sequences are duplicated)}
  \item{percBalancedTrans}{Percentage of balanced translocations (default: 1, i.e. all translocations are balanced)}
  \item{transInsert}{Maximum number of randomly inserted nucleotides at each translocation breakpoint in bp (default: 0bp, i.e. off)}
  \item{bpSeqSize}{Length of the breakpoint sequences in the output}
  \item{random}{If TRUE, the SVs will be placed randomly within the genome or the given regions; otherwise, the given regions will be used as SV coordinates (\code{random} can also be a vector of five elements with TRUE/FALSE for every SV in the following order: deletions, insertions, inversions, duplications, translocations)}
  \item{seed}{Fixed seed for generation of random SV positions}
}

\details{
  About the SV types supported:
  \itemize{
    \item Deletions: A segment is cut out from the genome.
    \item Insertions: A segment is cut or copied (see parameter \code{percCopiedIns}) from one chromosome A and inserted into another chromosome B.
    \item Inversions: A segment is cut out from one chromosome and its reverse complement is inserted at the same place without loss or a shift of sequence.
    \item Tandem duplication: A segment is duplicated one after the other at most \code{maxDups} times.
    \item Translocation: A segment from the 5' or 3' end of one chromosome A is exchanged with the 5' or 3' end of another chromosome B. If it is not balanced (see parameter \code{percBalancedTrans}), the segment from chromosome B will be lost, what results in a duplicated sequence from chromosome A. Segments translocated between two different ends (5'<->3' or 3'<->5') are always inverted.
  }
  About the function arguments:
  \itemize{
    \item The region arguments (\code{regionsDels},\code{regionsIns}, \code{regionsInvs},\code{regionsDups},\code{regionsTrans}) can be used in two ways: 1. as subsets of the genome where to place the SVs randomly or 2. as coordinates of SVs that shall be implemented at these exact positions. The latter can be useful to implement a predefined set of previously detected or known SVs. Set the parameter \code{random} to FALSE accordingly. In this case, the rownames of the region arguments will be used to name the SVs in the output.
    \item In case of insertions and translocations, where two genomic regions are involved, add columns "chrB" and "startB" for insertions and "chrB", "startB", "endB" for translocations to the \code{GRanges} objects in \code{regionsTrans} and \code{regionsIns}, respectively.
    \item It is recommended to set the size of an SV individually for every deletion, insertion, inversion or tandem duplication. See the function \command{estimateSVSizes} to estimate beta-distributed sizes from a training set of SVs (defaults for deletions, insertions, inversions and tandem duplications are available). Using the argument \code{size} to set the size for all SVs overrides all other size arguments.
    \item There is no size argument for translocations. After random generation of the breakpoint, the translocation spans the chromosome until the closest of both chromosome ends.
  }
  Misc:
  \itemize{
    \item By default, the human genome (hg19) will be used, which requires the package \pkg{BSgenome.Hsapiens.UCSC.hg19}.
    \item SVs will not be placed within unannotated regions or assembly gaps denoted by the letter "N". These regions are detected and filtered automatically.
  }
}

\value{
  The rerranged genome as a \code{DNAStringSet}. Its \code{metadata} slot contains a named list of \code{data.frames} with information about the simulated SVs:
  \item{deletions}{The coordinates of the implemented deletions and the breakpoint sequence}
  \item{insertions}{The coordinates of the origin (chrA) and destination (chrB) of the inserted sequence and the breakpoints sequences at both ends (5' and 3'). If the sequence is cut out from the original chromosome A, the sequence of this breakpoint is given as well.}
  \item{inversions}{The coordinates of the implemented inversions and the breakpoint sequences at both ends (5' and 3')}
  \item{tandemDuplications}{The coordinates of the duplicated sequence and the breakpoint sequence}
  \item{translocations}{The coordinates of the translocated sequences and the two breakpoint sequences (if balanced)}
  The coordinates in the tables refer to the "normal" reference genome. \cr
  All the list items can also be written to the specified output directory (which is the current directory by default). The genome will be saved in FASTA format and the SVs \code{data.frames} as CSV tables.
}

\author{Christoph Bartenhagen}

\seealso{\code{\link{estimateSVSizes}}}

\examples{

## Toy example: Artificial genome with two chromosomes
genome = DNAStringSet(c("AAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTT", "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCC"))
names(genome) = c("chr1","chr2")

## Three deletions of sizes 10bp each
sim = simulateSV(output=NA, genome=genome, dels=3, sizeDels=10, bpSeqSize=10, seed=246)
sim
metadata(sim)

## Three insertions of 5bp each; all cut-and-paste-like (default)
sim = simulateSV(output=NA, genome=genome, ins=3, sizeIns=5, bpSeqSize=10, seed=246)
sim
metadata(sim)
## Three insertions of 5bp each; all copy-and-paste-like (note the parameter \code{percCopiedIns})
sim = simulateSV(output=NA, genome=genome, ins=3, sizeIns=5, percCopiedIns=1, bpSeqSize=10, seed=246)
sim
metadata(sim)

## Three inversions of sizes 2bp, 4bp and 6bp
sim = simulateSV(output=NA, genome=genome, invs=3, sizeInvs=c(2,4,6), bpSeqSize=10, seed=246)
sim
metadata(sim)

## A tandem duplication of 4bp with at most ten duplications
## The duplication shall be placed somewhere within chr4:18-40
library(GenomicRanges)
region = GRanges(IRanges(10,30),seqnames="chr1")
sim = simulateSV(output=NA, genome=genome, dups=1, sizeDups=4, regionsDups=region, maxDups=10, bpSeqSize=10, seed=246)
sim
metadata(sim)

## A balanced translocation (default)
sim = simulateSV(output=NA, genome=genome,trans=1, bpSeqSize=6, seed=246)
sim
metadata(sim)
## The same translocation, but unbalanced (note the parameter \code{percBalancedTrans})
sim = simulateSV(output=NA, genome=genome, trans=1, percBalancedTrans=0, bpSeqSize=6, seed=246)
sim
metadata(sim)

## Simulate all four SV types at once:
## 2 deletions (5bp), 2 insertions (5bp),2 inversions (3bp), 1 tandem duplication (4bp), 1 translocations
sim = simulateSV(output=NA, genome=genome, dels=2, ins=2, invs=2, dups=1, trans=1, sizeDels=5, sizeIns=5, sizeInvs=3, sizeDups=4, maxDups=3, percCopiedIns=0.5, bpSeqSize=10, seed=246)
sim
metadata(sim)

## Avoid random generation of coordinates and implement a given deletion of 10bp on chr2:16-25
knownDeletion = GRanges(IRanges(16,25), seqnames="chr2")
names(knownDeletion) = "myDeletion"
knownDeletion
sim = simulateSV(output=NA, genome=genome, regionsDels=knownDeletion, bpSeqSize=10, random=FALSE)
sim
metadata(sim)

## Avoid random generation of coordinates and implement a given insertion from chr1:16:25 at chr2:26
knownInsertion = GRanges(IRanges(16,25), seqnames="chr1", chrB="chr2", startB=26)
names(knownInsertion) = "myInsertion"
knownInsertion
sim = simulateSV(output=NA, genome=genome, regionsIns=knownInsertion, bpSeqSize=10, random=FALSE)
sim
metadata(sim)

## This example simulates a translocation t(9;22) leading to the BCR-ABL fusion gene.
## It uses simple breakpoints within 9q34.1 and 22q11.2 for demonstration
## Take care to add coordinates of both chromosomes to the GRanges object:
trans_BCR_ABL = GRanges(IRanges(133000000,141213431), seqnames="chr9", chrB="chr22", startB=23000000, endB=51304566, reciprocal=TRUE)
names(trans_BCR_ABL) = "BCR_ABL"
trans_BCR_ABL
## This example requires the \pkg{BSgenome.Hsapiens.UCSC.hg19} which is used by default (hence, no genome argument)
## Setting the argument \code{transInsert=20} adds up to 20 random nucleotides at both breakpoints
\dontrun{sim = simulateSV(output=NA, chrs=c("chr9", "chr22"), regionsTrans=trans_BCR_ABL, transInsert=20, bpSeqSize=30, random=FALSE)}
}

\keyword{datagen}

