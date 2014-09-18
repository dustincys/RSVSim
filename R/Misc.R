setMethod("estimateSVSizes",
          signature(n="numeric", svSizes="missing", default="character"),
          function(n, svSizes, minSize, maxSize, default, hist){
            if(!(default %in% c("deletions", "insertions", "inversions", "tandemDuplications"))){
              stop("Invalid argument: No default values for SV type ", default)
            }
            ## Default shape parameters were estimated from DGV release 2012-03-29
            ## Only studies with "method=sequencing" and size between 500bp and 10kb were selected (for insertions the extreme peak between 6000bp and 6200bp has been removed)
            ## SVs left: 1129 Deletions, 490 Insertions, 202 Inversions, 145 Tandem Duplications
            if(default == "deletions"){
              shape1=0.3692374
              shape2=2.689951
            }
            if(default == "insertions"){
              shape1=0.7106429
              shape2=2.274061
            }
            if(default == "inversions"){
              shape1=0.5786958
              shape2=1.83066
            }
            if(default == "tandemDuplications"){
              shape1=0.4084596
              shape2=1.998082
            }
            if(any(is.na(c(minSize, maxSize)))){
              return(.getSVSizes(n=n, shape1=shape1, shape2=shape2, min=500, max=10000, hist=hist))
            }else{
              return(.getSVSizes(n=n, shape1=shape1, shape2=shape2, min=minSize, max=maxSize, hist=hist))
            }
          })

setMethod("estimateSVSizes",
          signature(n="numeric", svSizes="missing", default="missing"),
          function(n, svSizes, minSize, maxSize, default, hist){
            stop("Missing argument(s): Please give either a set of SVs or choose a set of default values.")
          })

setMethod("estimateSVSizes",
          signature(n="numeric", svSizes="numeric", default="missing"),
          function(n, svSizes, minSize, maxSize, default, hist){

            require(MASS)
            
            ## calculate the two shape parameters from given set of SVs
            min = min(svSizes)
            max = max(svSizes)            

            svSizes = (svSizes-min)/(max-min+1)

            ## make sure values are >0 and <1
            svSizes[svSizes == 0] = svSizes[svSizes == 0] + 1e-5
            svSizes[svSizes == 1] = svSizes[svSizes == 1] - 1e-5

#            m = mean(svSizes)
#            m = (m-min)/(max-min)
#            v = var(svSizes)
#            v = v/(max-min)^2
            
            ## shape parameters for beta distribution
#            shape1 = m*(((m*(1-m))/v)-1)
#            shape2 = (1-m)*(((m*(1-m))/v)-1)
            
            beta = suppressWarnings(fitdistr(svSizes, "beta", start=list(shape1=0.1, shape2=0.1)))
            shape1 = beta$estimate[1]
            shape2 = beta$estimate[2]
            
            if(any(c(is.na(minSize), is.na(maxSize)))){
              return(.getSVSizes(n=n, shape1=shape1, shape2=shape2, min=min, max=max, hist=hist))
            }else{
              return(.getSVSizes(n=n, shape1=shape1, shape2=shape2, min=minSize, max=maxSize, hist=hist))
            }
          })

.getSVSizes <- function(n, shape1, shape2, min, max, hist){

  simSizes = rbeta(n, shape1=shape1, shape2=shape2)
  simSizes = round((simSizes*(max-min+1))+min)
  if(hist == TRUE){
    hist(simSizes)
  }
  return(simSizes)
  
}

setMethod("compareSV",
          signature("data.frame", "data.frame"),
          function(querySVs, simSVs, tol){
            simSVs_overlap = .compareSV(querySVs, simSVs, tol)
            return(simSVs_overlap)
          })

setMethod("compareSV",
          signature("character", "data.frame"),
          function(querySVs, simSVs, tol){
            querySVs= read.table(querySVs, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
            simSVs_overlap = .compareSV(querySVs, simSVs, tol)
            return(simSVs_overlap)
          })

setMethod("compareSV",
          signature("character", "character"),
          function(querySVs, simSVs, tol){
            simSVs = read.table(simSVs, fill=TRUE, header=TRUE)
            querySVs= read.table(querySVs, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
            simSVs_overlap = .compareSV(querySVs, simSVs, tol)
            return(simSVs_overlap)
          })

## deletions: (1) BED with colums chr,start,end(,bpSeq) or (2) BEDPE with colums chr,start1,end1,chr,start2,end2(,bpSeq)
## insertions: BEDPE with columns chrA,startA,endA, chrB,startB,endB(, bpSeq); two rows for one SV (5' and 3' breakpoint respectively)
## inversions: (1) BED with colums chr,start,end(,bpSeq1, bpSeq2) or (2) BEDPE with colums chr,start1,end1,chr,start2,end2(,bpSeq1, bpSeq2)
## tandem duplications: (1) BED with colums chr,start,end or (2) BEDPE with colums chr,start1,end1,chr,start2,end2
## translocations: BEDPE with columns chrA,startA,endA, chrB,startB,endB, (bpSeq1, bpSeq2)
.compareSV <- function(querySVs, simSVs, tol){
  
  ## get the type of SV by checking the column names in simSVs
  type = NA
  if(all(c("Name", "Chr", "Start", "End", "Size",  "BpSeq") %in% colnames(simSVs))){
    type = "deletion"
  }
  if(all(c("Name", "Chr", "Start", "End", "Size", "BpSeq_5prime", "BpSeq_3prime") %in% colnames(simSVs))){
    type = "inversion"
  }
  if(all(c("Name", "ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Size", "Copied", "BpSeqA", "BpSeqB_5prime", "BpSeqB_3prime") %in% colnames(simSVs))){
    type = "insertion"
  }
  if(all(c("Name", "Chr", "Start", "End", "Size", "Duplications", "BpSeq") %in% colnames(simSVs))){
    type = "tandemDuplication"
  }
  if(all(c("Name", "ChrA", "StartA", "EndA", "SizeA", "ChrB", "StartB", "EndB", "SizeB", "Balanced", "BpSeqA", "BpSeqB") %in% colnames(simSVs))){
    type = "translocation"
  }
  if(is.na(type)){
    stop("Invalid input: Format of the second set of SVs differs from the output format of the simulator.")
  }
  
  ## load user-supplied set of svs (BED- or BEDPE-file)
  querySVs_Bp1 = querySVs_Bp1 = NULL
  ## split BED input into appropriate GRanges objects, which depends on the number of given columns (BED or BEDPE, with or without bpSeq)
  ## give every SV an id to remember which breakpoints belong together after GRanges sorted the coordinates
  ## inversions
  if(ncol(querySVs) == 8){
    querySVs_Bp1 = GRanges(IRanges(querySVs[, 2], querySVs[, 3]), seqnames=querySVs[, 1], bpSeq1=querySVs[, 7], bpSeq2=querySVs[, 8], id=1:nrow(querySVs))
    querySVs_Bp2 = GRanges(IRanges(querySVs[, 5], querySVs[, 6]), seqnames=querySVs[, 4], bpSeq1=querySVs[, 7], bpSeq2=querySVs[, 8], id=1:nrow(querySVs))
  }
  ## translocations, insertions or deletions
  if(ncol(querySVs) == 7){
    querySVs_Bp1 = GRanges(IRanges(querySVs[, 2], querySVs[, 3]), seqnames=querySVs[, 1], bpSeq=querySVs[, 7], id=1:nrow(querySVs))
    querySVs_Bp2 = GRanges(IRanges(querySVs[, 5], querySVs[, 6]), seqnames=querySVs[, 4], bpSeq=querySVs[, 7], id=1:nrow(querySVs))
  }
  ## translocations, inversions, insertions or deletions (no bpSeq)
  if(ncol(querySVs) == 6){
    querySVs_Bp1 = GRanges(IRanges(querySVs[, 2], querySVs[, 3]), seqnames=querySVs[, 1], bpSeq=NA, bpSeq2=NA, id=1:nrow(querySVs))
    querySVs_Bp2 = GRanges(IRanges(querySVs[, 5], querySVs[, 6]), seqnames=querySVs[, 4], bpSeq=NA, bpSeq2=NA, id=1:nrow(querySVs))
  }
  ## inversions
  if(ncol(querySVs) == 5){
    querySVs_Bp1 = GRanges(IRanges(querySVs[, 2], querySVs[, 2]), seqnames=querySVs[, 1], bpSeq1=querySVs[, 4], bpSeq2=querySVs[, 5], id=1:nrow(querySVs))
    querySVs_Bp2 = GRanges(IRanges(querySVs[, 3], querySVs[, 3]), seqnames=querySVs[, 1], bpSeq1=querySVs[, 4], bpSeq2=querySVs[, 5], id=1:nrow(querySVs))
  }
  ## deletions
  if(ncol(querySVs) == 4){
    querySVs_Bp1 = GRanges(IRanges(querySVs[, 2], querySVs[, 2]), seqnames=querySVs[, 1], bpSeq=querySVs[, 4], id=1:nrow(querySVs))
    querySVs_Bp2 = GRanges(IRanges(querySVs[, 3], querySVs[, 3]), seqnames=querySVs[, 1], bpSeq=querySVs[, 4], id=1:nrow(querySVs))
  }
  ## deletions, inversions (no bpSeq)
  if(ncol(querySVs) == 3){
    querySVs_Bp1 = GRanges(IRanges(querySVs[, 2], querySVs[, 2]), seqnames = querySVs[, 1], bpSeq=NA, bpSeq1=NA, bpSeq2=NA, id=1:nrow(querySVs))
    querySVs_Bp2 = GRanges(IRanges(querySVs[, 3], querySVs[, 3]), seqnames = querySVs[, 1], bpSeq=NA, bpSeq1=NA, bpSeq2=NA, id=1:nrow(querySVs))
  }
  if(all(is.null(querySVs_Bp1), is.null(querySVs_Bp2))){
    stop("Invalid SV file: Please make sure, that the SV file is in BED (chr,start,end) or BEDPE format (chr1,start1,end1,chr2,start2,end2) and contains only the breakpoint sequence(s) as additional column")
  }
  numOverlaps = 0
  
  ## For deletions and tandem duplications compare the 1 breakpoint and the 1 breakpoint sequence
  if(type == "deletion" | type == "tandemDuplication"){
    simSVs$Overlap = ""
    simSVs$OverlapBpSeq = NA
    for(v in 1:nrow(simSVs)){
      simSV = simSVs[v, ]
      ## check correct start and end coordinate
      chr = as.character(simSV$Chr)
      simSV_Bp1 = GRanges(IRanges(simSV$Start, simSV$Start), seqnames=chr)
      simSV_Bp2 = GRanges(IRanges(simSV$End, simSV$End), seqnames=chr)
      overlap = .getOverlap(list(simSV_Bp1, simSV_Bp2), list(querySVs_Bp1[seqnames(querySVs_Bp1) == chr], querySVs_Bp2[seqnames(querySVs_Bp2) == chr]), list(simSV$BpSeq, querySVs_Bp1[seqnames(querySVs_Bp1) == chr]$bpSeq), list(NA, NA), tol, type)
      if(!is.null(overlap)){
        simSVs$Overlap[v] = overlap[[1]]
        simSVs$OverlapBpSeq[v] = overlap[[2]]
        numOverlaps = numOverlaps + overlap[[4]]
      }
    }
  }
  ## For inversions compare the 1 breakpoint and the 2 breakpoint sequences
  ## seq|invseq|seq
  if(type == "inversion"){
    simSVs$Overlap = ""
    simSVs$OverlapBpSeq_5prime = NA
    simSVs$OverlapBpSeq_3prime = NA
    for(v in 1:nrow(simSVs)){
      simSV = simSVs[v, ]
      ## check correct start and end coordinate
      chr = as.character(simSV$Chr)
      simSV_Bp1 = GRanges(IRanges(simSV$Start, simSV$Start), seqnames=chr)
      simSV_Bp2 = GRanges(IRanges(simSV$End, simSV$End), seqnames=chr)
      overlap = .getOverlap(list(simSV_Bp1, simSV_Bp2), list(querySVs_Bp1[seqnames(querySVs_Bp1) == chr], querySVs_Bp2[seqnames(querySVs_Bp2) == chr]), list(simSV$BpSeq_5prime, querySVs_Bp1[seqnames(querySVs_Bp1) == chr]$bpSeq1), list(simSV$BpSeq_3prime, querySVs_Bp1[seqnames(querySVs_Bp1) == chr]$bpSeq2), tol, type)
      if(!is.null(overlap)){
        simSVs$Overlap[v] = overlap[[1]]
        simSVs$OverlapBpSeq_5prime[v] = overlap[[2]]
        simSVs$OverlapBpSeq_3prime[v] = overlap[[3]]
        numOverlaps = numOverlaps + overlap[[4]]
      }
    }
  }
  ## For insertions compare the 2 breakpoints and 2 breakpoint sequences
  ## chrB|chrA|chrB
  if(type == "insertion"){
    simSVs$Overlap_5prime = ""
    simSVs$Overlap_3prime = ""
    simSVs$OverlapBpSeq_5prime = NA
    simSVs$OverlapBpSeq_3prime = NA
    for(v in 1:nrow(simSVs)){
      simSV = simSVs[v, ]
      ## Compare breakpoint at chrB|chrA (5', i.e. start)
      simSV_Bp1 = GRanges(IRanges(simSV$StartB, simSV$StartB), seqnames=simSV$ChrB)
      simSV_Bp2 = GRanges(IRanges(simSV$StartA, simSV$StartA), seqnames=simSV$ChrA)
      chr1 = as.character(seqnames(simSV_Bp1))
      chr2 = as.character(seqnames(simSV_Bp2))
      querySVs_Bp12 = c(querySVs_Bp1, querySVs_Bp2)
      overlap = .getOverlap(list(simSV_Bp1, simSV_Bp2), list(querySVs_Bp12[seqnames(querySVs_Bp12) == chr1], querySVs_Bp12[seqnames(querySVs_Bp12) == chr2]), list(simSV$BpSeqB_5prime, querySVs_Bp12[seqnames(querySVs_Bp12) == chr1]$bpSeq), list(NA, NA), tol, type)
      if(!is.null(overlap)){
        simSVs$Overlap_5prime[v] = overlap[[1]]
        simSVs$OverlapBpSeq_5prime[v] = overlap[[2]]
        numOverlaps = numOverlaps + overlap[[4]]
      }
      ## Compare second breakpoint at chrA|chrB (3', i.e. end)
      simSV_Bp1 = GRanges(IRanges(simSV$EndA, simSV$EndA), seqnames=simSV$ChrA)
      simSV_Bp2 = GRanges(IRanges(simSV$StartB, simSV$StartB), seqnames=simSV$ChrB)
      chr1 = as.character(seqnames(simSV_Bp1))
      chr2 = as.character(seqnames(simSV_Bp2))
      overlap = .getOverlap(list(simSV_Bp1, simSV_Bp2), list(querySVs_Bp12[seqnames(querySVs_Bp12) == chr1], querySVs_Bp12[seqnames(querySVs_Bp12) == chr2]), list(simSV$BpSeqB_3prime, querySVs_Bp12[seqnames(querySVs_Bp12) == chr1]$bpSeq), list(NA, NA), tol, type)
      if(!is.null(overlap)){
        simSVs$Overlap_3prime[v] = overlap[[1]]
        simSVs$OverlapBpSeq_3prime[v] = overlap[[2]]
        numOverlaps = numOverlaps + overlap[[4]]
      }
      
    }
  }
  ## For translocations compare the 2 breakpoints and 1 breakpoint sequence
  if(type == "translocation"){
    simSVs$Overlap = ""
    simSVs$OverlapBpSeq_A = NA
    simSVs$OverlapBpSeq_B = NA
    for(v in 1:nrow(simSVs)){
      simSV = simSVs[v, ]
      ## Breakpoints can be either on 5' terminal or on 3' terminal
      ## Only check one breakpoint since the coordinates are only switched between chrB|chrA and chrA|chrB, but practically the same
      if(simSV$StartB == 1){
        simSV_Bp1 = GRanges(IRanges(simSV$EndB, simSV$EndB), seqnames=simSV$ChrB) ## 5', i.e. check the ends
      }else{
        simSV_Bp1 = GRanges(IRanges(simSV$StartB, simSV$StartB), seqnames=simSV$ChrB) ## 3', i.e. check the starts
      }
      if(simSV$StartA == 1){
        simSV_Bp2 = GRanges(IRanges(simSV$EndA, simSV$EndA), seqnames=simSV$ChrA) ## 5', i.e. check the ends
      }else{
        simSV_Bp2 = GRanges(IRanges(simSV$StartA, simSV$StartA), seqnames=simSV$ChrA) ## 3', i.e. check the starts
      }
      chr1 = as.character(seqnames(simSV_Bp1))
      chr2 = as.character(seqnames(simSV_Bp2))
      querySVs_Bp12 = c(querySVs_Bp1, querySVs_Bp2)
      ## If translocation is not balanced, the breakpoint sequence on chrA does not need to be checked
      if(simSV$Balanced == FALSE){
        overlap = .getOverlap(list(simSV_Bp1, simSV_Bp2), list(querySVs_Bp12[seqnames(querySVs_Bp12) == chr1], querySVs_Bp12[seqnames(querySVs_Bp12) == chr2]), list(simSV$BpSeqB, querySVs_Bp12[seqnames(querySVs_Bp12) == chr1]$bpSeq), list(NA, NA), tol, type)
        if(!is.null(overlap)){
          simSVs$Overlap[v] = overlap[[1]]
          simSVs$OverlapBpSeq_B[v] = overlap[[2]]
          numOverlaps = numOverlaps + overlap[[4]]
        }
      }else{
        overlap = .getOverlap(list(simSV_Bp1, simSV_Bp2), list(querySVs_Bp12[seqnames(querySVs_Bp12) == chr1], querySVs_Bp12[seqnames(querySVs_Bp12) == chr2]), list(simSV$BpSeqB, querySVs_Bp12[seqnames(querySVs_Bp12) == chr1]$bpSeq), list(simSV$BpSeqA, querySVs_Bp12[seqnames(querySVs_Bp12) == chr2]$bpSeq), tol, type)
        if(!is.null(overlap)){
          simSVs$Overlap[v] = overlap[[1]]
          simSVs$OverlapBpSeq_B[v] = overlap[[2]]
          simSVs$OverlapBpSeq_A[v] = overlap[[3]]
          numOverlaps = numOverlaps + overlap[[4]]
        }
      }
    }
  }

  ## report sensitivity and precision
  if(type == "insertion"){
    tp = numOverlaps
    fp = nrow(querySVs) - tp
    fn = sum(simSVs$Overlap_5prime == "") + sum(simSVs$Overlap_3prime == "")
  }
  else{
    tp = numOverlaps
    fp = nrow(querySVs) - tp
    fn = sum(simSVs$Overlap == "")
  }
  sensitivity = round((tp / (tp + fn)) * 100)
  precision = round((tp / (tp + fp)) * 100)
  message("Sensitivity: ", sensitivity, "%")
  message("Precision: ", precision, "%")

  return(simSVs)
  
}

                         
.getOverlap <- function(sim, query, bpSeqA, bpSeqB, tol, type){

  ## simple substitution matrix, which only counts matches
  substitutionMatrix = nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE, type = "DNA")

  start(query[[1]]) = start(query[[1]]) - tol
  end(query[[1]]) = end(query[[1]]) + tol
  start(query[[2]]) = start(query[[2]]) - tol
  end(query[[2]]) = end(query[[2]]) + tol
  overlap1 = as.matrix(findOverlaps(sim[[1]], query[[1]]))
  overlap2 = as.matrix(findOverlaps(sim[[2]], query[[2]]))
  if(nrow(overlap1) > 0 & nrow(overlap2) > 0){
    
    hit1 = overlap1[, "subjectHits"]
    hit2 = overlap2[, "subjectHits"]

#    if(query[[1]]$id[hit1] == query[[2]]$id[hit2]){
    ## for insertions and translocations, which have two breakpoints on two chromosomes, return two overlapping regions (separated by "; ")
    ## deletions, inversions and tandem duplications have only start and end, so return just one overlapping region
    ## returned regions include the given tolerance
    if(type == "insertion" | type == "translocation"){
      overlap = paste(paste(seqnames(query[[1]])[hit1], ":", start(query[[1]])[hit1]+tol, "-", end(query[[1]])[hit1]-tol, ", ", seqnames(query[[2]])[hit2], ":", start(query[[2]])[hit2]+tol, "-", end(query[[2]])[hit2]-tol, sep=""), collapse=", ")
    }else{
      overlap = paste(paste(seqnames(query[[1]])[hit1], ":", start(query[[1]])[hit1]+tol, "-", end(query[[2]])[hit2]-tol, sep=""), collapse=", ")
    }
    numOverlaps = length(hit1)
    bpSeqAlnScoreA = bpSeqAlnScoreB = NA
    if(!any(is.na(bpSeqA[[2]]))){
      ## align all overlaps and report the one with the maximum overlap
      bpSeqAlnScoreA = vector()
      for(h in hit1){
        bpSeqAlnScoreA = c(bpSeqAlnScoreA, pairwiseAlignment(bpSeqA[[1]], bpSeqA[[2]][h], type="local", substitutionMatrix=substitutionMatrix, gapOpening=0, gapExtension=0, scoreOnly=TRUE))
      }
      bpSeqAlnScoreA = max(bpSeqAlnScoreA)
      bpSeqAlnScoreA = round((bpSeqAlnScoreA / nchar(bpSeqA[[1]])) * 100)
    }
    if(!any(is.na(bpSeqB[[2]]))){
      bpSeqAlnScoreB = vector()
      for(h in hit2){
        bpSeqAlnScoreB = c(bpSeqAlnScoreB, pairwiseAlignment(bpSeqB[[1]], bpSeqB[[2]][h], type="local", substitutionMatrix=substitutionMatrix, gapOpening=0, gapExtension=0, scoreOnly=TRUE))
      }
      bpSeqAlnScoreB = max(bpSeqAlnScoreB)
      bpSeqAlnScoreB = round((bpSeqAlnScoreB / nchar(bpSeqB[[1]])) * 100)
    }
    return (list(overlap, bpSeqAlnScoreA, bpSeqAlnScoreB, numOverlaps))
  }
  return (NULL)

}


.readRepeatMaskerOutput <- function(file, save=TRUE){
  t = read.table(file, skip=2, fill=TRUE, stringsAsFactors=FALSE, colClasses=c("numeric","numeric","numeric","numeric","character","numeric","numeric","character","character","character","character","character","character","character","numeric"))

  t = t[t[, 5] %in% paste("chr", c(1:22,"X","Y"), sep=""), ]
  t = t[t[, 11] %in% c("LINE/L1","LINE/L2","SINE/Alu","SINE/MIR"), ]
  t[, 11] = gsub("LINE/","", t[, 11])
  t[, 11] = gsub("SINE/","", t[, 11])
  repeats = GRanges(IRanges(t[, 6], t[, 7]), seqnames=t[, 5], type=t[, 11])

  linesL1 = repeats[repeats$type == "L1"]
  linesL2 = repeats[repeats$type == "L2"]
  sinesAlu = repeats[repeats$type == "Alu"]
  sinesMIR = repeats[repeats$type == "MIR"]
  ## load previously saved segmental duplications to avoid import from UCSC browser
  if(file.exists(file.path(path.package("RSVSim"), "data", "segmentalDuplications.RData"))){
    data("segmentalDuplications", package="RSVSim", envir=environment())
  }else{
    require(rtracklayer)
    segDups = .loadFromUCSC_SegDups()
  }
  trs = .loadFromBSGenome_TandemRepeats()
  repeats = list(linesL1, linesL2, sinesAlu, sinesMIR, segDups, trs)
  if(save == TRUE){
    save(repeats, file=file.path(path.package("RSVSim"), "data", "repeats_hg19.RData"), compress="xz")
   }
  return(repeats)
}

## Loading the rmsk-track from UCSC may take up to 40 Minutes
.loadFromUCSC_RepeatMasks <- function(save=TRUE, verbose){
  
  ## from UCSV
  require(rtracklayer)
  chrs = paste("chr", c(1), sep="")
  mySession = browserSession("UCSC")
  genome(mySession) = "hg19"

  rmskTrack = ucscTableQuery(mySession, track="rmsk")
  chrs = paste("chr", c(1:22,"X","Y"), sep="")
  repeats = GRanges()

  if(verbose==TRUE) pb = txtProgressBar(min = 0, max = length(chrs), style = 3)

  for(i in 1:length(chrs)){
    c = chrs[i]
    rmskTrackChr = rmskTrack
    range(rmskTrackChr) = range(rmskTrackChr)[seqnames(range(rmskTrackChr)) == c]
    repeats_df = getTable(rmskTrackChr)
    repeats = c(repeats, GRanges(IRanges(repeats_df$genoStart, repeats_df$genoEnd), seqnames=factor(repeats_df$genoName, levels=chrs), type=as.character(repeats_df$repFamily)))
    
    if(verbose==TRUE) setTxtProgressBar(pb, i)

  }
  if(verbose==TRUE) close(pb)
  
  linesL1 = repeats[repeats$type == "L1"]
  linesL2 = repeats[repeats$type == "L2"]
  sinesAlu = repeats[repeats$type == "Alu"]
  sinesMIR = repeats[repeats$type == "MIR"]  
  ## load previously saved segmental duplications to avoid import from UCSC browser
  if(file.exists(file.path(path.package("RSVSim"), "data", "segmentalDuplications.RData"))){
    data("segmentalDuplications", package="RSVSim", envir=environment())
  }else{
    segDups = .loadFromUCSC_SegDups()
  }
  trs = .loadFromBSGenome_TandemRepeats()
  repeats = list(linesL1, linesL2, sinesAlu, sinesMIR, segDups, trs)
  if(save == TRUE){
    save(repeats, file=file.path(path.package("RSVSim"), "data", "repeats_hg19.RData"), compress="xz")
  }
  return(repeats)
}

.loadFromUCSC_SegDups <- function(){
  chrs = paste("chr", c(1:22,"X","Y"), sep="")

  mySession = browserSession("UCSC")
  genome(mySession) = "hg19"

  segDups = getTable(ucscTableQuery(mySession, track="genomicSuperDups"))
  segDups = segDups[segDups$chrom %in% chrs & segDups$otherChrom %in% chrs, ]
  segDups = c(
    GRanges(IRanges(segDups$chromStart, segDups$chromEnd), seqnames=as.character(segDups$chrom)),
    GRanges(IRanges(segDups$otherStart, segDups$otherEnd), seqnames=as.character(segDups$otherChrom))
    )
  segDups$type = "SD"
  
  return(segDups)
}

.loadFromBSGenome_TandemRepeats <- function(){
  ## this function requires BSgenome.Hsapiens.UCSC.hg19 to be loaded
  chrs = paste("chr", c(1:22,"X","Y"), sep="")
  tr = GRanges()
  seqlevels(tr) = chrs
  for(c in chrs){
    t = as(masks(Hsapiens[[c]])["TRF"], "data.frame")
    t = GRanges(IRanges(t$start, t$end), seqnames=c)
    seqlevels(t) = chrs
    tr = c(tr, t)
  }  
  tr$type = "TR"
  return(tr)
  
}


##########################################################################################################
## Methods for internal use only
##########################################################################################################


.getHG19 <- function(chrs){

  ## use library BSgenome.Hsapiens.UCSC.hg19 as default genome
  require(BSgenome.Hsapiens.UCSC.hg19)
  require(BSgenome.Hsapiens.UCSC.hg19.masked)
  Hsapiens_masked = BSgenome.Hsapiens.UCSC.hg19.masked

  if(any(is.na(chrs))){
    chrs = names(Hsapiens)[1:24]
  }
  genome = DNAStringSet()
  gaps = data.frame()
  for(c in chrs){
    genome = append(genome, DNAStringSet(DNAString(Hsapiens[[c]])))
    g = as(masks(Hsapiens_masked[[c]])["AGAPS"], "data.frame")
    g$seqnames = c
    gaps = rbind(gaps, g)
  }
  names(genome) = chrs
  gaps = GRanges(IRanges(gaps$start, gaps$end), seqnames=gaps$seqnames)
  return(list(genome, gaps))
}

.getDummyDataframe <- function(){
  return(data.frame(seqnames=0,start=0,end=0,mechanism="", bpRegion="", stringsAsFactors=FALSE)[-1, ])
}

## subtract interval list from other interval list (both GRanges objects)
## e.g. subtract gaps from genome ranges
.subtractIntervals <- function(subject, subtrahend){
  
#  diff = setdiff(disjoin(c(subject,subtrahend)), subtrahend)
  diff = setdiff(subject, subtrahend)
  seqlevels(diff) = unique(as.character(seqnames(diff)))
  return(diff)
  
}

subject = GRanges(IRanges(c(1,11,21),c(10,20,30)), seqnames=c("A","B","C"), type=c("t1","t2","t3"))
names(subject) = c("t1","t2")
subtrahend = GRanges(IRanges(c(4,6),c(14,16)), seqnames=c("A","B"), type=c("t4","t5"))
names(subtrahend) = c("t3","t4")

## These functions are only for internal purposes
.testSVSim <- function(){

  runs = 10
  
  for(i in 1:runs){
    message("+++++++++++++++ Test run ", i, " ++++++++++++++++++++++++++++++")
  
    size1 = 20000000
    size2 = 15000000
    size3 = 10000000
    size4 = 5000000
    size5 = 1000000
    size6 = 500000
    genome = DNAStringSet(c(
      paste(sample(c("A","G","T","C"), size1, replace=TRUE), collapse=""),
      paste(sample(c("A","G","T","C"), size2, replace=TRUE), collapse=""),
      paste(sample(c("A","G","T","C"), size3, replace=TRUE), collapse=""),
      paste(sample(c("A","G","T","C"), size4, replace=TRUE), collapse=""),
      paste(sample(c("A","G","T","C"), size5, replace=TRUE), collapse=""),
      paste(sample(c("A","G","T","C"), size6, replace=TRUE), collapse="")
      ))
    names(genome) = c("chr1","chr2","chr3","chr4","chr5","chr6")
    
    dels = 100
    ins = 100
    invs = 100
    dups = 100
    trans = 5
    sizeDels = sample(1000:5000, dels, replace=TRUE)
    sizeIns = sample(1000:5000, ins, replace=TRUE)
    sizeInvs = sample(1000:5000, invs, replace=TRUE)
    sizeDups = sample(1000:5000, dups, replace=TRUE)
    bpSeqSize = 300
    sim = simulateSV(output=NA, genome=genome, dels=dels, ins=ins, invs=invs, dups=dups, trans=trans, sizeDels=sizeDels, sizeIns=sizeIns, sizeInvs=sizeInvs, sizeDups=sizeDups, maxDups=10, bpSeqSize=bpSeqSize, random=TRUE, repeatBias=FALSE, percSNPs=0, indelProb=0)

    ## sim2 is for testing the non-random implementation of predefined sv regions
    regionsDels = GRanges(IRanges(metadata(sim)$deletions$Start, metadata(sim)$deletions$End), seqnames=metadata(sim)$deletions$Chr)
    regionsIns = GRanges(IRanges(metadata(sim)$insertions$StartA, metadata(sim)$insertions$EndA), seqnames=metadata(sim)$insertions$ChrA, chrB=metadata(sim)$insertions$ChrB, startB=metadata(sim)$insertions$StartB, endB=metadata(sim)$insertions$EndB)
    regionsInvs = GRanges(IRanges(metadata(sim)$inversions$Start, metadata(sim)$inversions$End), seqnames=metadata(sim)$inversions$Chr)
    regionsDups = GRanges(IRanges(metadata(sim)$tandemDuplications$Start, metadata(sim)$tandemDuplications$End), seqnames=metadata(sim)$tandemDuplications$Chr)
    regionsTrans = GRanges(IRanges(metadata(sim)$translocations$StartA, metadata(sim)$translocations$EndA), seqnames=metadata(sim)$translocations$ChrA, chrB=metadata(sim)$translocations$ChrB, startB=metadata(sim)$translocations$StartB, endB=metadata(sim)$translocations$EndB)
    sim2 = simulateSV(output=NA, genome=genome, dels=dels, ins=ins, invs=invs, dups=dups, trans=trans, sizeDels=sizeDels, sizeIns=sizeIns, sizeInvs=sizeInvs, sizeDups=sizeDups, regionsDels=regionsDels, regionsIns=regionsIns, regionsInvs=regionsInvs, regionsDups=regionsDups, regionsTrans=regionsTrans, maxDups=10, bpSeqSize=bpSeqSize, random=FALSE, repeatBias=FALSE, percSNPs=0, indelProb=0)
    sim = sim2
    
    bpSeqSize = bpSeqSize / 2
    maxMismatch = 0
    for(type in c("deletions","insertions","inversions","tandemDuplications","translocations")){
      sv = metadata(sim)[[type]]
      
      ## Deletions ################################
      if(type == "deletions"){
        iscorrect = c()
        for(i in 1:nrow(sv)){
          seq = DNAString(sv$BpSeq[i])
          aln = matchPattern(subseq(seq, 1, bpSeqSize), genome[[as.character(sv$Chr[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          s = (start(aln) == (sv$Start[i] - bpSeqSize))
          aln = matchPattern(subseq(seq, bpSeqSize+1, bpSeqSize*2), genome[[as.character(sv$Chr[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          e = (end(aln) == (sv$End[i] + bpSeqSize))
          iscorrect = c(iscorrect, s & e)
        }
        message(type, ": ", round(sum(iscorrect)/nrow(sv), digits=2)*100, " % correct ")
        if(any(!iscorrect)){
          message(type, paste(sv$Type[!iscorrect], collapse=","), " are not correct")
        }
      }
      ## Inversions ################################
      if(type == "inversions"){
        iscorrect = c()
        for(i in 1:nrow(sv)){
          
          ## 5' breakpoint sequence
          seq = DNAString(sv$BpSeq_5prime[i])
          aln = matchPattern(subseq(seq, 1, bpSeqSize), genome[[as.character(sv$Chr[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          s = (start(aln) == (sv$Start[i] - bpSeqSize))
          aln = matchPattern(reverseComplement(subseq(seq, bpSeqSize+1, bpSeqSize*2)), genome[[as.character(sv$Chr[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          e = (end(aln) == sv$End[i])
          iscorrect = c(iscorrect, s & e)
          
          ## 3' breakpoint sequence
          seq = DNAString(sv$BpSeq_3prime[i])
          aln = matchPattern(reverseComplement(subseq(seq, 1, bpSeqSize)), genome[[as.character(sv$Chr[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          s = (start(aln) == sv$Start[i])
          aln = matchPattern(subseq(seq, bpSeqSize+1, bpSeqSize*2), genome[[as.character(sv$Chr[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          e = (end(aln) == (sv$End[i]+ bpSeqSize))
          iscorrect = c(iscorrect, s & e)
        }
        message(type, ": ", round(sum(iscorrect)/(nrow(sv)*2), digits=2)*100, " % correct ")
        if(any(!iscorrect)){
          message(type, paste(sv$Type[!iscorrect], collapse=","), " are not correct")
        }
      }
      
      ## Insertions ################################
      if(type == "insertions"){
        iscorrect = c()
        for(i in 1:nrow(sv)){
          
          ## breakpoint sequence A (deleted segment)
          seq = DNAString(sv$BpSeqA[i])
          aln = matchPattern(subseq(seq, 1, bpSeqSize), genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          s = (start(aln) == (sv$StartA[i] - bpSeqSize))
          aln = matchPattern(subseq(seq, bpSeqSize+1, bpSeqSize*2), genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          e = (end(aln) == (sv$EndA[i] + bpSeqSize))
          iscorrect = c(iscorrect, s & e)
          
          ## breakpoint sequence B 5'
          seq = DNAString(sv$BpSeqB_5prime[i])
          aln = matchPattern(subseq(seq, 1, bpSeqSize), genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          s = (start(aln) == (sv$StartB[i] - bpSeqSize))
          aln = matchPattern(subseq(seq, bpSeqSize+1, bpSeqSize*2), genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          e = (start(aln) == sv$StartA[i])
          iscorrect = c(iscorrect, s & e)
          
          ## breakpoint sequence B 3'
          seq = DNAString(sv$BpSeqB_3prime[i])
          aln = matchPattern(subseq(seq, 1, bpSeqSize), genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          s = (end(aln) == sv$EndA[i])
          aln = matchPattern(subseq(seq, bpSeqSize+1, bpSeqSize*2), genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          e = (start(aln) == sv$StartB[i])
          
          iscorrect = c(iscorrect, s & e)
          
        }
        message(type, ": ", round(sum(iscorrect)/(nrow(sv)*3), digits=2)*100, " % correct ")
        if(any(!iscorrect)){
          message(type, paste(sv$Type[!iscorrect], collapse=","), " are not correct")
        }
      }    
      ## Tandem duplications ################################
      if(type == "tandemDuplications"){
        iscorrect = c()
        for(i in 1:nrow(sv)){
          seq = DNAString(sv$BpSeq[i])
          aln = matchPattern(subseq(seq, 1, bpSeqSize), genome[[as.character(sv$Chr[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          s = (end(aln) == sv$End[i])
          aln = matchPattern(subseq(seq, bpSeqSize+1, bpSeqSize*2), genome[[as.character(sv$Chr[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          e = (start(aln) == sv$Start[i])
          ## also test the number of duplication
          aln = matchPattern(seq, sim[[as.character(sv$Chr[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
          d = ((length(aln)+1) == sv$Duplications[i])
          iscorrect = c(iscorrect, s & e & d)
        }
        message(type, ": ", round(sum(iscorrect)/nrow(sv), digits=2)*100, " % correct ")
        if(any(!iscorrect)){
          message(type, paste(sv$Type[!iscorrect], collapse=","), " are not correct")
        }
      }
      ## Translocations ################################
      if(type == "translocations"){
        iscorrect = c()
        for(i in 1:nrow(sv)){
          
          ## breakpoint sequence A
          seq = DNAString(sv$BpSeqA[i])
          ## case xxx... <-> ...xxx
          if(sv$StartA[i] == 1 & sv$StartB[i] != 1){
            seq1 = reverseComplement(subseq(seq, 1, bpSeqSize))
            aln = matchPattern(seq1, genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
            s = (start(aln) == sv$StartB[i])
            seq2 = subseq(seq, bpSeqSize+1, bpSeqSize*2)
            aln = matchPattern(seq2, genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
            e = (end(aln) == sv$EndA[i] + bpSeqSize)
          }
          ## case ...xxx <-> xxx...
          if(sv$StartA[i] != 1 & sv$StartB[i] == 1){
            seq1 = subseq(seq, 1, bpSeqSize)
            aln = matchPattern(seq1, genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
            s = (start(aln) == sv$StartA[i] - bpSeqSize)
            seq2 = reverseComplement(subseq(seq, bpSeqSize+1, bpSeqSize*2))
            aln = matchPattern(seq2, genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
            e = (end(aln) == sv$EndB[i])
          }
          ## case xxx... <-> xxx...
          if(sv$StartA[i] == 1 & sv$StartB[i] == 1){
            seq1 = subseq(seq, 1, bpSeqSize)
            aln = matchPattern(seq1, genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
            s = (end(aln) == sv$EndB[i])
            seq2 = subseq(seq, bpSeqSize+1, bpSeqSize*2)
            aln = matchPattern(seq2, genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
            e = (end(aln) == (sv$EndA[i] + bpSeqSize))
          }
          ## case ...xxx <-> ...xxx
          if(sv$StartA[i] != 1 & sv$StartB[i] != 1){
            seq1 = subseq(seq, 1, bpSeqSize)
            aln = matchPattern(seq1, genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
            s = (start(aln) == sv$StartA[i] - bpSeqSize)
            seq2 = subseq(seq, bpSeqSize+1, bpSeqSize*2)
            aln = matchPattern(seq2, genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
            e = (start(aln) == sv$StartB[i])
          }
          iscorrect = c(iscorrect, s & e)
          
          ## breakpoint sequence B
          seq = DNAString(sv$BpSeqB[i])
          if(sv$Balanced[i] == TRUE){
            ## case xxx... <-> ...xxx
            if(sv$StartA[i] == 1 & sv$StartB[i] != 1){
              seq1 = subseq(seq, 1, bpSeqSize)
              aln = matchPattern(seq1, genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
              s = (start(aln) == sv$StartB[i] - bpSeqSize)
              seq2 = reverseComplement(subseq(seq, bpSeqSize+1, bpSeqSize*2))
              aln = matchPattern(seq2, genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
              e = (end(aln) == sv$EndA[i])
            }
            ## case ...xxx <-> xxx...
            if(sv$StartA[i] != 1 & sv$StartB[i] == 1){
              seq1 = reverseComplement(subseq(seq, 1, bpSeqSize))
              aln = matchPattern(seq1, genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
              s = (start(aln) == sv$StartA[i])
              seq2 = subseq(seq, bpSeqSize+1, bpSeqSize*2)
              aln = matchPattern(seq2, genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
              e = (end(aln) == sv$EndB[i] + bpSeqSize)
            }
            ## case xxx... <-> xxx...
            if(sv$StartA[i] == 1 & sv$StartB[i] == 1){
              seq1 = subseq(seq, 1, bpSeqSize)
              aln = matchPattern(seq1, genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
              s = (end(aln) == sv$EndA[i])
              seq2 = subseq(seq, bpSeqSize+1, bpSeqSize*2)
              aln = matchPattern(seq2, genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
              e = (end(aln) == (sv$EndB[i] + bpSeqSize))
            }
            ## case ...xxx <-> ...xxx
            if(sv$StartA[i] != 1 & sv$StartB[i] != 1){
              seq1 = subseq(seq, 1, bpSeqSize)
              aln = matchPattern(seq1, genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
              s = (start(aln) == sv$StartB[i] - bpSeqSize)
              seq2 = subseq(seq, bpSeqSize+1, bpSeqSize*2)
              aln = matchPattern(seq2, genome[[as.character(sv$ChrA[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
              e = (start(aln) == sv$StartA[i])
            }
            iscorrect = c(iscorrect, s & e)
          }else{
            seq = DNAString(sv$BpSeqB[i])
            aln = matchPattern(seq, genome[[as.character(sv$ChrB[i])]], with.indels=TRUE, max.mismatch=maxMismatch)
            s = (start(aln) == (sv$StartB[i] - bpSeqSize)) | (end(aln) == (sv$EndB[i] + bpSeqSize))
            iscorrect = c(iscorrect, s)
          }        
        }
        message(type, ": ", round(sum(iscorrect)/(nrow(sv)*2), digits=2)*100, " % correct ")
        if(any(!iscorrect)){
          message(type, paste(sv$Type[!iscorrect], collapse=","), " are not correct")
        }
      }
    }   
  }
}
