.addBreakpointMutations <- function(region, percSNPs, indelProb, maxIndelSize){

  regionSize = length(region)
  adjustBy = 0

  if(regionSize > 3 & percSNPs > 0){
    ## place SNPs
    idx = sample(x=c(TRUE,FALSE), size=regionSize, replace=TRUE, prob=c(percSNPs, 1-percSNPs))
    region = replaceLetterAt(region, idx, sample(x=c("A","C","G","T"), size=sum(idx), replace=TRUE))
  }
  
  if(regionSize > 3 & indelProb > 0){
    ## place maximum one indel into each region (50/50 if deletion or insertion)
    type = sample(c("del","ins"), 1)
    doIndel = sample(c(TRUE, FALSE), 1, prob=c(indelProb, 1-indelProb))
    if(doIndel & type == "del"){
      start = sample(1:(regionSize-maxIndelSize), 1)
      end = start + sample(1:maxIndelSize, 1)-1
      region = DNAString(paste(subseq(region, 1, start-1), subseq(region, end+1, length(region)), sep=""))
      adjustBy = adjustBy + (end-start+1)*(-1)
    }
    if(doIndel & type == "ins"){
      pos = sample(1:regionSize, 1)
      numNuc = sample(1:maxIndelSize, 1)
      insSeq = paste(sample(c("A","T","G","C"), numNuc, replace=TRUE), collapse="")
      region = DNAString(paste(subseq(region, 1, pos), insSeq, subseq(region, pos+1, length(region)), sep=""))
      adjustBy = adjustBy + numNuc
    }
  }

  return(list(region, adjustBy))

}

.execDeletion <- function(genomeSeq, c, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, bpSeqSize, i, bpFlankSize, percSNPs, indelProb, maxIndelSize){
  
  pos = posDel[i, ]
  start = as.integer(pos$start)
  end = as.integer(pos$end)

  ## add random mutations (SNPs an/or indels) at the breakpoints flanking regions (5' and 3')
  ## take care not to exceed the genome limits
  ## take care to adjust coordinates of other SVs correspondingly
  sf  = max(1, start-bpFlankSize)
  ef = min(length(genomeSeq), end+bpFlankSize)
  flankingRegion1 = subseq(genomeSeq, sf, start-1)
  flankingRegion2 = subseq(genomeSeq, end+1, ef)
  flankingRegion1 = .addBreakpointMutations(flankingRegion1, percSNPs, indelProb, maxIndelSize)
  adjustBy = flankingRegion1[[2]]
  flankingRegion1 = flankingRegion1[[1]]
  flankingRegion2 = .addBreakpointMutations(flankingRegion2, percSNPs, indelProb, maxIndelSize)
  adjustBy = adjustBy + flankingRegion2[[2]]
  flankingRegion2 = flankingRegion2[[1]]
  
  ## paste genome: fist part of chromosome - flanking bp region (5') - deleted region - flanking bp region (3') - rest of the chromosome
  genomeSeq = DNAString(paste(subseq(genomeSeq, 1, sf-1), flankingRegion1, flankingRegion2, subseq(genomeSeq, ef+1, length(genomeSeq)), sep=""))
  
  adjustBy = adjustBy + (end-start+1)*-1 # adjust positions by deletion size
  posDel = .adjustPositions(posDel, c, end, adjustBy)
  posIns_1= .adjustPositions(posIns_1, c, end, adjustBy)
  posIns_2= .adjustPositions(posIns_2, c, end, adjustBy)
  posInv = .adjustPositions(posInv, c, end, adjustBy)
  posDup = .adjustPositions(posDup, c, end, adjustBy)
  posTrans_1 = .adjustPositions(posTrans_1, c, end, adjustBy)
  posTrans_2 = .adjustPositions(posTrans_2, c, end, adjustBy)

  return(list(genomeSeq, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2))
}

.execInsertion <- function(genomeSeqA, genomeSeqB, c1, c2, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, i, bpSeqSize, bpFlankSize, percSNPs, indelProb, maxIndelSize){

  ## chrA
  pos = posIns_1[i, ]
  startA = as.integer(pos$start)
  endA = as.integer(pos$end)
  copied = pos$copied
  insSeq = subseq(genomeSeqA, startA, endA)
  
  ## remove or copy sequence on chrA
  if(copied == TRUE){
    ## nothing to do in this case
  }else{
    ## add random mutations (SNPs an/or indels) at the breakpoints flanking regions (5' and 3')
    ## take care not to exceed the genome limits
    ## take care to adjust coordinates of other SVs correspondingly
    sf  = max(1, startA-bpFlankSize)
    ef = min(length(genomeSeqA), endA+bpFlankSize)
    flankingRegion1 = subseq(genomeSeqA, sf, startA-1)
    flankingRegion2 = subseq(genomeSeqA, endA+1, ef)
    flankingRegion1 = .addBreakpointMutations(flankingRegion1, percSNPs, indelProb, maxIndelSize)
    adjustBy = flankingRegion1[[2]]
    flankingRegion1 = flankingRegion1[[1]]
    flankingRegion2 = .addBreakpointMutations(flankingRegion2, percSNPs, indelProb, maxIndelSize)
    adjustBy = adjustBy + flankingRegion2[[2]]
    flankingRegion2 = flankingRegion2[[1]]
    
    genomeSeqA = DNAString(paste(subseq(genomeSeqA, 1, sf-1), flankingRegion1, flankingRegion2, subseq(genomeSeqA, ef+1, length(genomeSeqA)), sep=""))
    ## adjust coordinates downstream of insertion
    ## (take care, that the end of ins1 is adjusted as well!)
    adjustBy = adjustBy + (endA-startA+1) * -1
    posIns_1$end[i] = posIns_1$start[i]
    posDel = .adjustPositions(posDel, c1, endA, adjustBy)
    posIns_1 = .adjustPositions(posIns_1, c1, endA, adjustBy)
    posIns_2 = .adjustPositions(posIns_2, c1, endA, adjustBy)
    posInv = .adjustPositions(posInv, c1, endA, adjustBy)
    posDup = .adjustPositions(posDup, c1, endA, adjustBy)
    posTrans_1 = .adjustPositions(posTrans_1, c1, endA, adjustBy)
    posTrans_2 = .adjustPositions(posTrans_2, c1, endA, adjustBy)
    ## update genome if insertion happens on the same chromosome
    if(c1 == c2){
     genomeSeqB = genomeSeqA
    }

  }  
  
  ## insert sequence from chrA into chrB
  pos = posIns_2[i, ]
  startB = as.integer(pos$start)
  endB = as.integer(pos$end)

  ## add random mutations (SNPs an/or indels) at the breakpoints flanking regions (5' and 3')
  ## take care not to exceed the genome limits
  ## take care to adjust coordinates of other SVs correspondingly
  sf  = max(1, startB-bpFlankSize)
  ef = min(length(genomeSeqB), startB+bpFlankSize-1)
  flankingRegion1 = subseq(genomeSeqB, sf, startB-1)
  flankingRegion2 = subseq(genomeSeqB, startB, ef)
  flankingRegion1 = .addBreakpointMutations(flankingRegion1, percSNPs, indelProb, maxIndelSize)
  adjustBy = flankingRegion1[[2]]
  flankingRegion1 = flankingRegion1[[1]]
  flankingRegion2 = .addBreakpointMutations(flankingRegion2, percSNPs, indelProb, maxIndelSize)
  adjustBy = adjustBy + flankingRegion2[[2]]
  flankingRegion2 = flankingRegion2[[1]]
  
  genomeSeqB = DNAString(paste(subseq(genomeSeqB, 1, sf-1),flankingRegion1, insSeq, flankingRegion2, subseq(genomeSeqB, ef+1, length(genomeSeqB)), sep=""))
    
  ## adjust coordinates downstream of insertion
  adjustBy = adjustBy + endA-startA+1
  posDel = .adjustPositions(posDel, c2, startB, adjustBy)
  posIns_1 = .adjustPositions(posIns_1, c2, startB, adjustBy)
  posIns_2 = .adjustPositions(posIns_2, c2, startB, adjustBy)
  posInv = .adjustPositions(posInv, c2, startB, adjustBy)
  posDup = .adjustPositions(posDup, c2, startB, adjustBy)
  posTrans_1 = .adjustPositions(posTrans_1, c2, startB, adjustBy)
  posTrans_2 = .adjustPositions(posTrans_2, c2, startB, adjustBy)
  ## update genome if insertion happens on the same chromosome
  if(c1 == c2){
    genomeSeqA = genomeSeqB
  }

  return(list(genomeSeqA, genomeSeqB, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2))
}

.execInversion <- function(genomeSeq, c, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, bpSeqSize, i, bpFlankSize, percSNPs, indelProb, maxIndelSize){
  
  pos = posInv[i, ]
  start = as.integer(pos$start)
  end = as.integer(pos$end)

  ## add random mutations (SNPs an/or indels) at the breakpoints flanking regions (5' and 3')
  ## take care not to exceed the genome limits
  ## take care to adjust coordinates of other SVs correspondingly
  sf  = max(1, start-bpFlankSize)
  ef = min(length(genomeSeq), end+bpFlankSize)
  flankingRegion1 = subseq(genomeSeq, sf, start-1)
  flankingRegion2 = subseq(genomeSeq, end+1, ef)
  flankingRegion1 = .addBreakpointMutations(flankingRegion1, percSNPs, indelProb, maxIndelSize)
  adjustBy = flankingRegion1[[2]]
  flankingRegion1 = flankingRegion1[[1]]
  flankingRegion2 = .addBreakpointMutations(flankingRegion2, percSNPs, indelProb, maxIndelSize)
  adjustBy = adjustBy + flankingRegion2[[2]]
  flankingRegion2 = flankingRegion2[[1]]

  invertedSeq = reverseComplement(subseq(genomeSeq, start, end))
  genomeSeq = DNAString(paste(subseq(genomeSeq, 1, sf-1), flankingRegion1, invertedSeq, flankingRegion2, subseq(genomeSeq, ef+1, length(genomeSeq)), sep=""))

  posDel = .adjustPositions(posDel, c, end, adjustBy)
  posIns_1= .adjustPositions(posIns_1, c, end, adjustBy)
  posIns_2= .adjustPositions(posIns_2, c, end, adjustBy)
  posInv = .adjustPositions(posInv, c, end, adjustBy)
  posDup = .adjustPositions(posDup, c, end, adjustBy)
  posTrans_1 = .adjustPositions(posTrans_1, c, end, adjustBy)
  posTrans_2 = .adjustPositions(posTrans_2, c, end, adjustBy)

  return(list(genomeSeq, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2))
}

.execTandemDuplication <- function(genomeSeq, c, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, bpSeqSize, times, i, bpFlankSize, percSNPs, indelProb, maxIndelSize){
  
  pos = posDup[i, ]
  start = as.integer(pos$start)
  end = as.integer(pos$end)

  ## add random mutations (SNPs an/or indels) at the breakpoints flanking regions (5' and 3')
  ## take care not to exceed the genome limits
  ## take care to adjust coordinates of other SVs correspondingly
  sf  = max(1, start-bpFlankSize)
  ef = min(length(genomeSeq), end+bpFlankSize)
  flankingRegion1 = subseq(genomeSeq, sf, start-1)
  flankingRegion2 = subseq(genomeSeq, end+1, ef)
  flankingRegion1 = .addBreakpointMutations(flankingRegion1, percSNPs, indelProb, maxIndelSize)
  adjustBy = flankingRegion1[[2]]
  flankingRegion1 = flankingRegion1[[1]]
  flankingRegion2 = .addBreakpointMutations(flankingRegion2, percSNPs, indelProb, maxIndelSize)
  adjustBy = adjustBy + flankingRegion2[[2]]
  flankingRegion2 = flankingRegion2[[1]]

  s = min(end-start+1, bpSeqSize)
  bpSeq = paste(subseq(genomeSeq, end-s+1, end), subseq(genomeSeq, start, start+s-1), sep="")
  dupSeq = subseq(genomeSeq, start, end)
  dupSeq = paste(rep(dupSeq, times + 1), collapse="")
  genomeSeq = DNAString(paste(subseq(genomeSeq, 1, sf-1), flankingRegion1, dupSeq, flankingRegion2, subseq(genomeSeq, ef+1, length(genomeSeq)), sep=""))
  
  adjustBy = adjustBy + (end-start+1)*(times-1) # adjust positions by length of additional sequence
  posDel = .adjustPositions(posDel, c, end, adjustBy)
  posIns_1= .adjustPositions(posIns_1, c, end, adjustBy)
  posIns_2= .adjustPositions(posIns_2, c, end, adjustBy)
  posInv = .adjustPositions(posInv, c, end, adjustBy)
  posDup = .adjustPositions(posDup, c, end, adjustBy)
  posTrans_1 = .adjustPositions(posTrans_1, c, end, adjustBy)
  posTrans_2 = .adjustPositions(posTrans_2, c, end, adjustBy)
  
  return(list(genomeSeq, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, bpSeq))
}

.execTranslocation <- function(genomeSeqA, genomeSeqB, c1, c2, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, i, bpSeqSize, bpFlankSize, percSNPs, indelProb, maxIndelSize){
  
  ## chrA
  pos = posTrans_1[i, ]
  startA = as.integer(pos$start)
  endA = as.integer(pos$end)
  balanced = pos$balanced
  inverted = pos$inverted
  ## chrB
  pos = posTrans_2[i, ]
  startB = as.integer(pos$start)
  endB = as.integer(pos$end)

  ## eventually invert the translocated sequence (only when exchange between different ends; 5'<->3' or 3'<->5')
  if(inverted == TRUE){
    transSeqA = reverseComplement(subseq(genomeSeqA, startA, endA))
    transSeqB = reverseComplement(subseq(genomeSeqB, startB, endB))
  }else{
    transSeqA = subseq(genomeSeqA, startA, endA)
    transSeqB = subseq(genomeSeqB, startB, endB)
  }
  
  ## chrA ----> chrB
  
  ## add random mutations (SNPs an/or indels) at the breakpoints flanking regions (5' and 3')
  ## take care not to exceed the genome limits
  ## take care to adjust coordinates of other SVs correspondingly
  sf  = max(1, startB-bpFlankSize)
  ef = min(length(genomeSeqB), endB+bpFlankSize)
  flankingRegion1 = subseq(genomeSeqB, sf, startB-1)
  flankingRegion2 = subseq(genomeSeqB, endB+1, ef)
  flankingRegion1 = .addBreakpointMutations(flankingRegion1, percSNPs, indelProb, maxIndelSize)
  adjustBy = flankingRegion1[[2]]
  flankingRegion1 = flankingRegion1[[1]]
  flankingRegion2 = .addBreakpointMutations(flankingRegion2, percSNPs, indelProb, maxIndelSize)
  adjustBy = adjustBy + flankingRegion2[[2]]
  flankingRegion2 = flankingRegion2[[1]]
 
  genomeSeqB = DNAString(paste(subseq(genomeSeqB, 1, sf-1), flankingRegion1, transSeqA, flankingRegion2, subseq(genomeSeqB, ef+1, length(genomeSeqB)), sep=""))  
  
  ## adjust coordinates downstream of translocation
  ## (take care, that the end of trans1 is adjusted as well!)
  adjustBy = adjustBy + (endA-startA+1) - (endB-startB+1)
  posTrans_2$end[i] = posTrans_2$end[i] + adjustBy
  posDel = .adjustPositions(posDel, c2, endB, adjustBy)
  posIns_1 = .adjustPositions(posIns_1, c2, endB, adjustBy)
  posIns_2 = .adjustPositions(posIns_2, c2, endB, adjustBy)
  posInv = .adjustPositions(posInv, c2, endB, adjustBy)
  posDup = .adjustPositions(posDup, c2, endB, adjustBy)
  posTrans_1 = .adjustPositions(posTrans_1, c2, endB, adjustBy)
  posTrans_2 = .adjustPositions(posTrans_2, c2, endB, adjustBy)
  
  ## if balanced: chrA <---- chrB
  startA = as.integer(posTrans_1$start[i])
  endA = as.integer(posTrans_1$end[i])
 
  if(balanced == TRUE){

    ## add random mutations (SNPs an/or indels) at the breakpoints flanking regions (5' and 3')
    ## take care not to exceed the genome limits
    ## take care to adjust coordinates of other SVs correspondingly
    sf  = max(1, startA-bpFlankSize)
    ef = min(length(genomeSeqA), endA+bpFlankSize)
    flankingRegion1 = subseq(genomeSeqA, sf, startA-1)
    flankingRegion2 = subseq(genomeSeqA, endA+1, ef)
    flankingRegion1 = .addBreakpointMutations(flankingRegion1, percSNPs, indelProb, maxIndelSize)
    adjustBy = flankingRegion1[[2]]
    flankingRegion1 = flankingRegion1[[1]]
    flankingRegion2 = .addBreakpointMutations(flankingRegion2, percSNPs, indelProb, maxIndelSize)
    adjustBy = adjustBy + flankingRegion2[[2]]
    flankingRegion2 = flankingRegion2[[1]]
    
    genomeSeqA = DNAString(paste(subseq(genomeSeqA, 1, sf-1), flankingRegion1, transSeqB, flankingRegion2, subseq(genomeSeqA, ef+1, length(genomeSeqA)), sep=""))
    
    ## adjust coordinates downstream of balanced translocation
    adjustBy = adjustBy + (endB-startB+1) - (endA-startA+1)
  ## if not balanced: chrA <--NOT-- chrB
  }else{
    ## code for deleting the segment
    ## genomeSeqA = DNAString(paste(subseq(genomeSeqA, 1, startA-1), subseq(genomeSeqA, endA+1, length(genomeSeqA)), sep=""))    
    ## adjustBy = (endA-startA+1) * -1
    adjustBy = 0
  }
  ## adjust coordinates downstream of balanced translocation
  ## (take care, that the end of trans2 is adjusted as well!)
  posTrans_1$end[i] = posTrans_1$end[i] + adjustBy
  posDel = .adjustPositions(posDel, c1, endA, adjustBy)
  posIns_1 = .adjustPositions(posIns_1, c1, endA, adjustBy)
  posIns_2 = .adjustPositions(posIns_2, c1, endA, adjustBy)
  posInv = .adjustPositions(posInv, c1, endA, adjustBy)
  posDup = .adjustPositions(posDup, c1, endA, adjustBy)
  posTrans_1 = .adjustPositions(posTrans_1, c1, endA, adjustBy)
  posTrans_2 = .adjustPositions(posTrans_2, c1, endA, adjustBy)
  
  return(list(genomeSeqA, genomeSeqB, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2))
}




## misc ################################################################

.adjustPositions <- function(pos, chr, beyondPos, adjustBy){
  idx = (pos$seqnames == chr) & (pos$start > beyondPos)
  if(any(idx)){
    pos$start[idx] = pos$start[idx] + adjustBy
    pos$end[idx] = pos$end[idx] + adjustBy
  }
  return(pos)
}
  
.getBpSeq <- function(genome, pos, bpSeqSize, i){
  pos = pos[i, ]
  c = as.character(pos$seqnames)
  start = as.integer(pos$start)
  end = as.integer(pos$end)
  seq =  genome[[c]]
  bpSeq = vector(mode="character", length=2)
  ## no 5' breakpoint sequences for svs at the beginning of the genome (terminal translocations)
  if(start == 1){
    bpSeq[1] = ""
  }else{
    s = max(1, start-bpSeqSize) ## take care that breakpoint sequence does not exceed the genome boundaries
    e = min(start+bpSeqSize-1, length(seq)) ## dito
    bpSeq[1] = as.character(subseq(seq, s, e))
  }
  ## no 3' breakpoint sequences for svs at the end of the genome (terminal translocations)
  if(end == length(seq)){
    bpSeq[2] = ""
  }else{
    s = max(1, end-bpSeqSize+1) ## take care that breakpoint sequence does not exceed the genome boundaries
    e = min(end+bpSeqSize, length(seq)) ## dito
    bpSeq[2] = as.character(subseq(seq, s, e))
  }
  return(bpSeq)
}
  




