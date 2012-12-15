.execDeletion <- function(genomeSeq, c, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, bpSeqSize, i){
  
  pos = posDel[i, ]
  start = as.integer(pos$start)
  end = as.integer(pos$end)
  genomeSeq = DNAString(paste(subseq(genomeSeq, 1, start-1), subseq(genomeSeq, end+1, length(genomeSeq)), sep=""))
  
  adjustBy = (end-start+1)*-1 # adjust positions by deletion size
  posDel = .adjustPositions(posDel, c, end, adjustBy)
  posIns_1= .adjustPositions(posIns_1, c, end, adjustBy)
  posIns_2= .adjustPositions(posIns_2, c, end, adjustBy)
  posInv = .adjustPositions(posInv, c, end, adjustBy)
  posDup = .adjustPositions(posDup, c, end, adjustBy)
  posTrans_1 = .adjustPositions(posTrans_1, c, end, adjustBy)
  posTrans_2 = .adjustPositions(posTrans_2, c, end, adjustBy)

  return(list(genomeSeq, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2))
}

.execInsertion <- function(genomeSeqA, genomeSeqB, c1, c2, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, i, bpSeqSize){
  
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
    genomeSeqA = DNAString(paste(subseq(genomeSeqA, 1, startA-1), subseq(genomeSeqA, endA+1, length(genomeSeqA)), sep=""))
    ## adjust coordinates downstream of insertion
    ## (take care, that the end of ins1 is adjusted as well!)
    adjustBy = (endA-startA+1) * -1
    posIns_1$end[i] = posIns_1$start[i]
    posDel = .adjustPositions(posDel, c1, endA, adjustBy)
    posIns_1 = .adjustPositions(posIns_1, c1, endA, adjustBy)
    posIns_2 = .adjustPositions(posIns_2, c1, endA, adjustBy)
    posInv = .adjustPositions(posInv, c1, endA, adjustBy)
    posDup = .adjustPositions(posDup, c1, endA, adjustBy)
    posTrans_1 = .adjustPositions(posTrans_1, c1, endA, adjustBy)
    posTrans_2 = .adjustPositions(posTrans_2, c1, endA, adjustBy)
  }  
  
  ## insert sequence from chrA intro chrB
  pos = posIns_2[i, ]
  startB = as.integer(pos$start)
  endB = as.integer(pos$end)
  genomeSeqB = DNAString(paste(subseq(genomeSeqB, 1, startB-1), insSeq, subseq(genomeSeqB, startB, length(genomeSeqB)), sep=""))
    
  ## adjust coordinates downstream of insertion
  adjustBy = endA-startA+1
  posDel = .adjustPositions(posDel, c2, endB, adjustBy)
  posIns_1 = .adjustPositions(posIns_1, c2, endB, adjustBy)
  posIns_2 = .adjustPositions(posIns_2, c2, endB, adjustBy)
  posInv = .adjustPositions(posInv, c2, endB, adjustBy)
  posDup = .adjustPositions(posDup, c2, endB, adjustBy)
  posTrans_1 = .adjustPositions(posTrans_1, c2, endB, adjustBy)
  posTrans_2 = .adjustPositions(posTrans_2, c2, endB, adjustBy)

  return(list(genomeSeqA, genomeSeqB, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2))
}

.execInversion <- function(genomeSeq, posInv, bpSeqSize, i){
  
  pos = posInv[i, ]
  start = as.integer(pos$start)
  end = as.integer(pos$end)
  invertedSeq = reverseComplement(subseq(genomeSeq, start, end))
  genomeSeq = DNAString(paste(subseq(genomeSeq, 1, start-1), invertedSeq, subseq(genomeSeq, end+1, length(genomeSeq)), sep=""))
  return(genomeSeq)
}

.execTandemDuplication <- function(genomeSeq, c, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, bpSeqSize, times, i){
  
  pos = posDup[i, ]
  start = as.integer(pos$start)
  end = as.integer(pos$end)
  s = min(end-start+1, bpSeqSize)
  bpSeq = paste(subseq(genomeSeq, end-s+1, end), subseq(genomeSeq, start, start+s-1), sep="")
  dupSeq = subseq(genomeSeq, start, end)
  dupSeq = paste(rep(dupSeq, times), collapse="")
  genomeSeq = DNAString(paste(subseq(genomeSeq, 1, start-1), dupSeq, subseq(genomeSeq, end+1, length(genomeSeq)), sep=""))
  
  adjustBy = (end-start+1)*(times-1) # adjust positions by length of additional sequence
  posDel = .adjustPositions(posDel, c, end, adjustBy)
  posIns_1= .adjustPositions(posIns_1, c, end, adjustBy)
  posIns_2= .adjustPositions(posIns_2, c, end, adjustBy)
  posInv = .adjustPositions(posInv, c, end, adjustBy)
  posDup = .adjustPositions(posDup, c, end, adjustBy)
  posTrans_1 = .adjustPositions(posTrans_1, c, end, adjustBy)
  posTrans_2 = .adjustPositions(posTrans_2, c, end, adjustBy)
  
  return(list(genomeSeq, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, bpSeq))
}

.execTranslocation <- function(genomeSeqA, genomeSeqB, c1, c2, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, i, bpSeqSize, transInsert){
  
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
  ## add 0-transInsert random nucleotides at each breakpoint
  if(startB == 1 | transInsert < 1){
    numNuc1 = 0
  }else{
    numNuc1 = sample(1:transInsert, 1)
  }
  bpInsSeq1 = paste(sample(c("A","T","G","C"), numNuc1, replace=TRUE), collapse="")
  if(endB == length(genomeSeqB) | transInsert < 1){
    numNuc2 = 0
  }else{
    numNuc2 = sample(1:transInsert, 1)
  }
  bpInsSeq2 = paste(sample(c("A","T","G","C"), numNuc2, replace=TRUE), collapse="")
  transSeqA = DNAString(paste(bpInsSeq1, transSeqA, bpInsSeq2, sep=""))
  
  genomeSeqB = DNAString(paste(subseq(genomeSeqB, 1, startB-1), transSeqA, subseq(genomeSeqB, endB+1, length(genomeSeqB)), sep=""))  
  
  ## adjust coordinates downstream of translocation
  ## (take care, that the end of trans1 is adjusted as well!)
  adjustBy = (endA-startA+1) - (endB-startB+1) + numNuc1 + numNuc2
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
    ## add 0-transInsert random nucleotides at each breakpoint
    if(startA == 1 | transInsert < 1){
      numNuc1 = 0
    }else{
      numNuc1 = sample(1:transInsert, 1)
    }
    bpInsSeq1 = paste(sample(c("A","T","G","C"), numNuc1, replace=TRUE), collapse="")
    if(endA == length(genomeSeqA) | transInsert < 1){
      numNuc2 = 0
    }else{
      numNuc2 = sample(1:transInsert, 1)
    }
    bpInsSeq2 = paste(sample(c("A","T","G","C"), numNuc2, replace=TRUE), collapse="")
    transSeqB = DNAString(paste(bpInsSeq1, transSeqB, bpInsSeq2, sep=""))

    genomeSeqA = DNAString(paste(subseq(genomeSeqA, 1, startA-1), transSeqB, subseq(genomeSeqA, endA+1, length(genomeSeqA)), sep=""))
    
    ## adjust coordinates downstream of balanced translocation
    adjustBy = (endB-startB+1) - (endA-startA+1) + numNuc1 + numNuc2
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
  




