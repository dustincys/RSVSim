## calculate SV breakpoint coordinates #################################################

.drawPos <- function(bpRegions, size){

  ## randomly select a region (larger regions have higher probability)
  ## regions have a width >= size
  idx = sample(x=1:length(bpRegions), size=1, prob=width(bpRegions))  
  start = sample(x=start(bpRegions[idx]):(end(bpRegions[idx])-size+1), size=1)
  end =  start + size -1
  return(data.frame(seqnames=unique(seqnames(bpRegions)), start=start, end=end))
}


.drawPos_trans <- function(genome, bpRegions, chr){  
 
  ## draw breakpoint position and set the start/end to the outmost coordinate (1 for 5', chromosome end for 3')
  pos = .drawPos(bpRegions, 1) ## here: size = 1 -> start = end
  ## choose start or end whether the start/end of translocation is closer to start/end of the genome
  center = round(length(genome[[chr]]) / 2) + 1
  if(nrow(pos) > 0){
    if(pos$end < center){
      pos$start = 1
    }else{
      pos$end = length(genome[[chr]])
    }
  }
  return(pos)
  
}


## 1. translocations  #####################################################################

.simTranslocationPositions <- function(n, bpRegions, genome, percBalancedTrans, sizes){
  
  posTrans_1 = .getDummyDataframe()
  posTrans_2 = .getDummyDataframe()
#  translocations = data.frame(Type="",ChrA=0,StartA=0,EndA=0,ChrB=0,StartB=0,EndB=0,Balanced=FALSE,BpSeqA_5prime="",BpSeqA_3prime="",BpSeqB_5prime="",BpSeqB_3prime="", stringsAsFactors=FALSE)[-1, ]
  invertedTrans = c()

  if(n > 0){
    for(i in 1:n){
      
      ## first translocation partner
      ## check if regions on two different chromosomes are available
      if(length(unique(seqnames(bpRegions))) < 2){
        stop("No two chromosomes available for translocations")
      }
      ## sample chromosome; larger chromosomes have a higher probability
      chrs = levels(seqnames(bpRegions))
      probs = sapply(chrs, function(x){return(sum(width(bpRegions[seqnames(bpRegions) == x])))})
      chr1 = as.character(sample(chrs, 1, prob=probs))
      pos1 = .drawPos_trans(genome, bpRegions[seqnames(bpRegions) == chr1], chr1)
      posTrans_1 = rbind(posTrans_1, pos1)
      
      ## make sure new SVs do not overlap with this one
      ## to do this, remove the half of the genome where the sv lies on
      pos = GRanges(IRanges(pos1$start, pos1$end), seqnames=pos1$seqnames)
      center = round(length(genome[[chr1]]) / 2) + 1
      if(start(pos) == 1){
        end(pos) = center
      }else{
        start(pos) = center
      }
      bpRegions = .subtractIntervals(bpRegions, pos)

      ## second translocation partner
      chrs = chrs[chrs != chr1] # make sure second translocated segments lies on different chromosome
      probs = sapply(chrs, function(x){return(sum(width(bpRegions[seqnames(bpRegions) == x])))})
      chr2 = as.character(sample(chrs, 1, prob=probs))
      pos2 = .drawPos_trans(genome, bpRegions[seqnames(bpRegions) == chr2], chr2)
      posTrans_2 = rbind(posTrans_2, pos2)
      
      ## make sure new SVs do not overlap with this one
      ## to do this, remove the half of the genome where the sv lies on
      pos = GRanges(IRanges(pos2$start, pos2$end), seqnames=pos2$seqnames)
      center = round(length(genome[[chr2]]) / 2) + 1
      if(start(pos) == 1){
        end(pos) = center
      }else{
        start(pos) = center
      }
      bpRegions = .subtractIntervals(bpRegions, pos)
      
      ## always invert translocated segments between different ends (5'<->3' or 3'<->5')
      if((pos1$start == 1 & pos2$start != 1) | (pos1$start != 1 & pos2$start == 1)){
        invertedTrans = c(invertedTrans, i)
      }
    }
    
    ## randomly select translocations to be balanced
    balancedTrans = sample(1:n, round(n*percBalancedTrans))  ## indices of balanced translocations
    posTrans_1$balanced = FALSE
    posTrans_1$balanced[balancedTrans] = TRUE
    posTrans_2$balanced = posTrans_1$balanced

    posTrans_1$inverted = FALSE
    posTrans_1$inverted[invertedTrans] = TRUE
    posTrans_2$inverted = posTrans_1$inverted
    
    ## requires: balanced is the same for both translocation partners
    size1 = posTrans_1[, 3]-posTrans_1[, 2]+1
    size2 = posTrans_2[, 3]-posTrans_2[, 2]+1
    translocations = cbind("", posTrans_1[, 1:3], size1, posTrans_2[, 1:3], size2, posTrans_2[, "balanced"])
    colnames(translocations) = c("Name", "ChrA", "StartA", "EndA", "SizeA", "ChrB", "StartB", "EndB", "SizeB", "Balanced")
    translocations$Name = paste("translocation_", 1:n, sep="")

  }

  return(list(translocations, posTrans_1, posTrans_2))
}

## 2. insertions  #####################################################################

.simInsertionPositions <- function(n, bpRegions, genome, sizes, percCopiedIns){
  
  posIns_1 = .getDummyDataframe()
  posIns_2 = .getDummyDataframe()
#  insertions = data.frame(Type="",ChrA=0,StartA=0,EndA=0,ChrB=0,StartB=0,EndB=0,BpSeqA_5prime="",BpSeqA_3prime="",BpSeqB_5prime="",BpSeqB_3prime="", stringsAsFactors=FALSE)[-1, ]

  if(n > 0){
    for(i in 1:n){
      
      ## first translocation partner
      bpRegionsGood = bpRegions[width(bpRegions) >= sizes[i]]
      ## check if regions on two different chromosomes are available
      if(length(unique(seqnames(bpRegionsGood))) < 2){
        stop("No two chromosomes available for insertions")
      }
      ## sample chromosome; larger chromosomes have a higher probability
      chrs = levels(seqnames(bpRegions))
      probs = sapply(chrs, function(x){return(sum(width(bpRegions[seqnames(bpRegions) == x])))})
      chr1 = as.character(sample(chrs, 1, prob=probs))
      pos1 = .drawPos(bpRegionsGood[seqnames(bpRegionsGood) == chr1], sizes[i])
      posIns_1 = rbind(posIns_1, pos1)
      
      ## make sure new SVs do not overlap with this one
      bpRegionsGood = .subtractIntervals(bpRegionsGood, GRanges(IRanges(pos1$start, pos1$end), seqnames=pos1$seqnames))
      
      ## second translocation partner
      chrs = chrs[chrs != chr1] # make sure second translocated segments lies on different chromosome
      probs = sapply(chrs, function(x){return(sum(width(bpRegions[seqnames(bpRegions) == x])))})
      chr2 = as.character(sample(chrs, 1, prob=probs))
      pos2 = .drawPos(bpRegionsGood[seqnames(bpRegionsGood) == chr2], sizes[i])
      posIns_2 = rbind(posIns_2, pos2)
      
      ## make sure new SVs do not overlap with this one
      bpRegions = .subtractIntervals(bpRegions, GRanges(IRanges(pos1$start, pos1$end), seqnames=pos1$seqnames))
      bpRegions = .subtractIntervals(bpRegions, GRanges(IRanges(pos2$start, pos2$end), seqnames=pos2$seqnames))
    }
    
    ## randomly select insertions to be copied
    copiedIns = sample(1:n, round(n*percCopiedIns))  ## indices of copied translocations
    posIns_1$copied = FALSE
    posIns_1$copied[copiedIns] = TRUE
    posIns_2$copied = posIns_1$copied
  
    ## requires: inverted and balanced are the same for both translocation partners
    insertions = cbind(paste("insertion_", 1:n, sep=""), posIns_1[, 1:3], posIns_2[, 1:3], sizes, posIns_2[, "copied"])
    colnames(insertions) = c("Name", "ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Size", "Copied")

  }

  return(list(insertions, posIns_1, posIns_2))
}  


## 3.  deletions, inversions and tandem duplications #################################
.simPositions <- function(n, bpRegions, sizes, type){
  
  pos = .getDummyDataframe()
#  svs = data.frame(Type="", Space=0,Start=0,End=0,Size=0,BpSeq="", stringsAsFactors=FALSE)[-1, ]
  if(n > 0){
    for(i in 1:n){
      bpRegionsGood = bpRegions[width(bpRegions) >= sizes[i]]
      if(length(bpRegionsGood) == 0){
        stop("There is no region large enough for ", type, " of size ", sizes[i], "bp")
      }
      ## sample chromosome; larger chromosomes have a higher probability
      chrs = levels(seqnames(bpRegionsGood))
      probs = sapply(chrs, function(x){return(sum(width(bpRegionsGood[seqnames(bpRegionsGood) == x])))})
      chr = sample(chrs, 1, prob=probs)
      p = .drawPos(bpRegionsGood[seqnames(bpRegionsGood) == chr], sizes[i])
      pos = rbind(pos, p)    
      bpRegions = .subtractIntervals(bpRegions, GRanges(IRanges(p$start, p$end), seqnames=p$seqnames))
    }
    
    svs = cbind("", pos, 0)
    colnames(svs) = c("Name", "Chr", "Start","End", "Size")
    svs$Size = sizes
    svs$Name = paste(type, 1:n, sep="")

  }

  return(list(svs, pos))
}
