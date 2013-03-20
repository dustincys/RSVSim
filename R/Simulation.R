## calculate SV breakpoint coordinates #################################################

.drawPos <- function(chr, bpRegionsList, weightsMechanisms, weightsRepeats, size){

  bpRegions = NULL
  sampleChr = FALSE
  chrs = seqlevels(bpRegionsList[["Random"]])
  ## a chromosome can be given (for translocations and insertions) or sampled randomly when set to NA (other SVs)
  if(is.na(chr)){
    sampleChr = TRUE
  }
    
  ## sample mechanism
  idx = weightsMechanisms[,1] > 0  # col-index is 1, because main function already passed the right subset of weights
  mechanism = sample(rownames(weightsMechanisms)[idx], 1, prob=weightsMechanisms[idx,1])
  weightsRepeats = weightsRepeats[, mechanism, drop=FALSE]
  ## sample region type (kind of repeat or any other region)
  idx = weightsRepeats > 0
  regionType = sample(rownames(weightsRepeats)[idx], 1, prob=weightsRepeats[idx, 1])
  
  ## while loop in case there no suitable bpregion can be found within any chromosome
  while(length(bpRegions) == 0 & length(chrs) > 0){

    ## sample chromosome if it is missing; larger chromosomes have a higher probability
    if(sampleChr == TRUE){
      probs = sapply(chrs, function(x){return(sum(width(bpRegionsList[["Random"]][seqnames(bpRegionsList[["Random"]]) == x])))})
      chr = sample(chrs, 1, prob=probs)
      chrs = chrs[chrs != chr]
    }
    
    bpRegions = bpRegionsList[[regionType]]
    bpRegions = bpRegions[seqnames(bpRegions) == chr]

    ## for NAHR, set the breakpoints within a repeat (plus some tolerance towards the repeat margins)
    ## make sure the regions do not overlap with valid "normal" regions
    if(mechanism == "NAHR" & regionType != "Random"){
      tol = 50
      bpRegions = bpRegions[width(bpRegions) >= size+(tol*2)]
      if(length(bpRegions) > 0){
        bpRegions = bpRegions - tol # tolerance to make sure, that breakpoints lie not too close to the repeat margins
        bpRegions = bpRegions[queryHits(findOverlaps(bpRegions, bpRegionsList[["Random"]], type="within"))]
      }
    }
    ## for NHR, TEI and VNTR set one of both breakpoints within the repeats
    ## but make sure they are large enough by extending the regions at the start or end (randomly)
    if((length(bpRegions) == 0 | mechanism %in% c("NHR", "TEI", "VNTR")) & regionType != "Random"){
      tol = 50
      tooSmall = width(bpRegions) < (size + tol)
      extendStart = sample(c(TRUE,FALSE), length(bpRegions), replace=TRUE)
      diff = size - width(bpRegions) + tol # plus some tolerance (50bp)
      start(bpRegions[tooSmall&extendStart]) = start(bpRegions[tooSmall&extendStart]) - diff[tooSmall&extendStart]
      end(bpRegions[tooSmall&!extendStart]) = end(bpRegions[tooSmall&!extendStart]) + diff[tooSmall&!extendStart]    
      bpRegions = bpRegions[queryHits(findOverlaps(bpRegions, bpRegionsList[["Random"]], type="within"))]
    }
    ## for any other "normal" mechanism and region just use the regions itself
    ## select only those which are large enough
    if(mechanism == "Other" | regionType == "Random"){
      bpRegions = bpRegions[width(bpRegions) >= size]
    }
    
    ## in case there is no repeat region large enough, try another chromosome or set the breakpoint randomly (of only one left)
    if(length(bpRegions) == 0){
      if(length(chrs) == 1){
        regionType = "Random"
      }
    }
  }

  ## abort, if there is no regions large enough for this SV, no matter which mechanism you choose
  if(length(bpRegions) == 0){
    stop("Sorry, not possible to find a spot on the genome, which is large enough for SV of size ", size, "bp")
  }
    
  ## randomly select a region (larger regions have higher probability)
  ## for NAHR and NHR, just use the flanking regions which already have the right size
  idx = sample(x=1:length(bpRegions), size=1, prob=width(bpRegions))
  if(mechanism %in% c("NAHR", "NHR")){
    start = start(bpRegions[idx])
    end = end(bpRegions[idx])
  }
  ## for any other regions (including TEIs and VNTRs), randomly select start and end within the region (regions are already large enough)
  if(mechanism %in% c("TEI", "VNTR", "Other")){
    start = sample(x=start(bpRegions[idx]):(end(bpRegions[idx])-size+1), size=1)
    end =  start + size -1
  }
#  return(data.frame(seqnames=unique(seqnames(bpRegions)), start=start, end=end, mechanism=mechanism, bpRegion=regionType))
  return(data.frame(seqnames=unique(seqnames(bpRegions)), start=start, end=end))
}


.drawPos_trans <- function(chr, genome, bpRegionsList, weightsMechanisms, weightsRepeats){  
 
  ## draw breakpoint position and set the start/end to the outmost coordinate (1 for 5', chromosome end for 3')
  pos = .drawPos(chr, bpRegionsList, weightsMechanisms, weightsRepeats, 1) ## here: size = 1 -> start = end
  
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

.simTranslocationPositions <- function(n, bpRegionsList, weightsMechanisms, weightsRepeats, genome, percBalancedTrans, sizes, verbose){
  
  posTrans_1 = .getDummyDataframe()
  posTrans_2 = .getDummyDataframe()
#  translocations = data.frame(Type="",ChrA=0,StartA=0,EndA=0,ChrB=0,StartB=0,EndB=0,Balanced=FALSE,BpSeqA_5prime="",BpSeqA_3prime="",BpSeqB_5prime="",BpSeqB_3prime="", stringsAsFactors=FALSE)[-1, ]
  invertedTrans = c()

  if(n > 0){
    
    if(verbose==TRUE) pb = txtProgressBar(min = 0, max = n, style = 3)
    
    for(i in 1:n){
      
      ## first translocation partner
      ## check if regions on two different chromosomes are available
#      if(length(unique(seqnames(bpRegions))) < 2){
#        stop("No two chromosomes available for translocations")
#      }
      ## sample chromosome; larger chromosomes have a higher probability
      chrs = seqlevels(bpRegionsList[["Random"]])
      probs = sapply(chrs, function(x){return(sum(width(bpRegionsList[["Random"]][seqnames(bpRegionsList[["Random"]]) == x])))})
      chr1 = as.character(sample(chrs, 1, prob=probs))
      
      pos1 = .drawPos_trans(chr1, genome, bpRegionsList, weightsMechanisms, weightsRepeats)
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
      bpRegionsList[["Random"]] = .subtractIntervals(bpRegionsList[["Random"]], pos)

      ## second translocation partner
      chrs = chrs[chrs != chr1] # make sure second translocated segments lies on different chromosome
      probs = sapply(chrs, function(x){return(sum(width(bpRegionsList[["Random"]][seqnames(bpRegionsList[["Random"]]) == x])))})
      chr2 = as.character(sample(chrs, 1, prob=probs))
      
      pos2 = .drawPos_trans(chr2, genome, bpRegionsList, weightsMechanisms, weightsRepeats)
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
      bpRegionsList[["Random"]] = .subtractIntervals(bpRegionsList[["Random"]], pos)
      
      ## always invert translocated segments between different ends (5'<->3' or 3'<->5')
      if((pos1$start == 1 & pos2$start != 1) | (pos1$start != 1 & pos2$start == 1)){
        invertedTrans = c(invertedTrans, i)
      }
      
      if(verbose==TRUE) setTxtProgressBar(pb, i)

    }
    
    if(verbose==TRUE) close(pb)
    
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

.simInsertionPositions <- function(n, bpRegionsList, weightsMechanisms, weightsRepeats, genome, sizes, percCopiedIns, verbose){
  
  posIns_1 = .getDummyDataframe()
  posIns_2 = .getDummyDataframe()

  if(n > 0){
    if(verbose==TRUE) pb = txtProgressBar(min = 0, max = n, style = 3)
    for(i in 1:n){
     
      ## first translocation partner

      ## sample chromosome; larger chromosomes have a higher probability
      chrs = seqlevels(bpRegionsList[["Random"]])
      probs = sapply(chrs, function(x){return(sum(width(bpRegionsList[["Random"]][seqnames(bpRegionsList[["Random"]]) == x])))})
      chr1 = as.character(sample(chrs, 1, prob=probs))
      
      pos1 = .drawPos(chr1, bpRegionsList, weightsMechanisms, weightsRepeats, sizes[i])
      posIns_1 = rbind(posIns_1, pos1)
      
      ## make sure new SVs do not overlap with this one
      bpRegionsList[["Random"]] = .subtractIntervals(bpRegionsList[["Random"]], GRanges(IRanges(pos1$start, pos1$end), seqnames=pos1$seqnames))
      
      ## second translocation partner
      chrs = chrs[chrs != chr1] # make sure second translocated segments lies on different chromosome
      probs = sapply(chrs, function(x){return(sum(width(bpRegionsList[["Random"]][seqnames(bpRegionsList[["Random"]]) == x])))})
      chr2 = as.character(sample(chrs, 1, prob=probs))
      
      pos2 = .drawPos(chr2, bpRegionsList, weightsMechanisms, weightsRepeats, sizes[i])
      posIns_2 = rbind(posIns_2, pos2)
      
      ## make sure new SVs do not overlap with this one
      bpRegionsList[["Random"]] = .subtractIntervals(bpRegionsList[["Random"]], GRanges(IRanges(pos2$start, pos2$end), seqnames=pos2$seqnames))
      
      if(verbose==TRUE) setTxtProgressBar(pb, i)
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

  if(verbose==TRUE) close(pb)

  return(list(insertions, posIns_1, posIns_2))
}  


## 3.  deletions, inversions and tandem duplications #################################
.simPositions <- function(n, bpRegionsList, weightsMechanisms, weightsRepeats, sizes, type, verbose){
  
  pos = .getDummyDataframe()

  if(n > 0){
    if(verbose==TRUE) pb = txtProgressBar(min = 0, max = n, style = 3)
    for(i in 1:n){      
 
      p = .drawPos(NA, bpRegionsList, weightsMechanisms, weightsRepeats, sizes[i])
      pos = rbind(pos, p)    
      bpRegionsList[["Random"]] = .subtractIntervals(bpRegionsList[["Random"]], GRanges(IRanges(p$start, p$end), seqnames=p$seqnames))
      if(verbose==TRUE) setTxtProgressBar(pb, i)
    }
    
    svs = cbind("", pos, 0)
    colnames(svs) = c("Name", "Chr", "Start","End", "Size")
    svs$Size = sizes
    svs$Name = paste(type, 1:n, sep="")

  }
  
  if(verbose==TRUE) close(pb)

  return(list(svs, pos))
}
