setMethod("simulateSV",
          signature(),

          function(output, genome, chrs, dels, ins, invs, dups, trans, size, sizeDels, sizeIns, sizeInvs, sizeDups, regionsDels, regionsIns, regionsInvs, regionsDups, regionsTrans, maxDups, percCopiedIns, percBalancedTrans, transInsert, bpSeqSize, random, seed){
            

  options(stringsAsFactors=FALSE, scipen=10)
  
  if(!missing(seed)){ 
    set.seed(seed)
  }

  ## default for random is TRUE
  ## is only one value given, set it for all SVs where regions were specified (only place, where a FALSE makes sense)
  if(length(random) == 1){
    r = rep(TRUE, 5)
    r[!c(missing(regionsDels), missing(regionsIns), missing(regionsInvs), missing(regionsDups), missing(regionsTrans))] = random
    random = r
  }
  
  .validateInput(output, genome, chrs, dels, ins, invs, dups, trans, sizeDels, sizeIns, sizeInvs, sizeDups, regionsDels, regionsIns, regionsInvs, regionsDups, regionsTrans, maxDups, percCopiedIns, percBalancedTrans, transInsert, bpSeqSize, random, seed)

  ## see if genome is given
  ## if not, read the hg19 by default
  if(missing(genome)){
    message("Loading human genome (hg19)")
    if(missing(chrs)){
      genome = .getHG19(NA)
    }else{
      genome = .getHG19(chrs)
    }
    gaps = genome[[2]]
    genome = genome[[1]]
    chrs = names(genome)
  ## if genome is given, there are two possibilities:
  ## 1. genome is a filename pointing to a FASTA-file
  ## 2. genome is a named DNAStringSet
  }else{
    if(class(genome) == "character"){
      genome = readDNAStringSet(filepath=genome)  ## read genome from FASTA-file
    }
    if(missing(chrs)){
      chrs = names(genome)
    }
    ## compute gaps within the given genome (N-sequences)
    gaps = GRanges()
    for(chr in chrs){
      Ns = whichAsIRanges(as.integer(genome[[chr]]) == as.integer(DNAString("N")))
      if(length(Ns) > 0){
        gaps = c(gaps, GRanges(Ns, seqnames=chr))
      }
    }
  }
  gaps = as.data.frame(gaps)[, 1:3]
  genomeCoords = GRanges(IRanges(start=1, end=width(genome)), seqnames=names(genome))

  ## put SVs anywhere in the genome if no regions were specified
  ## if regions were given, take care, that the regions lie within the available chromosomes
  if(missing(regionsDels)){
    regionsDels = genomeCoords
  }else{
    regionsDels = regionsDels[seqnames(regionsDels) %in% intersect(chrs, as.character(seqnames((regionsDels))))]
    if(length(regionsDels) == 0){
      stop("No regions on given chromosomes")
    }
  }
  if(missing(regionsIns)){
    regionsIns = genomeCoords
  }else{
    if(random[2] == TRUE){
      regionsIns = regionsIns[seqnames(regionsIns) %in% intersect(chrs, as.character(seqnames((regionsIns))))]
    }else{
      regionsIns = regionsIns[seqnames(regionsIns) %in% intersect(chrs, as.character(seqnames((regionsIns)))) & regionsIns$chrB %in% chrs]      
    }
    if(length(regionsIns) == 0){
      stop("No regions on given chromosomes")
    }
  }
  if(missing(regionsInvs)){
    regionsInvs = genomeCoords
  }else{
    regionsInvs = regionsInvs[seqnames(regionsInvs) %in% intersect(chrs, as.character(seqnames((regionsInvs))))]
    if(length(regionsInvs) == 0){
      stop("No regions on given chromosomes")
    }
  }
  if(missing(regionsDups)){
    regionsDups = genomeCoords
  }else{
    regionsDups = regionsDups[seqnames(regionsDups) %in% intersect(chrs, as.character(seqnames((regionsDups))))]
    if(length(regionsDups) == 0){
      stop("No regions on given chromosomes")
    }
  }
  if(missing(regionsTrans)){
    regionsTrans = genomeCoords
  }else{
    if(random[2] == TRUE){
      regionsTrans = regionsTrans[seqnames(regionsTrans) %in% intersect(chrs, as.character(seqnames((regionsTrans))))]
    }else{
      regionsTrans = regionsTrans[seqnames(regionsTrans) %in% intersect(chrs, as.character(seqnames((regionsTrans)))) & regionsTrans$chrB %in% chrs]
    }
    if(length(regionsTrans) == 0){
      stop("No regions on given chromosomes")
    }
  }

  ## restrict SV regions to given chromosomes

  ## inversion and deletion sizes
  if((length(sizeDels) > 1 & length(sizeDels) != dels) | (length(sizeIns) > 1 & length(sizeIns) != ins) | (length(sizeInvs) > 1 & length(sizeInvs) != invs) | (length(sizeDups) > 1 & length(sizeDups) != dups)){
    warning("Length of vectors with SV sizes vary from the number of SVs")
  }
  if(!missing(size)){
    sizeDels = sizeIns = sizeInvs = sizeDups = size
  }
  
  sizeDels = rep(sizeDels, dels)[1:dels]
  sizeIns = rep(sizeIns, ins)[1:ins]
  sizeInvs = rep(sizeInvs, invs)[1:invs]
  sizeDups = rep(sizeDups, dups)[1:dups]
  
  ## bpSeqSize will be used two times, in upstream and downstream direction; so divide it hear by two
  bpSeqSize = round(bpSeqSize / 2)
  

  ## Simulation ############################################################
  
  ## 1. simulate regions for translocations:
  type = "translocation"
  translocations = posTrans_1 = posTrans_2 = NULL
  if(random[5] == TRUE){
    ## 1.1 for chromosome terminals
    if(trans > 0){
      message("Calculating coordinates: ", trans, " translocations")
    
      subtrahend = gaps[, c("seqnames","start","end")]
      subtrahend = GRanges(IRanges(subtrahend$start, subtrahend$end), seqnames=subtrahend$seqnames)
      regionsTrans  = .subtractIntervals(regionsTrans, subtrahend)
      t = .simTranslocationPositions(trans, regionsTrans, genome, percBalancedTrans, c())
      translocations = t[[1]]
      posTrans_1 = t[[2]]
      posTrans_2 = t[[3]]
    }   
  }else{
    message("Using given coordinates: Translocations")    
    trans = length(regionsTrans)
    posTrans_1 = as.data.frame(regionsTrans)
    posTrans_1$names = rownames(posTrans_1)
    posTrans_2 = as.data.frame(regionsTrans)[, c("chrB","startB","endB")]
    colnames(posTrans_2) = c("seqnames","start","end")
    
    ## Add balanced information (default:TRUE) if not given
    if(!("balanced" %in% colnames(posTrans_1)) | !("balanced" %in% colnames(posTrans_2))){
      posTrans_1$balanced = posTrans_2$balanced = TRUE
    }
    ## Determine wich translocations need to be inverted (segments on different ends of chromosomes)
    idx = which((posTrans_1$start == 1 & posTrans_2$start != 1) | (posTrans_1$start != 1 & posTrans_2$start == 1))
    posTrans_1$inverted = posTrans_2$inverted = FALSE
    posTrans_1$inverted[idx] = posTrans_2$inverted[idx] = TRUE

    size1 = posTrans_1$end-posTrans_1$start+1
    size2 = posTrans_2$end-posTrans_2$start+1
    translocations = cbind(posTrans_1$names, posTrans_1[, c("seqnames","start","end")], size1, posTrans_2[, c("seqnames","start","end")], size2, posTrans_2$balanced)
    colnames(translocations) = c("Name", "ChrA", "StartA", "EndA", "SizeA", "ChrB", "StartB", "EndB", "SizeB","Balanced")
  }
  ## 2. simulate insertions
  type = "insertion"
  insertions = posIns_1 = posIns_2 = NULL
  if(random[2] == TRUE){
    if(ins > 0){
      message("Calculating coordinates: ", ins, " insertions")  
      
      subtrahend = rbind(gaps[, c("seqnames","start","end")], posTrans_1[, c("seqnames","start","end")], posTrans_2[, c("seqnames","start","end")])
      subtrahend = GRanges(IRanges(subtrahend$start, subtrahend$end), seqnames=subtrahend$seqnames)
      regionsIns  = .subtractIntervals(regionsIns, subtrahend)
      i = .simInsertionPositions(ins, regionsIns, genome, sizeIns, percCopiedIns)
      insertions = i[[1]]
      posIns_1 = i[[2]]
      posIns_2 = i[[3]]
    }
  }else{
    message("Using given coordinates: Insertions")    
    ins = length(regionsIns)
    posIns_1 = as.data.frame(regionsIns)
    size = posIns_1$end-posIns_1$start+1
    posIns_1$names = rownames(posIns_1)
    posIns_2 = as.data.frame(regionsIns)[, c("chrB","startB")]
    posIns_2$endB = posIns_2$startB + size - 1
    colnames(posIns_2) = c("seqnames","start","end")
    
    ## Add duplicate information (default:FALSE) if not given
    if(!("copied" %in% colnames(posIns_1)) | !("copied" %in% colnames(posIns_2))){
      posIns_1$copied = posIns_2$copied = FALSE
    }
 
    insertions = cbind(posIns_1$names, posIns_1[rownames(posIns_1), c("seqnames","start","end")], posIns_2[rownames(posIns_1), c("seqnames","start","end")], size, posIns_2$copied)
    colnames(insertions) = c("Name", "ChrA", "StartA", "EndA", "ChrB", "StartB", "EndB", "Size", "Copied")
  }

  ## 3. simulate deletions
  type = "deletion"
  deletions = posDel = NULL
  if(random[1] == TRUE){
    if(dels > 0){
      message("Calculating coordinates: ", dels, " deletions")
      
      subtrahend = rbind(gaps[, c("seqnames","start","end")], posIns_1[, c("seqnames","start","end")], posIns_2[, c("seqnames","start","end")], posTrans_1[, c("seqnames","start","end")], posTrans_2[, c("seqnames","start","end")])
      subtrahend = GRanges(IRanges(subtrahend$start, subtrahend$end), seqnames=subtrahend$seqnames)
      regionsDels  = .subtractIntervals(regionsDels, subtrahend)
      d = .simPositions(dels, regionsDels, sizeDels, "deletion")
      deletions = d[[1]]
      posDel = d[[2]]
    }
  }else{
    message("Using given coordinates: Deletions")
    dels = length(regionsDels)
    posDel = as.data.frame(regionsDels)
    posDel$names = rownames(posDel)
    size=posDel$end-posDel$start+1
    deletions = cbind(posDel$names, posDel[, c("seqnames","start","end")], size)
    colnames(deletions) = c("Name", "Chr", "Start","End","Size")
  }

  ## 4. simulate inversions
  type = "inversion"
  inversions = posInv = NULL
  if(random[3] == TRUE){
    if(invs > 0){
      message("Calculating coordinates: ", invs, " inversions")
      
      subtrahend = rbind(gaps[, c("seqnames","start","end")], posDel[, c("seqnames","start","end")], posIns_1[, c("seqnames","start","end")], posIns_2[, c("seqnames","start","end")], posTrans_1[, c("seqnames","start","end")], posTrans_2[, c("seqnames","start","end")])
      subtrahend = GRanges(IRanges(subtrahend$start, subtrahend$end), seqnames=subtrahend$seqnames)
      regionsInvs  = .subtractIntervals(regionsInvs, subtrahend)
      i = .simPositions(invs, regionsInvs, sizeInvs, type)
      inversions = i[[1]]
      posInv = i[[2]]
    }
  }else{
    message("Using given coordinates: Inversions")
    invs = length(regionsInvs)
    posInv = as.data.frame(regionsInvs)
    posInv$names = rownames(posInv)
    size=posInv$end-posInv$start+1
    inversions = cbind(posInv$names, posInv[, c("seqnames","start","end")], size)
    colnames(inversions) = c("Name","Chr", "Start","End", "Size")
  }
  
  ## 5. simulate tandem duplications
  type = "tandemDuplication"
  tandemDups = posDup = NULL
  if(random[4] == TRUE){
    if(dups > 0){
      message("Calculating coordinates: ", dups, " tandem duplications")
      
      subtrahend = rbind(gaps[, c("seqnames","start","end")], posDel[, c("seqnames","start","end")], posIns_1[, c("seqnames","start","end")], posIns_2[, c("seqnames","start","end")], posInv[, c("seqnames","start","end")], posTrans_1[, c("seqnames","start","end")], posTrans_2[, c("seqnames","start","end")])
      subtrahend = GRanges(IRanges(subtrahend$start, subtrahend$end), seqnames=subtrahend$seqnames)
      regionsDups  = .subtractIntervals(regionsDups, subtrahend)
      td = .simPositions(dups, regionsDups, sizeDups, type)
      tandemDups = td[[1]]
      posDup = td[[2]]
    }
  }else{
    message("Using given coordinates: Tandem Duplications")
    dups = length(regionsDups)
    posDup = as.data.frame(regionsDups)
    posDup$names = rownames(posDup)
    size=posDup$end-posDup$start+1
    tandemDups = cbind(posDup$names, posDup[, c("seqnames","start","end")], size)
    colnames(tandemDups) = c("Name","Chr", "Start","End", "Size")
  }
  
  ## Execution: implement regions into genome sequence ############################################################
    
  ## 1. inversions
  if(invs > 0){
    message("Rearranging genome: Inversions")
    for(i in 1:invs){
      chr = as.character(posInv$seqnames[i])
      genome[[chr]] = .execInversion(genome[[chr]], posInv, bpSeqSize, i)
    }
  }
  
  ## 2. translocations
  if(trans > 0){
    message("Rearranging genome: Translocations")
    for(i in 1:trans){
      chr1 = as.character(posTrans_1$seqnames[i])
      chr2 = as.character(posTrans_2$seqnames[i])
      rearrangement = .execTranslocation(genome[[chr1]], genome[[chr2]], chr1, chr2, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, i, bpSeqSize, transInsert)
      genome[[chr1]] = rearrangement[[1]]
      genome[[chr2]] = rearrangement[[2]]
      posDel = rearrangement[[3]]
      posIns_1 = rearrangement[[4]]
      posIns_2 = rearrangement[[5]]
      posInv = rearrangement[[6]]
      posDup = rearrangement[[7]]
      posTrans_1 = rearrangement[[8]]
      posTrans_2 = rearrangement[[9]]
    }

  }
  
  ## 3. insertions
  if(ins > 0){
    message("Rearranging genome: Insertions")
    for(i in 1:ins){
      chr1 = as.character(posIns_1$seqnames[i])
      chr2 = as.character(posIns_2$seqnames[i])      
      rearrangement = .execInsertion(genome[[chr1]], genome[[chr2]], chr1, chr2, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, i, bpSeqSize)
      genome[[chr1]] = rearrangement[[1]]
      genome[[chr2]] = rearrangement[[2]]
      posDel = rearrangement[[3]]
      posIns_1 = rearrangement[[4]]
      posIns_2 = rearrangement[[5]]
      posInv = rearrangement[[6]]
      posDup = rearrangement[[7]]
      posTrans_1 = rearrangement[[8]]
      posTrans_2 = rearrangement[[9]]
    }

  }
  
  ## 4. tandem duplications
  if(dups > 0){
    message("Rearranging genome: Tandem duplications")
    ## add a column for the number of duplications
    tandemDups = cbind(tandemDups[, c("Name", "Chr", "Start","End", "Size")], 0, "")
    names(tandemDups) = c("Name", "Chr", "Start","End", "Size", "Duplications", "BpSeq")
    for(i in 1:dups){
      times = sample(2:maxDups,1)  ## how many times the sequence is duplicated
      chr = as.character(posDup$seqnames[i])
      rearrangement = .execTandemDuplication(genome[[chr]], chr, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, bpSeqSize, times, i)
      genome[[chr]] = rearrangement[[1]]
      posDel = rearrangement[[2]]
      posIns_1 = rearrangement[[3]]
      posIns_2 = rearrangement[[4]]
      posInv = rearrangement[[5]]
      posDup = rearrangement[[6]]
      posTrans_1 = rearrangement[[7]]
      posTrans_2 = rearrangement[[8]]      
      tandemDups$BpSeq[i] = rearrangement[[9]]
      tandemDups$Duplications[i] = times
    }
  }

  ## 5. deletions
  if(dels > 0){
    message("Rearranging genome: Deletions")
    for(i in 1:dels){
      chr = as.character(posDel$seqnames[i])
      rearrangement = .execDeletion(genome[[chr]], chr, posDel, posIns_1, posIns_2, posInv, posDup, posTrans_1, posTrans_2, bpSeqSize, i)
      genome[[chr]] = rearrangement[[1]]
      posDel = rearrangement[[2]]
      posIns_1 = rearrangement[[3]]
      posIns_2 = rearrangement[[4]]
      posInv = rearrangement[[5]]
      posDup = rearrangement[[6]]
      posTrans_1 = rearrangement[[7]]
      posTrans_2 = rearrangement[[8]]
    }
  }

  ## Retreive breakpoint sequences for translocations, deletions and inversions (for duplications, it works during execution)
  if(invs > 0){
    inversions$BpSeq_5prime = inversions$BpSeq_3prime = ""
    for(i in 1:invs){
      bpSeq = .getBpSeq(genome, posInv, bpSeqSize, i)
      inversions$BpSeq_5prime[i] = bpSeq[1]
      inversions$BpSeq_3prime[i] = bpSeq[2]
    }
  }

  if(dels > 0){
    deletions$BpSeq = ""
    for(i in 1:dels){
      pos = posDel
      pos$end = pos$start
      bpSeq = .getBpSeq(genome, pos, bpSeqSize, i)
      deletions$BpSeq[i] = bpSeq[1]
    }
  }
  
  if(trans > 0){
    translocations$BpSeqB = translocations$BpSeqA = ""
    for(i in 1:trans){
      if(posTrans_1$balanced[i] == TRUE){
        bpSeq = .getBpSeq(genome, posTrans_1, bpSeqSize, i)
        translocations$BpSeqA[i] = bpSeq[which(bpSeq != "")][1]
      }
      bpSeq = .getBpSeq(genome, posTrans_2, bpSeqSize, i)
      translocations$BpSeqB[i] = bpSeq[which(bpSeq != "")][1]
    }
  }
  if(ins > 0){
    insertions$BpSeqB_3prime = insertions$BpSeqB_5prime = insertions$BpSeqA = ""
    for(i in 1:ins){
      if(posIns_1$copied[i] == FALSE){
        bpSeq = .getBpSeq(genome, posIns_1, bpSeqSize, i)
        insertions$BpSeqA[i] = bpSeq[1]
      }
      bpSeq = .getBpSeq(genome, posIns_2, bpSeqSize, i)
      insertions$BpSeqB_5prime[i] = bpSeq[1]
      insertions$BpSeqB_3prime[i] = bpSeq[2]
    }
  }

  ## Output ############################################################
  if(!is.na(output)){
    if(output == "."){
      message("Writing output to current directory")
    }else{
      message("Writing output to ", output)
    }
    
    if(dels > 0){
      write.table(deletions, file.path(output, "deletions.csv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
    }
    if(ins > 0){
      write.table(insertions, file.path(output, "insertions.csv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
    }
    if(invs > 0){
      write.table(inversions, file.path(output, "inversions.csv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
    }
    if(dups > 0){
      write.table(tandemDups, file.path(output, "tandemDuplications.csv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
    }
    if(trans > 0){
      write.table(translocations, file.path(output, "translocations.csv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
    }
    
    writeXStringSet(genome, file.path(output, "genome_rearranged.fasta"), append=FALSE, format="fasta")
  }
  
  idx = c(dels>0, ins>0, invs>0, dups>0, trans>0)
  
#  names(svs) = c("deletions", "insertions", "inversions","tandemDuplications","translocations")[idx]
  metadata(genome) = list(deletions=deletions, insertions=insertions, inversions=inversions, tandemDuplications=tandemDups, translocations=translocations)[idx]
  
  return(genome)
  
})


.validateInput <- function(output, genome, chrs, dels, ins, invs, dups, trans, sizeDels, sizeIns, sizeInvs, sizeDups, regionsDels, regionsIns, regionsInvs, regionsDups, regionsTrans, maxDups, percCopiedIns, percBalancedTrans, transInsert, bpSeqSize, random, seed){

  ## output
  if(!is.na(output)){
    if(!file.exists(output)){
      stop("Output directory does not exist")
    }
  }
  ## genome
  if(!missing(genome)){
    if(class(genome) == "DNAStringSet"){
      if(is.null(names(genome))){
        stop("Please enter chromosome names in your genome DNAStringSet")
      }
    }
  }
  ## chrs
  if(!missing(chrs) & !missing(genome)){
    if(!all(chrs %in% names(genome))){
      stop("Invalid argument: Specified chromosomes and chromosome names in the genome do not match")
    }
  }

  ## Number and size of SVs
  if(any(c(dels, ins, invs, dups, trans, sizeDels, sizeIns, sizeInvs, sizeDups) < 0)){
    stop("Invalid argument: Number of SVs and their size cannot be smaller than zero. Makes sense, doesn't it?")
  }
  
  ## percBalancedTrans and percduplicateIns
  if(percBalancedTrans < 0 | percBalancedTrans > 1 | percCopiedIns < 0 | percCopiedIns > 1){
    stop("Invalid argument: Percentages have to be given as value between 0 and 1.")
  }

  ## random and regions
  if(length(random) != 1 & length(random) != 5){
    stop("Invalid argument: Give either one value for argument random for all SVs or five values for each SV in this order: deletions, insertions, inversions, tandem duplications, translocations")
  }
  if(random[1] == FALSE & missing(regionsDels)){
    stop("Missing argument: Please specifiy the regions for deletions")
  }
  if(random[2] == FALSE & missing(regionsIns)){
    stop("Missing argument: Please specifiy the regions for insertions")
  }
  if(random[3] == FALSE & missing(regionsInvs)){
    stop("Missing argument: Please specifiy the regions for inversions")
  }
  if(random[4] == FALSE & missing(regionsDups)){
    stop("Missing argument: Please specifiy the regions for tandem duplications")
  }
  if(random[5] == FALSE & missing(regionsTrans)){
    stop("Missing argument: Please specifiy the regions for translocations")
  }
  
  ## regions
  if(!missing(regionsIns) & random[2] == FALSE){
    if(!all(c("chrB","startB") %in% colnames(mcols(regionsIns)))){
      stop("Invalid argument: regionsIns is missing columns chrB and startB")
    }
  }
  if(!missing(regionsTrans) & random[5] == FALSE){
    if(!all(c("chrB","startB","endB") %in% colnames(mcols(regionsTrans)))){
      stop("Invalid argument: regionsTrans is missing columns chrB, startB and endB")
    }
  }

}


## Deprecated: Transposons (maybe useful for future release)
#
#    if(trans[3] > 0){
#      ## 1.3 for intra-chromosomal translocations
#    
#      message("Calculating coordinates: ", trans[3], " intra-chromosomal translocations")  
#      
#      intrachrom=TRUE
#      subtrahend = RangedData(rbind(gaps[, c("space","start","end")], posTrans_1[, c("space","start","end")], posTrans_2[, c("space","start","end")]))
#      regionsTrans  = .subtractIntervals(regionsTrans, subtrahend)
#      t = .simTranslocationPositions(trans[3], regionsTrans, genome, percInvertedTrans, percBalancedTrans, sizeTrans2, "intrachrom")
#      translocations = rbind(translocations, t[[1]])
#      posTrans_1 = rbind(posTrans_1, t[[2]])
#      posTrans_2 = rbind(posTrans_2, t[[3]])
#    }


## Maybe add gene annotation in future release (if it makes sense)
#    require(biomaRt)
#    ensembl=useMart("ensembl")
#    dataset="hsapiens_gene_ensembl"
#    ensembl=useDataset(dataset, mart=ensembl)
#    bmAttributes = c(
#      "hgnc_symbol",
#      "chromosome_name",  
#      "start_position",
#      "end_position"
#      )
#    bmFilter=c("chromosomal_region")
#    bmValues = as.list(paste(substr(d$Chr, 4, nchar(d$Chr)), d$Start, d$End, sep=":"))
#    genes = getBM(attributes=bmAttributes, filter=bmFilter, values=bmValues, mart=ensembl)
#    d$Genes = ""
#    for(i in 1:nrow(d)){
#      gene_overlap = subset(genes, (chromosome_name == substr(d$Chr[i], 4, nchar(d$Chr))) & (IRanges(start_position, end_position) %in% IRanges(d$Start[i], d$End[i])) & hgnc_symbol != "")
#      d$Genes[i] = paste(unique(gene_overlap$hgnc_symbol), collapse=",")
#    }
