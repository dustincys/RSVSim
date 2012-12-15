setGeneric("simulateSV",
           function(output=".", genome, chrs, dels=0, ins=0, invs=0, dups=0, trans=0, size, sizeDels=10, sizeIns=10, sizeInvs=10, sizeDups=10, regionsDels, regionsIns, regionsInvs, regionsDups, regionsTrans, maxDups=10, percCopiedIns=0, percBalancedTrans=1, transInsert=0, bpSeqSize=100, random=TRUE, seed) {
             standardGeneric("simulateSV")
           })

setGeneric("estimateSVSizes",
           function(n, svSizes, minSize=NA, maxSize=NA, default="deletions", hist=TRUE){
             standardGeneric("estimateSVSizes")
           })

setGeneric("compareSV",
           function(querySVs, simSVs, tol=200){
             standardGeneric("compareSV")
           })
