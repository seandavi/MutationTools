#' Read in the results of running GASV
#'
#' Read in the results of running GASV
#'
#' @param fname The file name of the GASV "clusters" file
#'
#' @importFrom GenomicRanges GRanges,DataFrame,IRanges
#' 
#' @return A data frame
#'
readGASV = function(fname) {
    splitRanges = function(col) {
        # returns IRanges from a comma-separated column, like in GASV
        tmp = do.call(rbind,strsplit(as.character(col),','))
        return(IRanges(start=as.numeric(tmp[,1]),
                       end  =as.numeric(tmp[,2])))
    }
    d = read.delim(fname)
    seqnms = paste0('chr',as.character(d[,2]))
    seqnms[seqnms=='chr23']='chrX'
    seqnms[seqnms=='chr24']='chrY'
    leftbp = GRanges(seqnames = seqnms,
        ranges = splitRanges(d[,3])
    )
    rightbp = GRanges(seqnames = seqnms,
        ranges = splitRanges(d[,5])
    )
    df = d[,6:8]
    return(new('StrucVar',leftBP=leftbp,rightBP=rightbp,mcols=DataFrame(df)))
}
