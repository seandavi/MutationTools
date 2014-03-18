#'
#'
#'
#' @importFrom Biostrings complement
kataegisPlot <- function(vcf) {
    vcfsimple = vcf[isSNV(vcf)]
    d = distanceToPreviousVariant(vcfsimple)
    vardf = GRanges(seqnames=seqnames(vcfsimple),ranges=ranges(vcfsimple),
        ref = DNAStringSet(ref(vcfsimple)),
        alt = DNAStringSet(unlist(alt(vcfsimple))[start(PartitioningByEnd(alt(vcfsimple)))]))
    mcols(vardf)$alt[mcols(vardf)$ref %in% c('G','T')] = complement(mcols(vardf)$alt[mcols(vardf)$ref %in% c('G','T')])
    mcols(vardf)$ref[mcols(vardf)$ref %in% c('G','T')] = complement(mcols(vardf)$ref[mcols(vardf)$ref %in% c('G','T')])
    x = data.frame(mutation=paste(mcols(vardf)$ref,mcols(vardf)$alt,sep=" > "),
        distance=d)
    return(x)
}
