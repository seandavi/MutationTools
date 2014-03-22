#' Read a MAF file
#'
#' Reads in a MAF format file
#'
#' @param MAF The MAF file name
#' @param genome A string giving the genome name (eg., 'hg19')
#'
#' @return A \code\link{VRanges}} object
#' 
#' @import VariantAnnotation
#' @importClassesFrom VariantAnnotation VRanges
#' @importFrom VariantAnnotation VRanges
#' 
readMAF <- function(MAF,genome,...) {
    maf = read.delim(MAF,comment.char='#',...)
    
    rowRanges = VRanges(seqnames=Rle(maf$Chromosome),
        ranges=IRanges(start=maf$Start_position,end=maf$End_position),
        ref=as.character(maf$Reference_Allele),
        alt=as.character(maf$Tumor_Seq_Allele1),
        sampleNames=maf$Tumor_Sample_Barcode)
    return(rowRanges)
}
