#' Find the distance between each variant and the previous one
#'
#' Given a \code{\link{VCF-class}} object, return a vector of distances
#' between neighboring variants.  This is useful for plotting defining
#' Kataegis phenomenon in a set of variants.  Distance for the first
#' variant on each chromosome is NA.
#'
#' @param vcf A \code{\link{VCF-class}}
#' @export
#' @return A vector of nrow(vcf)-1 giving the distance between
#' each variant and the previous one.  
#' @examples
#' require(VariantAnnotation)
#' fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' param <- ScanVcfParam(fixed="ALT", geno=c("GT", "GL"), info=c("LDAF"))
#' vcf = readVcf(fl,"hg19",param=param)
#' d = distanceToPreviousVariant(vcf)
#' summary(d)
distanceToPreviousVariant <- function(vcf) {
    d = c(NA,start(rowData(vcf)[2:nrow(vcf)])-start(rowData(vcf))[1:nrow(vcf)-1])
    d[as.character(seqnames(rowData(vcf))[2:nrow(vcf)])!=as.character(seqnames(rowData(vcf))[1:nrow(vcf)-1])]=NA
    return(d)
}
