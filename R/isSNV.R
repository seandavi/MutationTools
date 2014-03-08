#' Determine if variants in a VCF object are SNPs
#'
#' Returns TRUE for variants that are SNVs and FALSE otherwise.
#' For variants with multiple ALT alleles, only the FIRST is used.
#'
#' @param VCF an object inheriting from the VCF class
#' @return logical vector with the same length as \code{vcf}
#' @keywords character
#' @seealso \code{\link{VCF-class}}
#' @export
#' @examples
#' library(VariantAnnotation)
#' fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' param <- ScanVcfParam(fixed="ALT", geno=c("GT", "GL"), info=c("LDAF"))
#' vcf = readVcf(fl,"hg19",param=param)
#' table(isSNV(vcf))
isSNV <- function(vcf) {
    if(!inherits(vcf,'VCF')) {
        stop('The isSimpleVariant function needs an object that inherits from the VCF class')
    }
    altall = unlist(alt(vcf))[start(PartitioningByEnd(alt(vcf)))]
    return((elementLengths(ref(vcf))==1) & (nchar(altall)==1))
}
