#' Determine if variants in a VCF object are SNPs
#'
#' Returns TRUE for variants that are SNVs and FALSE otherwise.
#' For variants with multiple ALT alleles, only the FIRST is used.
#'
#' @param variants an object inheriting from the \code\link{VCF}} or \code{\link{VRanges}} classes
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
isSNV <- function(variants) {

    if(inherits(variants,'VCF')) {
        refall = as.character(ref(variants))
        altall = as.character(unlist(alt(variants))[start(PartitioningByEnd(alt(variants)))])
        return((nchar(refall)==1) & (nchar(altall)==1))
    }

    if(inherits(variants,'VRanges')) {
        refall = as.character(ref(variants))
        altall = as.character(alt(variants))
        return((nchar(refall)==1) & (nchar(altall)==1))
    }

    stop('parameter variants must be a VRanges or VCF object')

}
