#' Calculate GC content in window around a set of variants
#'
#' @param vcf A \code{\link{VCF-class}} object
#' @param fasta A \code{\link{FaFile}} object
#' @param windowSize An integer representing the window size around each variant
#' to calculate GC content
#' @return A numeric vector giving the GC content for each variant window in the
#' input
#' @export
#' @examples \dontrun{
#' require(VariantAnnotation)
#' require(AnnotationHub)
#' fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' param <- ScanVcfParam(fixed="ALT", geno=c("GT", "GL"), info=c("LDAF"))
#' vcf = readVcf(fl,"hg19",param=param)
#' ah = AnnotationHub()
#' hg19fa = ah$ensembl.release.72.fasta.homo_sapiens.dna.Homo_sapiens.GRCh37.72.dna.toplevel.fa.rz
#' varGC = variantGCContent(vcf,hg19fa)
#' plot(density(varGC))
#' }
variantGCContent = function(vcf,fasta,windowSize=100) {
    if(!inherits(vcf,'VCF')) {
        stop('This function requires a VCF-class object')
    }
    if(!inherits(fasta,'FaFile')) {
        stop('This function requires an FaFile object')
    }
    oligoFreq = oligonucleotideFrequency(scanFa(fasta,resize(rowData(vcf),windowSize,fix='center')),1)
    return((oligoFreq[,'C']+oligoFreq[,'G'])/rowSums(oligoFreq))
}
