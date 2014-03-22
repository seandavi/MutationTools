#' Gather variant context around single nucleotide variants
#'
#' The variant context is defined by gathering the nucleotides
#' upstream and downstream of a variant.  This if often the first
#' step in building a mutation spectrum analysis.
#'
#' 
#' @param variants A \code{\link{VCF-class}} object
#' @param fasta A \code{\link{FaFile}} object that describes the
#' fasta file location of the reference sequence
#' @param extendBy The number of bases on either side of the
#' altered base to capture.  The altered base will be placed in the
#' center of the context.  For example, the default value of 1
#' will result in contexts of length 3.  Specifying 2 will result
#' in contexts of length 5, etc.
#'
#' @export
#'
#' @return A \code{\link{VariantContext}} object.
#'
#' @importClassesFrom VariantAnnotation VCF
#' @importClassesFrom Rsamtools FaFile
#'
#' @examples
#' library(VariantAnnotation)
#' fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' param <- ScanVcfParam(info=c("DP"))
#' vcf = readVcf(fl,"hg19",param=param)
#' seqlevels(vcf) = paste0('chr',seqlevels(vcf))
#' variantContext(vcf,fasta=FaFile('/data/CCRBioinfo/public/GATK/bundle/2.3/hg19/ucsc.hg19.fasta'))
setGeneric('variantContext',signature=c('VARIANTS','FASTA'),
           function(VARIANTS,FASTA,...) standardGeneric('variantContext'))



#' @describeIn variantContext
setMethod('variantContext',c('VCF','FaFile'),
          function(VARIANTS,FASTA,...) {
              #clear out info fields before converting to VRanges
              info(VARIANTS)=DataFrame(rep(NA,nrow(VARIANTS)))
              vrange = as(VARIANTS,'VRanges')
              return(.variantContext(vrange,FASTA,...))
          })
              
#' @describeIn variantContext
#'
setMethod('variantContext',c('VRanges','FaFile'),
          function(VARIANTS,FASTA,...) {
              return(.variantContext(VARIANTS,FASTA,...))
          })


.variantContext <- function(variants,fasta,extendBy=1) {
    simpleIdx = which((nchar(ref(variants))==1) & (nchar(alt(variants))==1))
    rowRanges = GRanges(seqnames=seqnames(variants),ranges=ranges(variants))
    context = scanFa(fasta,resize(rowRanges[simpleIdx],extendBy*2+1,'center'))
    tmpc = as.character(context)
    names(tmpc)=names(context)
    substr(tmpc,2,2) = alt(variants)[simpleIdx]
    context2 = DNAStringSet(as.character(tmpc))
    GTidx = which(as.character(ref(variants[simpleIdx])) %in% c('G','T'))
    context[GTidx]=reverseComplement(context[GTidx])
    context2[GTidx]=reverseComplement(context2[GTidx])
    rowRanges = GenomicRanges:::newGRanges('VariantContext',seqnames=seqnames(rowRanges[simpleIdx]),
        ranges=ranges(rowRanges[simpleIdx]),strand=Rle('*',length(context)),
        mcols=DataFrame(refContext = context,
            altContext = context2,
            sampleName = factor(as.character(sampleNames(variants)[simpleIdx]))),
        seqlengths=seqlengths(variants),seqinfo=seqinfo(variants))
    return(rowRanges)
}

doMatrix = function(vc,sampleName=NULL) {
    if(!is.null(sampleName)) {
        vc=vc[vc$sampleName %in% sampleName]
    }
    allDNA=c('A','C','G','T')
    tmp = expand.grid(allDNA,allDNA,allDNA,c("A","C"),stringsAsFactors=FALSE)
    tmp = tmp[tmp[,2]!=tmp[,4],]
    tocolumn = apply(tmp[,1:3],1,paste,collapse="")
    fromcolumn = apply(tmp[,c(1,4,3)],1,paste,collapse="")
    df = data.frame(from = tmp[,4],to=tmp[,2],context=apply(tmp[,c(1,3)],1,paste,collapse="_"),
        fromcolumn,tocolumn)
    rownames(df)=NULL
    t1 = table(as.character(vc$refContext),as.character(vc$altContext))
    df$count = apply(df,1,function(x) {
        if((x[4] %in% rownames(t1)) & (x[5] %in% colnames(t1))) {
            return(t1[x[4],x[5]])
        } else {
            return(0)
        }
    })
    return(df)
}

plotVariantContext = function(vc) {
    require(ggplot2)
    ggplot(data=doMatrix(vc),
           aes(x=context,y=count,fill=context)) +
               geom_bar(stat='identity') +
               facet_grid(to ~ from) +
               theme(axis.text.x=element_text(angle=45,hjust=1))
}
