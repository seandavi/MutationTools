
.valid.VariantContext = function(object) {
    return(TRUE)}

#' Contains variant context
#'
#' This is testing
#'
#' @importClassesFrom GenomicRanges GRanges
#' 
setClass("VariantContext",
#         representation(refContext = "DNAStringSet", # or XStringSet?
#                        altContext = "DNAStringSet",
#                        sampleName = "factor"),
         contains = "GRanges",
         validity = .valid.VariantContext)
