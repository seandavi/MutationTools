#'
#'
setClass("StrucVar",
         representation=list(leftBP='GRanges',
             rightBP='GRanges',
             mcols='DataFrame'))

setGeneric('leftBP',function(object) {standardGeneric('leftBP')})
setMethod('leftBP',c('StrucVar'),definition = function(object) {return(object@leftBP)})

setGeneric('rightBP',function(object) {standardGeneric('rightBP')})
setMethod('rightBP',c('StrucVar'),definition = function(object) {return(object@rightBP)})

