#' wrapper for WGBS settings appropriate to dmrseq
#' 
#' @param bs              a bsseq object
#' @param testCovariate   the pData column to test on 
#' @param bpSpan          span of smoother AND 2x max gap in DMR CpGs (1000)
#' @param ...             other arguments to pass along to dmrseq
#'
#' @return                a GRanges object (same as from dmrseq)
#' 
#' @import dmrseq
#'
#' @export
WGBSeq <- function(bs, testCovariate, bpSpan=1000, ...) { 

  dmrseq(filterLoci(bs, testCovariate), testCovariate=testCovariate, 
         bpSpan=bpSpan, maxGap=(bpSpan/2), maxGapSmooth=bpSpan*2, ...)

}

