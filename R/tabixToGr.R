
#' Read from tabix-indexed bed file to GRanges or tibble
#' The fourth column is DNA methylation level and the fifth is depth.
#' This is also the default output of Biscuit BED.
#' 
#' @param paths paths to the bed files
#' @param chrm  chromosome name
#' @param beg   start coordinate of CpG
#' @param end   end coordinate of CpG
#' @param as_gr whether to output GRanges
#' @param min_depth mask all depth < min_depth as NA. This is useful as depths
#' are not recorded in the output
#' @param sample_names sample names, just use paths if not specified
#' @param BPPARAM how to parallelize
#' @return a tibble or GRanges object with DNA methylation level
#'
#' @import GenomicRanges
#' @import BiocParallel
#' @import dplyr
#'
#' @export
tabixToGr <- function(
    paths, chrm, beg = 1, end = 9999999,
    as_gr = FALSE, min_depth = 0, sample_names = NULL,
    BPPARAM = SerialParam()) {

    input_range <- GRanges(chrm, IRanges(beg, end))
    df_list <- bplapply(paths, function(path) {
        df <- as.tibble(t(simplify2array(strsplit(scanTabix(
            path, param=input_range)[[1]], '\t'))), stringsAsFactors = FALSE)
        
        colnames(df) <- c('chrm','beg','end','beta','depth')
        df$beg <- as.numeric(df$beg)
        df$beg <- df$beg + 1
        df$end <- as.numeric(df$end)
        df$beta[df$beta == '.'] <- NA
        df$beta <- as.numeric(df$beta)
        df$depth <- as.integer(df$depth)
        df$depth[is.na(df$beta)] <- 0
        df
    })
    
    ## make sure the coordinates are the same
    same_coordinates <- sapply(seq_len(length(df_list)-1),
        function(i) identical(df_list[[i]][,1:3], df_list[[i+1]][,1:3]))
    stopifnot(all(same_coordinates))

    ## collate data
    data <- lapply(df_list, function(df) {
        ifelse(df$depth >= min_depth, df$beta, NA)
    })

    ## set sample names
    if (is.null(sample_names)) {
        sample_names <- paths
    }
    stopifnot(length(sample_names) == length(paths))
    names(data) <- sample_names

    ## output either as GR or a tibble
    if (as_gr) {
        res <- GRanges(
            seqnames = df_list[[1]]$chrm,
            ranges = IRanges(df_list[[1]]$beg, df_list[[1]]$end))
        mcols(res) <- as.data.frame(data)
    } else {
        res <- bind_cols(df_list[[1]][,1:3], as.tibble(data))
    }
    res
}
