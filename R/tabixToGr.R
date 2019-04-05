
#' Read from tabix-indexed bed file to list objects
#' 
#' @param paths paths to the bed files
#' @param chrm  chromosome name
#' @param beg   start coordinate of CpG
#' @param end   end coordinate of CpG 2^29 is the max taken
#' @param sample_names sample names, just use paths if not specified
#' @param BPPARAM how to parallelize
#' @return a list object with DNA methylation level and depth
#'
#' @import BiocParallel
#' @import dplyr
#'
#' @export
tabixRetrieve <- function(
    paths, chrm, beg = 1, end = 2^28,
    min_depth = 0, sample_names = NULL,
    BPPARAM = SerialParam()) {

    input_range <- GRanges(chrm, IRanges(beg, end))
    df_list <- bplapply(paths, function(path) {
        df <- as_tibble(t(simplify2array(strsplit(scanTabix(
            path, param=input_range)[[1]], '\t'))), stringsAsFactors = FALSE)
        
        colnames(df) <- c('chrm','beg','end','beta','depth')
        df$beg <- as.numeric(df$beg)
        df$beg <- df$beg + 1
        df$end <- as.numeric(df$end)
        ## in case the tabix is mal-formed
        df$beta[df$beta == '.' | df$beta == 'NA'] <- NA
        df$beta <- as.numeric(df$beta)
        df$depth <- as.integer(df$depth)
        df$depth[is.na(df$beta)] <- 0
        df
    })
    
    ## make sure the coordinates are the same
    same_coordinates <- sapply(seq_len(length(df_list)-1),
        function(i) identical(df_list[[i]][,1:3], df_list[[i+1]][,1:3]))
    stopifnot(all(same_coordinates))

    ## set sample names
    if (is.null(sample_names)) {
        sample_names <- paths
    }
    stopifnot(length(sample_names) == length(paths))
    names(df_list) <- sample_names

    df_list
}

#' Read from tabix-indexed bed file into a tibble
#'
#' @param ... tabixRetrieve arguments
#' @param min_depth minimum depth as depth is discarded in GRanges
#' @return a tibble of DNA methylation level
#' @import dplyr
#' @seealso tabixRetrieve
#'
#' @export
tabixToTable <- function(..., min_depth = 0) {
    df_list <- tabixRetrieve(...)

    ## collate data (methylation levels)
    data <- lapply(df_list, function(df) {
        ifelse(df$depth >= min_depth, df$beta, NA)
    })
    bind_cols(df_list[[1]][,1:3], as.tibble(data))
}

#' Read from tabix-indexed bed file into a matrix
#' The row names are the CpG coordinates. The column names are
#' sample names.
#'
#' @param ... tabixRetrieve and tabixToTable arguments
#' @param min_depth minimum depth as depth is discarded in GRanges
#' @return a matrix of DNA methylation level
#' @seealso tabixRetrieve
#' @seealso tabixToTable
#'
#' @export
tabixToMatrix <- function(..., min_depth = 0) {

    df_list <- tabixRetrieve(...)

    ## collate data (methylation levels)
    mat <- do.call(cbind, lapply(df_list, function(df) {
        as.numeric(ifelse(df$depth >= min_depth, df$beta, NA))
    }))
    colnames(mat) <- names(df_list)
    rownames(mat) <- paste0(
        df_list[[1]]$chrm, ':',
        df_list[[1]]$beg, '-',
        df_list[[1]]$end)
    mat
}

#' Read from tabix-indexed bed file into a GRanges
#'
#' @param ... tabixRetrieve arguments
#' @param min_depth minimum depth as depth is discarded in GRanges
#' @return a GenomicRanges of DNA methylation level
#' @seealso tabixRetrieve
#' @import GenomicRanges
#'
#' @export
tabixToGr <- function(..., min_depth = 0) {
    df_list <- tabixRetrieve(...)

    ## collate data (methylation levels)
    data <- lapply(df_list, function(df) {
        ifelse(df$depth >= min_depth, df$beta, NA)
    })
    
    res <- GRanges(
        seqnames = df_list[[1]]$chrm,
        ranges = IRanges(df_list[[1]]$beg, df_list[[1]]$end))
    mcols(res) <- as.data.frame(data)
    res
}
