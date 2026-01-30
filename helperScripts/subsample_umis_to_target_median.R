library(Matrix)

subsample_umis_to_target_median <- function(count_matrix, target_median_umi) {
  if (!inherits(count_matrix, "dgCMatrix")) {
    stop("Input count_matrix must be a sparse matrix of class 'dgCMatrix'")
  }
  
  message("Calculating UMI counts per cell...")
  cell_umis <- Matrix::colSums(count_matrix)
  
  message("Computing downsampling fractions...")
  downsample_factors <- pmin(1, target_median_umi / cell_umis)
  
  message("Subsampling...")
  downsample_column <- function(counts_col, factor) {
    if (factor >= 1) {
      return(counts_col)
    }
    if (sum(counts_col) == 0) {
      return(counts_col)
    }
    unrolled <- rep(seq_along(counts_col), counts_col)
    subsampled <- sample(unrolled, size = round(length(unrolled) * factor))
    tab <- table(factor(subsampled, levels = seq_along(counts_col)))
    return(as.integer(tab))
  }
  
  # Initialize new sparse matrix
  downsampled_matrix <- Matrix(0, nrow = nrow(count_matrix), ncol = ncol(count_matrix), sparse = TRUE)
  
  for (i in seq_len(ncol(count_matrix))) {
    downsampled_matrix[, i] <- downsample_column(count_matrix[, i], downsample_factors[i])
    if (i %% 100 == 0) message("Processed ", i, " / ", ncol(count_matrix), " cells...")
  }
  
  message("Done. Final median UMI per cell: ", median(Matrix::colSums(downsampled_matrix)))
  return(downsampled_matrix)
}