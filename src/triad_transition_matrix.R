# ================================================================
# File: triad_transition_matrix.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Purpose:
#   Compute a 4×4 transition count matrix of triad-type changes between
#   consecutive time points, aggregated over all triangles (i<j<k).
#
# Input:
#   - dconn  : 3D dynamic connectivity array [time × nROI × nROI]
#   - dtriad : 4D dynamic triad type tensor [time × nROI × nROI × nROI]
#
# Output:
#   - A 4×4 matrix of transition counts (no normalization, no diagonal edits).
#
# Notes:
#   - Triad types are assumed to be {+3, -1, +1, -3} in this fixed order.
# ================================================================

triad_transition_equivalent <- function(dtriad) {
  types <- c(+3, -1, +1, -3)
  K <- length(types)
  nROI <- dim(dtriad)[2]
  
  trans <- matrix(0L, nrow = K, ncol = K)
  
  for (i in 1:(nROI - 2L)) {
    for (j in (i + 1L):(nROI - 1L)) {
      for (k in (j + 1L):nROI) {
        ts <- dtriad[, i, j, k]            # triad labels over time (length T)
        # exact intent of ts1 <- ts[1:4749]; ts2 <- ts[-1] generalized:
        from <- ts[-length(ts)]
        to   <- ts[-1L]
        
        # Vectorized tabulation (assumes no NAs; identical to your sums if none)
        from_idx <- match(from, types)
        to_idx   <- match(to,   types)
        lin <- from_idx + (to_idx - 1L) * K
        counts <- tabulate(lin, nbins = K * K)
        trans <- trans + matrix(counts, nrow = K, ncol = K)
      }
    }
  }
  
  # Match original post-processing:
  diag(trans) <- NA
  trans <- trans / sum(trans, na.rm = TRUE)
  
  colnames(trans) <- paste0(types)
  rownames(trans) <- paste0(types)
  return(trans)
}
