# ================================================================
# File: dyn_triad_type.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance Theory in Functional Brain Networks.
#
# Purpose:
#   Compute the raw signed triad sums over time from a dynamic
#   signed connectivity tensor. This reproduces the behavior of
#   the original DynTriad function (returns sij + sik + sjk).
#
# Input:
#   dconn : 3D numeric array [n_win × n_ROI × n_ROI]
#           (e.g., dynamic correlations). Values are binarized by sign.
#
# Sign binarization:
#   - values > 0  ->  +1
#   - values < 0  ->  -1
#   - values == 0 ->   0   (zeros are kept and contribute 0 to the sum)
#
# Output:
#   dtriad : 4D integer array [n_win × n_ROI × n_ROI × n_ROI]
#            Only upper-triangle triplets (i < j < k) are filled.
#            Possible values: {-3, -1, 1, 3}.
#
# Notes:
#   - This function DOES NOT classify or remap values; it returns raw sums.
#   - No NA is introduced for zero edges; zeros remain and affect the sum.
#   - Assumes dconn[t,,] is a (possibly symmetric) ROI×ROI matrix per window.
#   - Time/space complexity scales with O(n_win * n_ROI^3).
#
# ================================================================

dyn_triad_type <- function(dconn) {
  # --- validate input ---
  if (length(dim(dconn)) != 3) {
    stop("dconn must be a 3D array: [n_win × n_ROI × n_ROI].")
  }
  if (!is.numeric(dconn)) {
    stop("dconn must be numeric.")
  }
  
  # --- binarize signs (neg->-1, pos->+1, zero stays 0) ---
  x <- dconn
  x[x > 0] <-  1L
  x[x < 0] <- -1L
  # zeros remain 0L
  
  # --- dimensions ---
  dm    <- dim(x)              # [n_win, n_ROI, n_ROI]
  n_win <- dm[1]
  n_roi <- dm[2]
  
  # --- allocate output: only (i<j<k) entries are filled ---
  dtriad <- array(NA_integer_, dim = c(n_win, n_roi, n_roi, n_roi))
  
  # --- compute triad sums ---
  for (t in 1:n_win) {
    Xt <- x[t, , ]             # ROI × ROI at window t
    for (i in 1:(n_roi - 2)) {
      for (j in (i + 1):(n_roi - 1)) {
        sij <- Xt[i, j]
        for (k in (j + 1):n_roi) {
          dtriad[t, i, j, k] <- sij + Xt[i, k] + Xt[j, k]
        }
      }
    }
  }
  
  return(dtriad)
}
