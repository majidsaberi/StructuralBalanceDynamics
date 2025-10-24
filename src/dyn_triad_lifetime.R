# ================================================================
# File: dyn_triad_lifetime.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance Theory in Functional Brain Networks.
#
# Purpose:
#   Compute run-length encoded lifetimes of triad types over time.
#
# Input:
#   dtriad : 4D integer array [n_win × n_ROI × n_ROI × n_ROI]
#            Triad-type codes per time window (e.g., from DynTriad):
#              +3  : balanced  + + +
#              -2  : balanced  - + -
#              +1  : imbalanced + - +
#              -3  : imbalanced - - -
#            NA   : undefined (if any edge had zero sign)
#
# Parameters:
#   drop_na_runs : logical (default TRUE)
#                  If TRUE, remove runs where rle$values == NA.
#                  If FALSE, NA runs (if any) are kept in the sequence.
#
# Output:
#   lf : 3D array of lists [n_ROI × n_ROI × n_ROI]
#        For i<j<k, lf[[i, j, k]] is an rle object with:
#          - $values  : triad-type codes in order of appearance
#          - $lengths : number of consecutive windows in each run
#        (Entries where i>=j or j>=k remain NULL.)
#
# Notes:
#   - This function does not aggregate; it stores full run-level detail.
#   - Use a separate summarizer to compute mean lifetimes per triad type.
# ================================================================

dyn_triad_lifetime <- function(dtriad, drop_na_runs = TRUE) {
  # --- validate input ---
  if (length(dim(dtriad)) != 4) {
    stop("dtriad must be a 4D array: [n_win × n_ROI × n_ROI × n_ROI].")
  }
  if (!is.numeric(dtriad) && !is.integer(dtriad)) {
    stop("dtriad must be numeric/integer triad codes.")
  }
  
  dm    <- dim(dtriad)         # [n_win, n_ROI, n_ROI, n_ROI]
  n_win <- dm[1]
  n_roi <- dm[2]
  
  # Pre-allocate a 3D array of lists
  lf <- array(vector("list", n_roi * n_roi * n_roi),
              dim = c(n_roi, n_roi, n_roi))
  
  # Compute run-length encoding for each (i,j,k) with i<j<k
  for (i in 1:(n_roi - 2)) {
    for (j in (i + 1):(n_roi - 1)) {
      for (k in (j + 1):n_roi) {
        ts <- dtriad[, i, j, k]
        
        # rle on integer/numeric vector (NAs allowed)
        r <- rle(ts)
        
        if (drop_na_runs) {
          keep <- !is.na(r$values)
          r$values  <- r$values[keep]
          r$lengths <- r$lengths[keep]
        }
        
        lf[[i, j, k]] <- r
      }
    }
  }
  
  return(lf)
}
