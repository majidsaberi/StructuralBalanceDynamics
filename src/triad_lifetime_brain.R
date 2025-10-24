# ================================================================
# File: triad_lifetime_brain.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance Theory in Functional Brain Networks.
#
# Purpose:
#   Compute whole-brain MEAN triad lifetimes (in windows) for each
#   triad SUM type using a POOLED aggregation strategy:
#     -> Pool ALL runs from ALL triangles for a given code
#     -> Take the mean of those run-lengths
#
#   This estimates the expected lifetime of each triad state across
#   the whole brain by weighting every observed run equally.
#
# Input:
#   lf : 3D array [n_ROI × n_ROI × n_ROI] of run-length-encoded (RLE)
#        triad-state sequences per triangle (i,j,k). Each cell may be:
#          - list(values = <vector>, lengths = <vector>)
#          - list(list(<lengths_vec>, <values_vec>), ...)
#          - list(<lengths_vec>, <values_vec>)
#
# Triad code convention (RAW SUMS, not remapped):
#      +3 : + + +     (balanced)
#      -1 : - + -     (balanced)
#      +1 : + - +     (imbalanced)
#      -3 : - - -     (imbalanced)
#
# Aggregation:
#   - Pooled mean: mean( length_r | value_r == code ) over ALL runs r
#   - Triangles with no runs for a code do not contribute (ignored)
#
# Output:
#   Named numeric vector of pooled mean lifetimes 
#       c(Triad_+++ = mean_energy(+3),
#         Triad_-+- = mean_energy(-1),
#         Triad_+-+ = mean_energy(+1),
#         Triad_--- = mean_energy(-3))
#
# Notes:
#   - Multiply by TR or window step to convert to seconds if desired.
#
# ================================================================

triad_lifetime_brain <- function(lf) {
  # --- helper: extract (lengths, values) from a cell in lf ---
  .get_len_vals <- function(x) {
    # Case A: named list with $lengths and $values
    if (is.list(x) && !is.null(x$lengths) && !is.null(x$values)) {
      return(list(len = x$lengths, val = x$values))
    }
    # Case B: nested list as in original (x[[1]][[1]], x[[1]][[2]])
    if (is.list(x) && length(x) >= 1 && is.list(x[[1]]) && length(x[[1]]) >= 2) {
      return(list(len = x[[1]][[1]], val = x[[1]][[2]]))
    }
    # Case C: simple 2-element list (lengths, values)
    if (is.list(x) && length(x) >= 2) {
      return(list(len = x[[1]], val = x[[2]]))
    }
    # Fallback: empty
    return(list(len = numeric(0), val = numeric(0)))
  }
  
  # --- accumulate ALL runs across ALL triangles ---
  all_vals <- integer(0)
  all_len  <- numeric(0)
  
  dm <- dim(lf)
  for (i in seq_len(dm[1])) {
    for (j in seq_len(dm[2])) {
      for (k in seq_len(dm[3])) {
        lv <- .get_len_vals(lf[[i, j, k]])
        if (length(lv$len) > 0 && length(lv$val) > 0) {
          # append, keeping alignment between lengths and values
          all_vals <- c(all_vals, as.integer(lv$val))
          all_len  <- c(all_len,  as.numeric(lv$len))
        }
      }
    }
  }
  
  # --- helper: pooled mean lifetime for a specific triad code ---
  mean_life <- function(code) {
    idx <- which(all_vals == code)
    if (length(idx) == 0) return(NA_real_)
    mean(all_len[idx], na.rm = TRUE)
  }
  
  # --- compute pooled means in your preferred label order ---
  out <- c(
    mean_life(+3),  # Pos3:  +++
    mean_life(-1),  # Neg1:  -+-
    mean_life(+1),  # Pos1:  +-+
    mean_life(-3)   # Neg3:  ---
  )
  
  names(out) <- c("Triad_+++", "Triad_-+-", "Triad_+-+", "Triad_---")
  return(out)
}
