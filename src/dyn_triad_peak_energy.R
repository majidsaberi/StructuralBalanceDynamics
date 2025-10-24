# ================================================================
# File: dyn_triad_peak_energy.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance Theory in Functional Brain Networks.
#
# Purpose:
#   For each triangle (i<j<k) and each contiguous lifetime segment
#   (as returned by LifeTime()/rle on triad type), compute the
#   peak absolute triadic energy within that segment.
#
# Inputs:
#   dtriad_energy : 4D array [time × nROI × nROI × nROI]
#                   Output of DynTriadEnergy(), containing the
#                   momentary triadic energy for (i,j,k) at time t.
#
#   lf            : 3D list [nROI × nROI × nROI]
#                   Output of LifeTime(); for each (i,j,k),
#                   lf[[i,j,k]] is an RLE object with:
#                     - values  : triad-type code over segments
#                     - lengths : segment lengths (in time windows)
#                   This function OVERWRITES 'lengths' with the
#                   peak absolute energy of that segment.
#
# Output:
#   A 3D list with the same structure as 'lf', where for each
#   (i,j,k), lf[[i,j,k]]$lengths holds the peak |energy| within
#   each corresponding lifetime segment.
#
# Notes:
#   - Uses max(abs(energy), na.rm = TRUE) within each segment.
#   - If a segment is all NA, the peak is set to NA_real_.
#   - Keeps the original segment boundaries from the RLE.
#
# ================================================================

dyn_triad_peak_energy <- function(dtriad_energy, lf) {
  # Helper: safe peak of absolute values with NA handling
  peak_abs <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    return(max(abs(x), na.rm = TRUE))
  }
  
  dm <- dim(dtriad_energy)  # [time × nROI × nROI × nROI]
  if (length(dm) != 4) {
    stop("dtriad_energy must be a 4D array: [time × nROI × nROI × nROI].")
  }
  
  # Copy structure so we don't mutate the input 'lf' object by reference
  de <- lf
  
  nROI <- dm[2]
  
  for (i in 1:(nROI - 2)) {
    for (j in (i + 1):(nROI - 1)) {
      for (k in (j + 1):nROI) {
        # Time series of triadic energy for this triangle
        e_t <- dtriad_energy[, i, j, k]
        
        # RLE of triad types (as produced by LifeTime())
        rle_obj <- lf[[i, j, k]]
        if (is.null(rle_obj) || !is.list(rle_obj) ||
            !all(c("lengths", "values") %in% names(rle_obj))) {
          # If malformed, skip gracefully
          de[[i, j, k]] <- rle_obj
          next
        }
        
        # Walk through each lifetime segment and compute peak |energy|
        start_idx <- 1L
        for (seg in seq_along(rle_obj$lengths)) {
          seg_len <- rle_obj$lengths[seg]
          end_idx <- start_idx + seg_len - 1L
          
          # Guard for bounds
          if (end_idx > length(e_t)) end_idx <- length(e_t)
          if (start_idx <= end_idx) {
            seg_vals <- e_t[start_idx:end_idx]
            de[[i, j, k]]$lengths[seg] <- peak_abs(seg_vals)
          } else {
            de[[i, j, k]]$lengths[seg] <- NA_real_
          }
          
          start_idx <- end_idx + 1L
        }
      }
    }
  }
  
  return(de)
}
