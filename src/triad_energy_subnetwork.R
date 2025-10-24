# ================================================================
# File: triad_energy_subnetwork.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Project:
#   Empirical Evidence for Structural Balance Theory in Functional Brain Networks
#
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance Theory in Functional Brain Networks.
#
# Purpose:
#   Compute subnetwork MEAN peak ABSOLUTE triadic energy for each
#   triad SUM code using a POOLED aggregation strategy:
#     -> Consider ONLY triangles fully within `rois`
#     -> Pool ALL segments from those triangles for a given code
#     -> Take the mean of those segment energies
#
#   Output order matches lifetime summaries:
#     Pos3 (+++), Neg1 (-+-), Pos1 (+-+), Neg3 (---)
#
# Inputs:
#   de   : 3D array/list [n_ROI × n_ROI × n_ROI]
#          Output of DynTriadPeakEnergy(). Each cell is an RLE-like object:
#            values  = triad type code per segment (3, -1, 1, -3)
#            lengths = peak ABSOLUTE energy per segment
#          Supported shapes:
#            - list(values = <codes>, lengths = <energies>)
#            - list(list(<energies_vec>, <codes_vec>), ...)
#            - list(<energies_vec>, <codes_vec>)
#   rois : Integer vector of ROI indices defining the subnetwork
#
# Triad code convention (RAW SUMS, not remapped):
#      +3 : + + +     (balanced)
#      -1 : - + -     (balanced)
#      +1 : + - +     (imbalanced)
#      -3 : - - -     (imbalanced)
#
# Aggregation:
#   - POOLED mean within subnetwork:
#       mean( energy_s | code_s == code, triangle_s ⊆ rois )
#   - Triangles with no segments for a code do not contribute to that code.
#
# Output:
#   Named numeric vector of pooled mean peak energies within subnetwork:
#       c(Triad_+++ = mean_energy(+3),
#         Triad_-+- = mean_energy(-1),
#         Triad_+-+ = mean_energy(+1),
#         Triad_--- = mean_energy(-3))
#
# Notes:
#   - If a code never appears in the subnetwork, its mean is NA.
#   - Ensure upstream energies are ABSOLUTE if you expect non-negative values.
#
# Example:
#   # de produced by DynTriadPeakEnergy(...)
#   out <- triad_energy_subnetwork(de, rois = c(1,2,3,4,5))
#   print(out)
# ================================================================

triad_energy_subnetwork <- function(de, rois) {
  # ---- basic checks ---------------------------------------------------------
  if (missing(de) || is.null(de)) stop("`de` (DynTriadPeakEnergy output) is required.")
  if (missing(rois) || length(rois) < 3) stop("`rois` must contain at least 3 ROI indices.")
  
  rois <- sort(unique(as.integer(rois)))
  
  dm <- dim(de)
  if (length(dm) != 3) stop("`de` must be a 3D array/list [n_ROI × n_ROI × n_ROI].")
  if (!(dm[1] == dm[2] && dm[2] == dm[3])) {
    stop("`de` spatial dims must be a cube: dim(de) should be [n_ROI, n_ROI, n_ROI].")
  }
  nroi <- dm[1]
  if (any(rois < 1 | rois > nroi)) stop("`rois` contains indices out of range.")
  
  # ---- helper: extract (energies, codes) from a cell in 'de' ----------------
  .get_energy_codes <- function(x) {
    # Case A: named list with $lengths (energies) and $values (codes)
    if (is.list(x) && !is.null(x$lengths) && !is.null(x$values)) {
      return(list(en = as.numeric(x$lengths), code = as.integer(x$values)))
    }
    # Case B: nested list as in earlier pipelines: x[[1]][[1]] (energies), x[[1]][[2]] (codes)
    if (is.list(x) && length(x) >= 1 && is.list(x[[1]]) && length(x[[1]]) >= 2) {
      return(list(en = as.numeric(x[[1]][[1]]), code = as.integer(x[[1]][[2]])))
    }
    # Case C: simple 2-element list (energies, codes)
    if (is.list(x) && length(x) >= 2) {
      return(list(en = as.numeric(x[[1]]),     code = as.integer(x[[2]])))
    }
    # Fallback: empty
    return(list(en = numeric(0), code = integer(0)))
  }
  
  # ---- accumulate ALL segments across ALL triangles in the subnetwork -------
  all_codes  <- integer(0)
  all_energy <- numeric(0)
  
  L <- length(rois)
  for (a in 1:(L - 2)) {
    i <- rois[a]
    for (b in (a + 1):(L - 1)) {
      j <- rois[b]
      for (c in (b + 1):L) {
        k <- rois[c]
        
        ec <- .get_energy_codes(de[[i, j, k]])
        if (length(ec$en) > 0 && length(ec$code) > 0) {
          # append while preserving alignment
          all_energy <- c(all_energy, ec$en)
          all_codes  <- c(all_codes,  ec$code)
        }
      }
    }
  }
  
  # ---- pooled mean energy for a specific triad code -------------------------
  mean_energy <- function(code) {
    idx <- which(all_codes == code)
    if (length(idx) == 0) return(NA_real_)
    mean(all_energy[idx], na.rm = TRUE)
  }
  
  # ---- compute pooled means in the requested order --------------------------
  out <- c(
    Pos3 = mean_energy(+3),  # +++
    Neg1 = mean_energy(-1),  # -+-
    Pos1 = mean_energy(+1),  # +-+
    Neg3 = mean_energy(-3)   # ---
  )
  names(out) <- c("Triad_+++", "Triad_-+-", "Triad_+-+", "Triad_---")
  return(out)
}
