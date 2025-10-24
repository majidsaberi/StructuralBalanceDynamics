# ================================================================
# File: triad_energy_brain.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance Theory in Functional Brain Networks.
#
# Purpose:
#   Compute whole-brain MEAN peak ABSOLUTE triadic energy for each
#   triad SUM code using a POOLED aggregation strategy:
#     -> Pool ALL segments from ALL triangles for a given code
#     -> Take the mean of those segment energies
#
#   Order of the output matches your lifetime summaries:
#     Pos3 (+++) , Neg1 (-+-) , Pos1 (+-+) , Neg3 (---)
#
# Input:
#   de : 3D array/list [n_ROI × n_ROI × n_ROI] where each cell is an
#        RLE-like object for DynTriadPeakEnergy(), supporting any of:
#          - list(values = <triad codes>, lengths = <peak energies>)
#          - list(list(<lengths_vec>, <values_vec>), ...)
#          - list(<lengths_vec>, <values_vec>)
#        Here, "lengths" holds the peak ABSOLUTE energy per segment and
#        "values" holds the triad code for that segment.
#
# Triad code convention (RAW SUMS, not remapped):
#      +3 : + + +     (balanced)
#      -1 : - + -     (balanced)
#      +1 : + - +     (imbalanced)
#      -3 : - - -     (imbalanced)
#
# Aggregation:
#   - Pooled mean: mean( energy_s | code_s == code ) over ALL segments s
#   - Triangles with no segments for a code do not contribute to that code
#
# Output:
#   Named numeric vector (in the order below) of pooled mean peak energies:
#       c(Triad_+++ = mean_energy(+3),
#         Triad_-+- = mean_energy(-1),
#         Triad_+-+ = mean_energy(+1),
#         Triad_--- = mean_energy(-3))
#
# Notes:
#   - If a code never appears in the entire dataset, its mean is NA.
#   - Ensure DynTriadPeakEnergy() produced ABSOLUTE energies upstream
#     if you expect strictly non-negative values here.
#
# ================================================================

triad_energy_brain <- function(de) {
  # --- helper: extract (energies, codes) from a cell in 'de' ---
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
  
  # --- accumulate ALL segments across ALL triangles ---
  all_codes <- integer(0)
  all_energy <- numeric(0)
  
  dm <- dim(de)
  for (i in seq_len(dm[1])) {
    for (j in seq_len(dm[2])) {
      for (k in seq_len(dm[3])) {
        ec <- .get_energy_codes(de[[i, j, k]])
        if (length(ec$en) > 0 && length(ec$code) > 0) {
          # append while preserving alignment
          all_energy <- c(all_energy, ec$en)
          all_codes  <- c(all_codes,  ec$code)
        }
      }
    }
  }
  
  # --- pooled mean energy for a specific triad code ---
  mean_energy <- function(code) {
    idx <- which(all_codes == code)
    if (length(idx) == 0) return(NA_real_)
    mean(all_energy[idx], na.rm = TRUE)
  }
  
  # --- compute pooled means in the requested order ---
  out <- c(
    mean_energy(+3),  # Pos3: +++
    mean_energy(-1),  # Neg1: -+-
    mean_energy(+1),  # Pos1: +-+
    mean_energy(-3)   # Neg3: ---
  )
  
  names(out) <- c("Triad_+++", "Triad_-+-", "Triad_+-+", "Triad_---")
  return(out)
}
