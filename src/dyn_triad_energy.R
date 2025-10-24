# ================================================================
# File: dyn_triad_energy.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance Theory in Functional Brain Networks.
#
# Purpose:
#   Compute the dynamic triadic energy for every i<j<k triangle
#   across sliding-window time points. Energy is defined as:
#
#       E = sign(a) * |a|^(1/3),
#       where a = - S_ij * S_ik * S_jk
#
#   This formulation follows structural balance theory, where
#   lower (more negative) energy indicates higher structural tension
#   and positive energy reflects internally consistent interactions.
#
# Input:
#   dconn : 3D array [time × nROI × nROI]
#           produced by dyn_connectivity(), containing dynamic
#           connectivity matrices (not yet thresholded to ±1).
#
# Output:
#   4D array [time × nROI × nROI × nROI]
#       dtriad_energy[t,i,j,k] = triadic energy of triangle (i,j,k)
#       at time t. Entries are NA for i≥j≥k (since only i<j<k triangles exist).
#
# Notes:
#   - This dynamic energy reflects momentary triadic tension/consistency.
#   - Works with weighted correlations before sign-thresholding.
#   - Used later for computing peak triad energy per lifetime.
#
# ================================================================

dyn_triad_energy <- function(dconn) {
  dm <- dim(dconn)  # [time × nROI × nROI]
  dtriad_energy <- array(NA, dim = c(dm[1], dm[2], dm[2], dm[2]))
  
  for (t in 1:dm[1]) {
    for (i in 1:(dm[2] - 2)) {
      for (j in (i + 1):(dm[2] - 1)) {
        for (k in (j + 1):dm[2]) {
          
          # Compute energy term: a = - S_ij * S_ik * S_jk
          a <- - dconn[t, i, j] * dconn[t, i, k] * dconn[t, j, k]
          
          # Normalize using cubic root, preserving sign
          dtriad_energy[t, i, j, k] <- sign(a) * abs(a)^(1/3)
        }
      }
    }
  }
  
  return(dtriad_energy)
}
