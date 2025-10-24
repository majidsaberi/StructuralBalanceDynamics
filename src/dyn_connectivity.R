# ================================================================
# File: dyn_connectivity.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance Theory in Functional Brain Networks.
#
# Purpose:
#   Compute dynamic (time-resolved) functional connectivity (FC)
#   matrices using a sliding-window correlation approach.
#
# Description:
#   Given a regional time series matrix (ROI × Time), this function
#   divides the time series into overlapping windows of length WL
#   and computes a correlation matrix for each window.
#
#   The output is a 3D array of shape [n_windows × n_ROI × n_ROI],
#   where each slice along the first dimension corresponds to the
#   FC matrix at that time window.
#
# Input:
#   TS : Numeric matrix of regional time series (rows = ROIs, columns = timepoints)
#   WL : Integer, sliding window length in timepoints (default = 50)
#   method : Correlation method ("pearson", "spearman", etc.; default = "pearson")
#   use : NA handling strategy for correlation ("pairwise.complete.obs" recommended)
#
# Output:
#   A 3D array containing dynamic FC matrices:
#       [time_window × ROI × ROI]
#
# ================================================================

dyn_connectivity <- function(TS, WL = 50, method = "pearson", use = "pairwise.complete.obs") {
  
  # --- Validate input ---------------------------------------------------------
  if (!is.matrix(TS)) {
    stop("Input TS must be a numeric matrix: rows = ROIs, columns = timepoints.")
  }
  if (!is.numeric(TS)) stop("TS must be numeric.")
  if (WL < 2) stop("Window length WL must be at least 2 timepoints.")
  
  n_roi  <- nrow(TS)
  n_time <- ncol(TS)
  n_win  <- n_time - WL 
  
  if (n_win < 1) {
    stop("Window length WL is longer than the available time series length.")
  }
  
  # --- Initialize output array -----------------------------------------------
  dconn <- array(NA_real_, dim = c(n_win, n_roi, n_roi))
  
  # --- Compute windowed correlations -----------------------------------------
  message("Computing dynamic connectivity across ", n_win, " windows...")
  for (i in seq_len(n_win)) {
    # extract window segment
    seg <- TS[, i:(i + WL ), drop = FALSE]
    
    # compute ROI×ROI correlation matrix
    dconn[i, , ] <- stats::cor(t(seg), method = method, use = use)
    
    # optional progress message every 50 windows
    if (i %% 50 == 0) message("Processed window ", i, "/", n_win)
  }
  
  # --- Return results --------------------------------------------------------
  message("Dynamic connectivity computation completed.")
  return(dconn)
}
