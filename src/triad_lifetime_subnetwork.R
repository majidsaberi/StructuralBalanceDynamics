# ================================================================
# File: triad_lifetime_subnetwork.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance Theory in Functional Brain Networks.
#
# Purpose:
#   Compute subnetwork MEAN lifetimes for each triad type using a
#   POOLED aggregation strategy:
#     -> Consider ONLY triangles fully within `rois`
#     -> Pool ALL runs from those triangles for a given code
#     -> Take the mean of those run-lengths
#
# Inputs:
#   lf : 3D array/list [n_ROI × n_ROI × n_ROI] of RLE triad-state sequences,
#        e.g. from LifeTime(dtriad, ...).
#        Each lf[[i,j,k]] must provide (values = codes, lengths = run-lengths).
#        Supported shapes:
#          - list(values = <codes>, lengths = <runs>)
#          - list(list(<lengths_vec>, <values_vec>), ...)
#          - list(<lengths_vec>, <values_vec>)
#
#   rois : integer vector of ROI indices (1-based). Only triangles with all
#          three vertices in `rois` are included.
#
#   window_step_seconds : numeric or NULL. If provided, run-lengths (in windows)
#          are multiplied by this factor to report lifetimes in seconds.
#
#   code_scheme : "remapped" or "raw"
#          - "remapped": {+3, -2, +1, -3}  with -2 for (-+-)
#          - "raw"     : {+3, -1, +1, -3}  with -1 for (-+-)
#
# Output:
#   Named numeric vector (Triad_+++, Triad_-+-, Triad_+-+, Triad_---) of pooled
#   mean lifetimes (in windows, or seconds if window_step_seconds provided).
#
# Notes:
#   - Returns NA for any triad type not present in the subnetwork.
#   - This is the pooled mean across ALL runs, not an average of per-triangle means.
#
# ================================================================

triad_lifetime_subnetwork <- function(lf,
                                      rois,
                                      window_step_seconds = NULL,
                                      code_scheme = c("remapped","raw")) {
  code_scheme <- match.arg(code_scheme)
  
  # --- basic checks ---
  if (!is.array(lf) || length(dim(lf)) != 3) {
    stop("lf must be a 3D array/list of RLE entries [n_ROI × n_ROI × n_ROI].")
  }
  if (length(rois) < 3L) {
    out <- rep(NA_real_, 4)
    names(out) <- c("Triad_+++", "Triad_-+-", "Triad_+-+", "Triad_---")
    return(out)
  }
  
  n_roi <- dim(lf)[1]
  if (any(rois < 1L | rois > n_roi)) {
    stop("All 'rois' must be valid indices within [1, n_ROI].")
  }
  rois <- sort(unique(as.integer(rois)))
  
  # --- small extractor to tolerate common RLE shapes ---
  .get_len_vals <- function(x) {
    # Case A: named list
    if (is.list(x) && !is.null(x$lengths) && !is.null(x$values)) {
      return(list(len = as.numeric(x$lengths), val = as.integer(x$values)))
    }
    # Case B: nested list as in some pipelines: x[[1]][[1]] (len), x[[1]][[2]] (val)
    if (is.list(x) && length(x) >= 1 && is.list(x[[1]]) && length(x[[1]]) >= 2) {
      return(list(len = as.numeric(x[[1]][[1]]), val = as.integer(x[[1]][[2]])))
    }
    # Case C: simple 2-element list (len, val)
    if (is.list(x) && length(x) >= 2) {
      return(list(len = as.numeric(x[[1]]),     val = as.integer(x[[2]])))
    }
    # Fallback: empty
    list(len = numeric(0), val = integer(0))
  }
  
  # --- accumulate ALL runs across ALL triangles in the subnetwork ---
  all_vals <- integer(0)
  all_len  <- numeric(0)
  
  L <- length(rois)
  for (a in 1:(L - 2)) {
    i <- rois[a]
    for (b in (a + 1):(L - 1)) {
      j <- rois[b]
      for (c in (b + 1):L) {
        k <- rois[c]
        r <- lf[[i, j, k]]
        lv <- .get_len_vals(r)
        if (length(lv$len) > 0 && length(lv$val) > 0) {
          all_vals <- c(all_vals, lv$val)
          all_len  <- c(all_len,  lv$len)
        }
      }
    }
  }
  
  if (length(all_vals) == 0L) {
    out <- rep(NA_real_, 4)
    names(out) <- c("Triad_+++", "Triad_-+-", "Triad_+-+", "Triad_---")
    return(out)
  }
  
  # --- optional conversion to seconds ---
  if (!is.null(window_step_seconds)) {
    if (!is.numeric(window_step_seconds) || length(window_step_seconds) != 1L || window_step_seconds <= 0) {
      stop("'window_step_seconds' must be a single positive number if provided.")
    }
    all_len <- all_len * window_step_seconds
  }
  
  # --- code set per scheme ---
  code_negposneg <- if (code_scheme == "remapped") -2L else -1L
  
  # --- pooled mean lifetime for a specific triad code ---
  mean_life <- function(code) {
    idx <- which(all_vals == code)
    if (length(idx) == 0L) return(NA_real_)
    mean(all_len[idx], na.rm = TRUE)
  }
  
  out <- c(
    mean_life(+3L),           # +++
    mean_life(code_negposneg),# -+-
    mean_life(+1L),           # +-+
    mean_life(-3L)            # ---
  )
  names(out) <- c("Triad_+++", "Triad_-+-", "Triad_+-+", "Triad_---")
  return(out)
}
