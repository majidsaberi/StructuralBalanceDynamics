# ================================================================
# File: main.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
# 
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance
# Theory in Functional Brain Networks.
#
# Purpose:
#   Demonstrate the full dynamic triad pipeline on synthetic data and
#   phase-randomized surrogate time series, including:
#     - Dynamic connectivity
#     - Triad types, lifetimes, and energy
#     - Subnetwork summaries
#     - Triadic transition matrices
#
# ================================================================

message("== Dynamic triad demo with surrogates + transitions: start ==")
rm(list = ls())
setwd("/Users/majid/Projects/Balance-Validation/Paper/Submission/FNP/Codes/")

# ---- Source files ------------------------------------------------------------
source("src/dyn_connectivity.R")          # compute dynamic FC matrices (windowed correlations)
source("src/dyn_triad_type.R")            # compute signed triad types over time
source("src/dyn_triad_lifetime.R")        # compute RLE lifetimes for each triangle
source("src/triad_lifetime_brain.R")      # whole-brain pooled mean lifetimes for triadic types
source("src/triad_lifetime_subnetwork.R") # subnetwork-only lifetime summary for triadic types
source("src/dyn_triad_energy.R")          # instantaneous triadic energy over time
source("src/dyn_triad_peak_energy.R")     # peak |energy| per lifetime segment
source("src/triad_energy_brain.R")        # whole-brain pooled mean peak energies for triadic types
source("src/triad_energy_subnetwork.R")   # subnetwork-only mean peak energies for triadic types
source("src/surrogates.R")                # generate phase–randomized surrogate time series
source("src/triad_transition_matrix.R")   # compute 4×4 triadic transition probability matrix

# ---- Create a dummy ROI×TIME matrix -----------------------------------------
set.seed(42)
nROI  <- 20
nTIME <- 180
TS <- matrix(rnorm(nROI * nTIME), nrow = nROI, ncol = nTIME)

# ---- Parameters --------------------------------------------------------------
WL <- 30  # note: your dyn_connectivity uses window i:(i+WL) -> effective length WL+1

# ---- Run dynamic connectivity ------------------------------------------------
message("\n== Step 1: dyn_connectivity ==")
dconn <- dyn_connectivity(TS, WL = WL, method = "pearson", use = "pairwise.complete.obs")
message("dconn dims [n_win × nROI × nROI] = ", paste(dim(dconn), collapse = " × "))

# ---- Triad type over time ----------------------------------------------------
message("\n== Step 2: dyn_triad_type ==")
dtriad <- dyn_triad_type(dconn)
message("dtriad dims [n_win × nROI × nROI × nROI] = ", paste(dim(dtriad), collapse = " × "))

# ---- Lifetime RLE per triangle ----------------------------------------------
message("\n== Step 3: dyn_triad_lifetime ==")
lf <- dyn_triad_lifetime(dtriad, drop_na_runs = TRUE)
message("lf is a 3D array of lists with dims: ", paste(dim(lf), collapse = " × "))

# ---- Whole-brain pooled mean lifetimes --------------------------------------
message("\n== Step 4: triad_lifetime_brain ==")
life_brain <- triad_lifetime_brain(lf)
print(round(life_brain, 3))

# ---- Triad energy over time --------------------------------------------------
message("\n== Step 5: dyn_triad_energy ==")
de_time <- dyn_triad_energy(dconn)
message("triad energy dims: ", paste(dim(de_time), collapse = " × "))

# ---- Peak |energy| per lifetime segment -------------------------------------
message("\n== Step 6: dyn_triad_peak_energy ==")
lf_peak <- dyn_triad_peak_energy(de_time, lf)
message("lf_peak shares lf's structure; 'lengths' now hold peak |energy| per segment.")

# ---- Whole-brain pooled mean peak |energy| ----------------------------------
message("\n== Step 7: triad_energy_brain ==")
energy_brain <- triad_energy_brain(lf_peak)
print(round(energy_brain, 3))

# ---- Subnetwork examples -----------------------------------------------------
sub_rois <- 1:6
message("\n== Step 8: Subnetwork pooled means (ROIs: ", paste(sub_rois, collapse = ","), ") ==")
life_sub   <- triad_lifetime_subnetwork(lf, rois = sub_rois, window_step_seconds = NULL, code_scheme = "raw")
energy_sub <- triad_energy_subnetwork(lf_peak, rois = sub_rois)

message("\nSubnetwork pooled mean lifetimes (windows):")
print(round(life_sub, 3))
message("\nSubnetwork pooled mean peak |energy|:")
print(round(energy_sub, 3))

# ---- SURROGATE: build surrogate matrix and run pipeline ---------------------
message("\n== SURROGATE: building phase-randomized matrix ==")
set.seed(4242)  # separate seed for surrogates
TS_surr <- TS
for (r in 1:nROI) TS_surr[r, ] <- surrogates(TS[r, ])

message("\n== SURROGATE: dyn_connectivity -> triad metrics ==")
dconn_s <- dyn_connectivity(TS_surr, WL = WL, method = "pearson", use = "pairwise.complete.obs")
dtriad_s <- dyn_triad_type(dconn_s)
lf_s <- dyn_triad_lifetime(dtriad_s, drop_na_runs = TRUE)
life_brain_s <- triad_lifetime_brain(lf_s)
de_time_s <- dyn_triad_energy(dconn_s)
lf_peak_s <- dyn_triad_peak_energy(de_time_s, lf_s)
energy_brain_s <- triad_energy_brain(lf_peak_s)

# ---- Subnetwork example (same subset for both) -------------------------------
sub_rois <- 1:6
life_sub     <- triad_lifetime_subnetwork(lf,    rois = sub_rois, window_step_seconds = NULL, code_scheme = "raw")
life_sub_s   <- triad_lifetime_subnetwork(lf_s,  rois = sub_rois, window_step_seconds = NULL, code_scheme = "raw")
energy_sub   <- triad_energy_subnetwork(lf_peak,   rois = sub_rois)
energy_sub_s <- triad_energy_subnetwork(lf_peak_s, rois = sub_rois)

# ---- Triadic transition matrices (final example) -----------------------------
# Triad types assumed order: {+3, -1, +1, -3}; function returns normalized counts with NA diagonal
message("\n== Triadic transition matrices (row=from, col=to; diagonal set to NA) ==")
trans_orig <- triad_transition_equivalent(dtriad)
trans_surr <- triad_transition_equivalent(dtriad_s)

print(round(trans_orig, 4))
cat("\n")
print(round(trans_surr, 4))

# ---- Compact comparisons for lifetimes & energies ----------------------------
fmt <- function(x) sprintf("%7.3f", x)

message("\n== Whole-brain pooled MEAN lifetimes (windows) ==")
lab <- names(life_brain)
cat(paste0("Triad      |  Original  | Surrogate |  Δ (S - O)\n",
           "-----------+------------+-----------+-----------\n"))
for (i in seq_along(lab)) {
  o <- life_brain[i]; s <- life_brain_s[i]
  cat(sprintf("%-10s| %s | %s | %s\n", lab[i], fmt(o), fmt(s), fmt(s - o)))
}

message("\n== Whole-brain pooled MEAN peak |energy| ==")
labE <- names(energy_brain)
cat(paste0("Triad      |  Original  | Surrogate |  Δ (S - O)\n",
           "-----------+------------+-----------+-----------\n"))
for (i in seq_along(labE)) {
  o <- energy_brain[i]; s <- energy_brain_s[i]
  cat(sprintf("%-10s| %s | %s | %s\n", labE[i], fmt(o), fmt(s), fmt(s - o)))
}

message("\n== SUBNETWORK (ROIs 1–6) pooled MEAN lifetimes (windows) ==")
labS <- names(life_sub)
cat(paste0("Triad      |  Original  | Surrogate |  Δ (S - O)\n",
           "-----------+------------+-----------+-----------\n"))
for (i in seq_along(labS)) {
  o <- life_sub[i]; s <- life_sub_s[i]
  cat(sprintf("%-10s| %s | %s | %s\n", labS[i], fmt(o), fmt(s), fmt(s - o)))
}

message("\n== SUBNETWORK (ROIs 1–6) pooled MEAN peak |energy| ==")
labSE <- names(energy_sub)
cat(paste0("Triad      |  Original  | Surrogate |  Δ (S - O)\n",
           "-----------+------------+-----------+-----------\n"))
for (i in seq_along(labSE)) {
  o <- energy_sub[i]; s <- energy_sub_s[i]
  cat(sprintf("%-10s| %s | %s | %s\n", labSE[i], fmt(o), fmt(s), fmt(s - o)))
}

