# ================================================================
# File: surrogates.R
# Author: Majid Saberi
# Email: majidsa@umich.edu
# Date: October 2025
#
# Citation:
# Saberi, M. et al. (2025). Empirical Evidence for Structural Balance
# Theory in Functional Brain Networks.
#
# Purpose:
#   Generate phase-randomized surrogate time series for null-model
#   construction in dynamic brain network analysis. The surrogate
#   preserves the original amplitude spectrum while destroying temporal
#   phase structure, enabling comparison against networks expected by chance.
#
# Method Summary:
#   1) Apply FFT to original time series.
#   2) Retain the amplitude spectrum (preserves power distribution).
#   3) Randomize phases of positive-frequency components uniformly
#      in (-pi, pi], and mirror them to negative frequencies to enforce
#      Hermitian symmetry (ensuring a real-valued inverse FFT).
#   4) Keep DC component (index 1) unchanged.
#   5) For even-length signals, randomize Nyquist component phase.
#   6) Apply inverse FFT to obtain a real surrogate time series.
#
# Input:
#   x - Numeric vector representing the original brain time series
#
# Output:
#   Numeric vector of same length as input, representing the phase-
#   randomized surrogate time series.
#
# Properties:
#   - Preserves the power spectrum and autocorrelation structure
#   - Eliminates meaningful temporal phase relationships
#   - Maintains statistical properties of the original signal
#   - Produces real-valued output suitable for further analysis
#
# Notes:
#   - Fully consistent with the Methods section in the manuscript.
#   - Based on established surrogate data methods for nonlinear time series analysis
#   - Suitable for testing significance of functional connectivity measures
#
# ================================================================

surrogates <- function(x) {
  # Convert input to numeric vector to ensure proper processing
  x <- as.numeric(x)
  Tn <- length(x)
  
  # Handle very short time series (edge case)
  if (Tn < 2) return(x)
  
  # Step 1: Compute Fast Fourier Transform of original time series
  X <- fft(x)
  amp <- Mod(X)  # Extract amplitude spectrum
  phi <- Arg(X)  # Extract phase spectrum (for reference)
  
  # Identify frequency indices
  half <- floor(Tn / 2)  # Index for positive frequency boundary
  
  # Step 2: Randomize phases for positive frequency components
  if (half > 1) {
    # Generate random phases uniformly distributed in (-pi, pi]
    rand_phi <- runif(half - 1, min = -pi, max = pi)
    
    # Apply random phases to positive frequencies (indices 2 to half)
    for (k in 2:half) {
      X[k] <- amp[k] * exp(1i * rand_phi[k-1])
    }
    
    # Step 3: Enforce Hermitian symmetry by mirroring to negative frequencies
    # This ensures the inverse FFT produces a real-valued signal
    for (k in 2:half) {
      X[Tn - k + 2] <- Conj(X[k])  # Negative frequency components
    }
  }
  
  # Step 4: Handle Nyquist component for even-length signals
  if (Tn %% 2 == 0) {
    nyq_phi <- runif(1, min = -pi, max = pi)
    X[half + 1] <- amp[half + 1] * exp(1i * nyq_phi)
  }
  
  # Note: DC component (index 1) remains unchanged to preserve mean
  
  # Step 5: Apply inverse FFT and return real part
  # Division by Tn is required for proper scaling in R's FFT implementation
  Re(fft(X, inverse = TRUE) / Tn)
}
