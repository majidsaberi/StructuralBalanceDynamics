# Dynamic Structural Balance Analysis in Functional Brain Networks

This repository implements the computational framework developed for investigating *Structural Balance Theory (SBT)* in **dynamic brain networks**  (**Saberi, M. et al., 2025**).

The code estimates how patterns of positive and negative interactions among brain regions form and evolve over time,  
quantifying the **stability**, **tension**, and **transitions** of triadic relationships.

This framework extends SBT to brain networks to address a fundamental open question that, until now, had not been empirically tested in the brain:  
> Are balanced triads genuinely stable, and are imbalanced triads transient?

This project introduces and implements a **dynamic framework** that directly tests these assumptions using resting-state fMRI data.  
It provides all computational tools required to reproduce the analyses described in the manuscript, including computation of **triadic lifetimes**, **peak energies**, **surrogate modeling**, and **triadic transition mapping**.

---

## Background

Structural Balance Theory, originally developed in social network science, can be applied to brain connectivity to characterize  
how cooperative and antagonistic interactions among brain regions organize the network over time.  
Each triplet of regions (a *triad*) can be classified as follows:

| Triad | Configuration | Interpretation |
|-------|----------------|----------------|
| +++ | all positive | Balanced |
| −+− | two negative, one positive | Balanced |
| +−+ | two positive, one negative | Imbalanced |
| −−− | all negative | Imbalanced |

This framework provides a dynamic approach for testing whether **balanced triads persist longer** than **imbalanced ones** —  
an empirical test of a central prediction of Structural Balance Theory.
---

## Overview of the Computational Pipeline

```text
Time Series (ROI × Time)
        │
        ▼
Sliding-Window Correlation → Dynamic Connectivity (DFC)
        │
        ▼
Triad Classification (+3, −1, +1, −3)
        │
        ├──► Triad Lifetimes (temporal stability)
        ├──► Triadic Absolute Peak Energy (momentary tension)
        ├──► Whole-Brain and Subnetwork Aggregates
        ├──► Surrogate Analysis (phase-randomized control)
        └──► Triad Transition Matrix (4×4 state dynamics)



## 📌 Repository Structure

```text
.
├── main.R                         # Main demo pipeline (synthetic + surrogate data)
├── README.md                      # Project documentation
└── src
    ├── dyn_connectivity.R         # Compute sliding-window dynamic functional connectivity
    ├── dyn_triad_type.R           # Derive triad-type tensors (+3, −1, +1, −3)
    ├── dyn_triad_lifetime.R       # Run-length encoding of triad lifetimes
    ├── triad_lifetime_brain.R     # Compute whole-brain mean lifetimes across all triads
    ├── triad_lifetime_subnetwork.R# Lifetime summary restricted to specific ROI subnetworks
    ├── dyn_triad_energy.R         # Calculate dynamic triadic energy (structural tension)
    ├── dyn_triad_peak_energy.R    # Extract peak |energy| within each lifetime segment
    ├── triad_energy_brain.R       # Compute whole-brain mean peak energies
    ├── triad_energy_subnetwork.R  # Subnetwork-level mean peak energies
    ├── surrogates.R               # Generate phase-randomized surrogate time series
    └── triad_transition_matrix.R  # Build 4×4 triad-type transition probability matrix



Please cite the below paper if you use them:

> **Saberi, M. et al. (2025). _Empirical Evidence for Structural Balance Theory in Functional Brain Networks._**

✉️ Contact

Author: Majid Saberi
Email: majidsa@umich.edu
