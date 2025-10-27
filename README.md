# Dynamic Structural Balance Analysis in Functional Brain Networks

This repository implements the computational framework developed by  
**Majid Saberi et al. (2025)** for investigating *structural balance theory (SBT)* in dynamic brain networks.

The code estimates how patterns of positive and negative interactions among brain regions form and evolve over time,  
quantifying the **stability**, **tension**, and **transitions** of triadic relationships.

---

## 🧠 Background

Structural Balance Theory, originating in social network science, can be applied to brain connectivity to characterize  
how cooperative and antagonistic interactions organize the network over time.  
In this framework, each triplet of regions (a *triad*) can be:

| Triad | Configuration | Interpretation |
|-------|----------------|----------------|
| +++ | all positive | strongly coordinated (balanced) |
| −+− | two negative, one positive | balanced antagonism (“enemy of my enemy”) |
| +−+ | two positive, one negative | unstable configuration |
| −−− | all negative | mutual opposition (imbalanced) |

The project provides a dynamic approach for testing whether balanced triads are indeed more persistent than imbalanced ones—  
an empirical validation of a central prediction of SBT.

---

## ⚙️ Overview of the Computational Pipeline

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
        ├──► Triadic Energy (momentary tension)
        ├──► Peak Energy (within-lifetime maximum)
        ├──► Whole-Brain and Subnetwork Aggregates
        ├──► Surrogate Analysis (phase-randomized control)
        └──► Triad Transition Matrix (4×4 state dynamics)



# Dynamic Strutural Balance Analysis for Functional Brain Networks

This repository provides a complete and reproducible pipeline for computing **dynamic structural balance metrics** from functional brain time series. The workflow includes:

- Sliding-window **dynamic functional connectivity**
- **Signed triad determination** over time
- **Triad Lifetimes** and **Peak Aboslute Energy**
- **Subnetwork-restricted summaries**
- **Phase-randomized surrogate analysis** for null modeling
- **Triadic transition probability matrices**

This code reproduces the analyses used in:

> **Saberi, M. et al. (2025). _Empirical Evidence for Structural Balance Theory in Functional Brain Networks._**


## 📌 Repository Structure

```text
.
├── main.R
├── README.md
└── src
    ├── dyn_connectivity.R
    ├── dyn_triad_type.R
    ├── dyn_triad_lifetime.R
    ├── triad_lifetime_brain.R
    ├── triad_lifetime_subnetwork.R
    ├── dyn_triad_energy.R
    ├── dyn_triad_peak_energy.R
    ├── triad_energy_brain.R
    ├── triad_energy_subnetwork.R
    ├── surrogates.R
    └── triad_transition_matrix.R



✉️ Contact

Author: Majid Saberi
Email: majidsa@umich.edu
