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

---

## 📌 Repository Structure

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
