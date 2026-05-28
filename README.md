# Censored Data Simulation Framework for Estimator and CUSUM performance

This repository contains the R files for the simulations. The code evaluates both estimation accuracy and CUSUM control charts performance under varying levels of left-censoring.

---

## Overview

The simulation study evaluates how different approaches for handling censored observations affect:

- Phase I estimation of process parameters (mean and standard deviation)
- Phase II CUSUM monitoring performance
- 
Two main simulation components are included:

1. Estimator comparison study (Bias & RMSE analysis)
2. CUSUM performance study (ARL analysis)

---

## Methods Compared

### Phase I estimation methods

- Maximum Likelihood Estimation (MLE) via `survreg` (Gaussian censored model)
- Regression on Order Statistics (ROS) using `NADA`
- Kaplan–Meier (KM) estimator
- Substitution methods:
  - LOD substitution
  - LOD / √2 substitution

---

### Phase II monitoring methods

- Standard CUSUM (`qcr`)
- Likelihood-based CUSUM (log-likelihood formulation)
- Imputation-based methods:
  - LOD substitution
  - LOD / √2 substitution
  -  Random stochastic imputation

---

## Simulation Design

### Data generation

- Normal distribution on latent scale
- Left censoring at fixed LOD
- Censoring levels:
  0%, 25%, 50%, 75%, 85%, 95%

### Sample sizes

- Phase I: 200 observations
- Phase II:
  - 25 observations (in-control)
  - 3000 observations (post-shift ARL evaluation)

### Process shifts

- Δ = 0 (in-control)
- Δ = 0.5σ
- Δ = 1σ

---

## Key Outputs

### 1. Estimation performance

The simulation evaluates:

- Bias and RMSE of mean estimators
- Bias and RMSE of standard deviation estimators

Methods compared:
- MLE
- ROS
- KM
- Substitution methods

---

### 2. CUSUM performance

The CUSUM study evaluates:

- Average Run Length (ARL)
- Median run length
- Standard deviation of run length (SDRL)
- Empirical quantiles (10%, 25%, 75%, 90%)

across:
- censoring levels
- estimation methods
- process shifts

---


## Dependencies

- NADA
- survival
- qcr
- ggplot2
- dplyr
- tidyr

---

## Purpose

This code supports methodological research on:

- SPC monitoring under incomplete data
- Comparison of classical vs likelihood-based approaches
- Performance evaluation of CUSUM under censoring

---

## Author note

This framework is intended for academic research and thesis reproducibility. All parameters (censoring, shifts, sample sizes) can be modified for extended analysis.
