# Censored Data Simulation Framework for SPC, CUSUM, and Estimator Performance

This repository contains an R-based simulation framework developed for a master’s thesis on statistical methods for censored data. The code evaluates both estimation accuracy and statistical process control (SPC) performance under varying levels of left-censoring.

The focus is on applications in pharmaceutical, microbiological, and environmental data, where observations fall below a limit of detection (LOD).

---

## Overview

The simulation study evaluates how different approaches for handling censored observations affect:

- Phase I estimation of process parameters (mean and standard deviation)
- Bias and RMSE of estimators under censoring
- Phase II CUSUM monitoring performance
- Robustness of classical vs likelihood-based methods

Two main simulation components are included:

1. Estimator comparison study (Bias & RMSE analysis)
2. CUSUM / SPC performance study (ARL analysis)

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
  - Stochastic truncated normal imputation

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

- Bias of mean estimators
- Bias of standard deviation estimators
- RMSE of all methods

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

## Key assumptions

- Data are normally distributed on latent scale
- Censoring is left censoring at fixed LOD
- Phase I and Phase II share identical parameters unless a shift is introduced
- Censoring is non-informative

---

## Important methodological notes

- ROS is applied on observed censored data (no log transform unless explicitly stated)
- MLE assumes Gaussian censored likelihood
- KM is used for survival-based estimation
- Substitution methods become unstable at high censoring levels (≥85%)
- Likelihood-based CUSUM does not require imputation

---

## Outputs generated

- Bias tables
- RMSE tables
- ARL summary tables
- Plots:
  - Bias vs censoring
  - RMSE vs censoring
  - ARL vs shift

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

- Censored data inference
- Robust estimation under LOD constraints
- SPC monitoring under incomplete data
- Comparison of classical vs likelihood-based approaches
- Performance evaluation of CUSUM under censoring

---

## Author note

This framework is intended for academic research and thesis reproducibility. All parameters (censoring, shifts, sample sizes) can be modified for extended analysis.
