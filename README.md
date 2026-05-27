# Censored Data Simulation for SPC and CUSUM Performance Evaluation

This repository contains an R-based simulation framework for evaluating statistical methods for censored data in the context of Statistical Process Control (SPC), with a focus on Phase I estimation and Phase II CUSUM monitoring.

## Overview

The simulation study compares different approaches for handling left-censored data, commonly encountered in food safety and environmental microbiology datasets where values fall below a limit of detection (LOD).

The main objective is to assess how different estimation and imputation strategies affect:

- Phase I estimation of process parameters (mean and standard deviation)
- Phase II CUSUM chart performance
- Detection performance under varying censoring levels and process shifts

## Methods Compared

The following approaches are evaluated:

### Phase I estimation
- Maximum Likelihood Estimation (MLE) via censored Gaussian model (`survreg`)
- Regression on Order Statistics (ROS) using the `NADA` package

### Phase II imputation / monitoring methods
- LOD substitution
- LOD / √2 substitution
- Stochastic Random Imputation (truncated normal approach)
- Likelihood-based CUSUM (no imputation required)

## Simulation Design

- Data are generated from a normal distribution with left censoring at a fixed LOD
- Censoring levels range from very low (≈0%) to 95%
- Sample size per simulation:
  - Phase I: 200 observations
  - Phase II (in-control): 25 observations
  - Phase II (post-shift): 3000 observations
- Process shifts:
  - 0σ (in-control)
  - 0.5σ shift
  - 1σ shift
- Number of simulations: 200 per scenario

The mean is adjusted according to censoring proportion to maintain consistent censoring levels across scenarios.

## CUSUM Implementation

- Standard CUSUM is implemented using the `qcr` package
- Likelihood-based CUSUM is implemented using log-likelihood ratios for censored and uncensored observations
- Performance is evaluated using:
  - Average Run Length (ARL)
  - Median run length
  - Standard deviation of run length (SDRL)
  - Quantiles (10%, 25%, 75%, 90%)

## Output

The script produces:

- Simulation results for all combinations of:
  - censoring level
  - method
  - process shift
- Summary performance metrics (ARL, quantiles, variability measures)
- Visualization of ARL performance across methods and censoring levels

## Visualization

The final plot shows:

- ARL as a function of process shift (Δ)
- Comparison of methods across censoring levels
- Uncertainty bands (interquartile range)

## Key Assumptions

- Data are normally distributed on the transformed scale
- Censoring is left-censoring at a fixed LOD
- Phase I and Phase II processes share the same underlying distribution except for introduced shifts

## Dependencies

- `NADA`
- `survival`
- `qcr`
- `ggplot2`
- `dplyr`
- `sn`

## Purpose

This code is intended for methodological research on:
- Handling censored microbiological data
- Robust estimation of process parameters
- SPC performance under censoring
- Comparison of classical and likelihood-based monitoring methods
