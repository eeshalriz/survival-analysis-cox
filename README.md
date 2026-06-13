# Survival Analysis: Kaplan-Meier & Cox Proportional Hazards

Analysing patient recovery times under two therapy types using non-parametric and semi-parametric survival models in R.

## Overview

This project applies survival analysis techniques to a clinical dataset comparing unmedicated and medicated therapy recovery times. It covers the full pipeline from Kaplan-Meier estimation through Cox regression and hypothesis testing.

## Methods

| Method | Purpose |
|--------|---------|
| Kaplan-Meier estimator | Non-parametric survival curve estimation |
| Survival probability at t=2 | Recovery probability by therapy group |
| Confidence intervals | Uncertainty quantification via KM bounds |
| Bayes' theorem | Conditional probability of therapy type given recovery |
| Cox proportional hazards | Semi-parametric regression for β estimation |
| Wald / LRT / Score tests | Three-way hypothesis testing for β = 0 |
| Parametric hazard model | Survival probability under assumed h₀(t) = 0.27 + 0.04t |

## Key Results

- **Medicated therapy** showed higher recovery probability within 2 years
- **Cox model β**: estimated from data; hazard ratio exp(β) quantifies the treatment effect
- All three hypothesis tests (Wald, LRT, Score) consistently assessed whether β differs significantly from 0
- Parametric survival estimates computed under a linear baseline hazard assumption

## Tech Stack

- **Language:** R
- **Package:** `survival`

## How to Run

```r
install.packages("survival")
source("survival_analysis.R")
```

The script prints a full summary to console and generates a KM survival curve plot.

## Files

```
survival_analysis.R   # Full annotated R script
README.md
```

## Author

**Eeshal Rizwan**  
BSc Statistical Data Science — Heriot-Watt University Dubai  
[linkedin.com/in/eeshalrizwan](https://linkedin.com/in/eeshalrizwan)
