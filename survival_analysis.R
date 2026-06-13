# =============================================================================
# Survival Analysis: Kaplan-Meier Estimation & Cox Proportional Hazards Model
# Author:  Eeshal Rizwan
# Course:  Statistical Modelling — Heriot-Watt University
# =============================================================================
# This script analyses patient recovery times under two therapy types:
#   - Unmedicated therapy (control group)
#   - Medicated therapy (treatment group)
#
# Methods covered:
#   1. Kaplan-Meier (KM) survival estimation
#   2. Survival probability estimation at t = 2 years
#   3. Confidence intervals for survival estimates
#   4. Bayes' theorem applied to recovery probabilities
#   5. Cox proportional hazards model — MLE of regression coefficient β
#   6. Hypothesis testing: Wald test, Likelihood Ratio Test, Score test
#   7. Parametric hazard function — survival probability calculation
# =============================================================================

library(survival)

# -----------------------------------------------------------------------------
# Data generation (reproducible synthetic dataset)
# -----------------------------------------------------------------------------
# Simulates 30 patients split equally between two therapy types.
# Time: gamma-distributed recovery times
# Censor: binary censoring indicator (1 = observed, 0 = censored)
# Med: 0 = unmedicated therapy, 1 = medicated therapy

set.seed(42)
n_group <- 15

TimeA   <- rgamma(n_group, shape = 2, rate = 1.0)   # unmedicated
TimeB   <- rgamma(n_group, shape = 3, rate = 1.7)   # medicated
CensorA <- rbinom(n_group, 1, 0.6)
CensorB <- rbinom(n_group, 1, 0.8)

data <- data.frame(
  Time   = c(TimeA, TimeB),
  Censor = c(CensorA, CensorB),
  Med    = c(rep(0, n_group), rep(1, n_group))  # 0 = unmedicated, 1 = medicated
)

cat("=== Dataset summary ===\n")
cat(sprintf("  Total patients:   %d\n", nrow(data)))
cat(sprintf("  Unmedicated:      %d\n", sum(data$Med == 0)))
cat(sprintf("  Medicated:        %d\n", sum(data$Med == 1)))
cat(sprintf("  Censoring rate:   %.1f%%\n\n", 100 * mean(data$Censor == 0)))

# -----------------------------------------------------------------------------
# 1. Kaplan-Meier survival estimation
# -----------------------------------------------------------------------------
km_unmed <- survfit(Surv(Time, Censor) ~ 1, data = subset(data, Med == 0))
km_med   <- survfit(Surv(Time, Censor) ~ 1, data = subset(data, Med == 1))
km_all   <- survfit(Surv(Time, Censor) ~ 1, data = data)

# Plot KM curves
plot(km_unmed,
     col  = "steelblue",
     lwd  = 2,
     xlab = "Time (years)",
     ylab = "Survival probability S(t)",
     main = "Kaplan-Meier Survival Curves by Therapy Type",
     conf.int = TRUE)

lines(km_med, col = "darkorange", lwd = 2, conf.int = TRUE)
legend("topright",
       legend = c("Unmedicated therapy", "Medicated therapy"),
       col    = c("steelblue", "darkorange"),
       lwd    = 2, bty = "n")

# -----------------------------------------------------------------------------
# 2. Survival and recovery probabilities at t = 2 years
# -----------------------------------------------------------------------------
t_eval <- 2   # evaluation time point

S_unmed <- summary(km_unmed, times = t_eval)$surv
S_med   <- summary(km_med,   times = t_eval)$surv
S_all   <- summary(km_all,   times = t_eval)$surv

# Recovery probability = 1 - S(t)
P_recovery_unmed <- 1 - S_unmed
P_recovery_med   <- 1 - S_med
P_recovery_all   <- 1 - S_all

cat("=== Recovery probabilities at t = 2 years ===\n")
cat(sprintf("  Unmedicated therapy: %.4f\n", P_recovery_unmed))
cat(sprintf("  Medicated therapy:   %.4f\n", P_recovery_med))
cat(sprintf("  Overall:             %.4f\n\n", P_recovery_all))

# -----------------------------------------------------------------------------
# 3. Confidence interval for overall recovery probability
# -----------------------------------------------------------------------------
km_summary    <- summary(km_all, times = t_eval)
CI_lower      <- 1 - km_summary$upper   # note: survival CI is inverted for recovery
CI_upper      <- 1 - km_summary$lower

cat("=== 95% CI for overall recovery probability ===\n")
cat(sprintf("  Lower bound: %.4f\n", CI_lower))
cat(sprintf("  Upper bound: %.4f\n\n", CI_upper))

# -----------------------------------------------------------------------------
# 4. Bayes' theorem — probability a recovered patient was unmedicated
# -----------------------------------------------------------------------------
P_unmed       <- mean(data$Med == 0)
P_med_group   <- mean(data$Med == 1)

# Total probability of recovery (law of total probability)
P_recovery_total <- P_recovery_unmed * P_unmed + P_recovery_med * P_med_group

# P(unmedicated | recovered) via Bayes' theorem
P_unmed_given_recovered <- (P_recovery_unmed * P_unmed) / P_recovery_total

cat("=== Bayes' Theorem ===\n")
cat(sprintf("  P(recovered | unmedicated): %.4f\n", P_recovery_unmed))
cat(sprintf("  P(recovered | medicated):   %.4f\n", P_recovery_med))
cat(sprintf("  P(unmedicated | recovered): %.4f\n\n", P_unmed_given_recovered))

# -----------------------------------------------------------------------------
# 5. Cox Proportional Hazards Model
# -----------------------------------------------------------------------------
# Model: h(t | Med) = h_0(t) * exp(β * Med)
# β > 0 → medicated group has higher hazard (faster recovery if time = failure)

cox_fit <- coxph(Surv(Time, Censor) ~ Med, data = data)
cox_summary <- summary(cox_fit)

beta_hat <- cox_fit$coefficients["Med"]

cat("=== Cox Proportional Hazards Model ===\n")
cat(sprintf("  β (MLE):          %.4f\n", beta_hat))
cat(sprintf("  exp(β) [HR]:      %.4f\n", exp(beta_hat)))
cat(sprintf("  Interpretation:   Medicated group has %.1f× the hazard of unmedicated\n\n",
            exp(beta_hat)))

# -----------------------------------------------------------------------------
# 6. Hypothesis tests for β = 0
# -----------------------------------------------------------------------------
# H_0: β = 0 (medication has no effect)
# H_1: β ≠ 0 (medication affects recovery rate)

# Wald test
wald_stat   <- cox_summary$waldtest["test"]
wald_pvalue <- cox_summary$waldtest["pvalue"]

# Likelihood Ratio Test
lrt_results  <- anova(cox_fit)
lrt_stat     <- lrt_results$Chisq[2]
lrt_pvalue   <- lrt_results$`Pr(>|Chi|)`[2]

# Score test
score_stat   <- cox_summary$sctest["test"]
score_pvalue <- cox_summary$sctest["pvalue"]

cat("=== Hypothesis Tests (H_0: β = 0) ===\n")
cat(sprintf("  Wald test:   statistic = %.4f, p-value = %.4f\n", wald_stat,  wald_pvalue))
cat(sprintf("  LRT:         statistic = %.4f, p-value = %.4f\n", lrt_stat,   lrt_pvalue))
cat(sprintf("  Score test:  statistic = %.4f, p-value = %.4f\n\n", score_stat, score_pvalue))

alpha <- 0.05
if (wald_pvalue < alpha) {
  cat("  → Reject H_0: β is significantly different from 0 at the 5% level.\n\n")
} else {
  cat("  → Fail to reject H_0: insufficient evidence that β ≠ 0 at the 5% level.\n\n")
}

# -----------------------------------------------------------------------------
# 7. Parametric survival probabilities using a given baseline hazard
# -----------------------------------------------------------------------------
# Baseline hazard: h_0(t) = 0.27 + 0.04t (assumed parametric form)
# Cumulative hazard: H_0(t) = integral of h_0(t) from 0 to t
# Survival:          S(t | Med) = exp(-H_0(t) * exp(β * Med))

baseline_hazard  <- function(t) 0.27 + 0.04 * t

H0_2 <- integrate(baseline_hazard, 0, 2)$value   # baseline cumulative hazard to t=2

S_parametric_unmed <- exp(-H0_2)                          # Med = 0
S_parametric_med   <- exp(-H0_2 * exp(beta_hat))          # Med = 1

cat("=== Parametric Survival Probabilities (h_0(t) = 0.27 + 0.04t) ===\n")
cat(sprintf("  P(recovery by t=2 | unmedicated): %.4f\n", 1 - S_parametric_unmed))
cat(sprintf("  P(recovery by t=2 | medicated):   %.4f\n", 1 - S_parametric_med))
