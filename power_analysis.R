################################################################################
# POWER ANALYSIS: A Priori Sample Size Determination
# Reproduces the manuscript's required N = 133
################################################################################

library(pwr)

cat("================================================================================\n")
cat("POWER ANALYSIS\n")
cat("================================================================================\n\n")

# ============================================================================
# ANALYTICAL POWER ANALYSIS (Cohen, 1988)
# ============================================================================
# Critical test: Model 3 from Script 02 (H1 Comparative Regression)
#   Meaning ~ Autonomous + GSC + Depression + Anxiety + Age + Sex
#   6 predictors total; testing single-predictor effect of Autonomous
#
# Parameters:
#   u  = 1        (numerator df: testing 1 predictor)
#   f2 = 0.0625   (medium effect: beta = .25, so f2 = .25^2 / (1 - .25^2) ~ .0667,
#                   or equivalently R2-change = .0588; Cohen's medium f2 = .0625)
#   sig.level = 0.05
#   power = 0.80

cat("Analytical Power Analysis (pwr::pwr.f2.test)\n")
cat("----------------------------------------------\n")
cat("  Effect size: f2 = 0.0625 (medium, corresponding to beta ~ .25)\n")
cat("  Numerator df (u): 1 (single predictor test)\n")
cat("  Significance level: 0.05\n")
cat("  Target power: 0.80\n\n")

result <- pwr.f2.test(u = 1, f2 = 0.0625, sig.level = 0.05, power = 0.80)
print(result)

v <- ceiling(result$v)  # residual df (denominator df)
k <- 6                  # number of predictors in Model 3
N_required <- v + k + 1

cat("\n  Residual df (v):", v, "\n")
cat("  Number of predictors (k):", k, "\n")
cat("  Required N = v + k + 1 =", N_required, "\n\n")

# ============================================================================
# MONTE CARLO SIMULATION (Confirmatory)
# ============================================================================
cat("Monte Carlo Simulation (10,000 iterations)\n")
cat("--------------------------------------------\n")
cat("  Simulating at N = 133 to verify power ~ 0.80\n\n")

set.seed(42)
n_sim <- 10000
n_target <- 133
beta_target <- 0.25

sig_count <- 0

for (i in 1:n_sim) {
  # Generate 6 predictors (uncorrelated for simplicity)
  X <- matrix(rnorm(n_target * 6), ncol = 6)
  colnames(X) <- c("Autonomous", "GSC", "Depression", "Anxiety", "Age", "Sex")

  # Generate outcome: only Autonomous has true effect (beta = 0.25)
  y <- beta_target * X[, 1] + rnorm(n_target)

  # Fit full model and test Autonomous
  dat <- data.frame(y = y, X)
  fit <- lm(y ~ Autonomous + GSC + Depression + Anxiety + Age + Sex, data = dat)
  p_auto <- summary(fit)$coefficients["Autonomous", "Pr(>|t|)"]

  if (p_auto < 0.05) sig_count <- sig_count + 1
}

simulated_power <- sig_count / n_sim

cat("  Simulated power at N =", n_target, ":", round(simulated_power, 3), "\n")
cat("  (", sig_count, "out of", n_sim, "simulations detected the effect)\n\n")

# ============================================================================
# SUMMARY
# ============================================================================
cat("================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

cat("Required sample size: N =", N_required, "\n")
cat("Actual sample size:   N = 429 (after listwise deletion from 514)\n")
cat("Adequacy ratio:       ", round(429 / N_required, 1), "x the required N\n\n")

cat("Manuscript wording:\n")
cat("  \"An a priori power analysis (Cohen, 1988; implemented via the pwr\n")
cat("  package in R) indicated that a sample size of N = 133 is required\n")
cat("  to achieve 80% power to detect a single-predictor effect of\n")
cat("  f2 = 0.0625 (corresponding to beta = .25) in a six-predictor\n")
cat("  regression model at alpha = .05.\"\n\n")

cat("================================================================================\n")
cat("POWER ANALYSIS COMPLETE\n")
cat("================================================================================\n")
