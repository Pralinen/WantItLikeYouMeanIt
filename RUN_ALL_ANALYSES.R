################################################################################
# MASTER SCRIPT: Run ALL 10 Analyses on ALL 3 Datasets (IMPROVED VERSION)
# Author: Liam Le Guellaff Pallin
# Date: 2025-11-20
#
# Purpose:
#   Runs ALL 10 scripts on:
#   - Study 1 Only (N = 168, CSUS)
#   - Study 2 Only (N = 269, WSU)
#   - Combined (N = 437, Both)
#
# Improvements:
#   - Error handling: keeps running if one analysis fails
#   - Progress tracking: saves progress after each script
#   - All output in one organized folder: ALL_STUDIES_ANALYSIS/
#   - Detailed error log
#   - Time tracking per script
#
# ESTIMATED TIME: 6-8 hours (mostly Bayesian models)
#
# USAGE:
#   cd "/Users/pralinen/Desktop/HT25/PSYK12/DK3 - Kandidatarbete/DataAnalysis"
#   nohup Rscript ALL_STUDIES_ANALYSIS/RUN_ALL_ANALYSES.R > ALL_STUDIES_ANALYSIS/output.log 2>&1 &
#
# MONITOR PROGRESS:
#   tail -f ALL_STUDIES_ANALYSIS/output.log
#   cat ALL_STUDIES_ANALYSIS/PROGRESS.txt
#   cat ALL_STUDIES_ANALYSIS/ERROR_LOG.txt
#
################################################################################

library(haven)
library(brms)
library(lavaan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ppcor)  # For semi-partial correlations in Script 02

# Set working directory
setwd("/Users/pralinen/Desktop/HT25/PSYK12/DK3 - Kandidatarbete/DataAnalysis")

# Create organized folder structure
base_dir <- "ALL_STUDIES_ANALYSIS"
dir.create(paste0(base_dir, "/Results"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(base_dir, "/Models"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(base_dir, "/Figures"), showWarnings = FALSE, recursive = TRUE)

# Initialize tracking
error_log <- character(0)
progress_log <- character(0)
script_times <- list()

cat("################################################################\n")
cat("################################################################\n")
cat("RUNNING ALL 10 SCRIPTS ON ALL 3 DATASETS\n")
cat("################################################################\n")
cat("################################################################\n\n")

cat("Output folder:", base_dir, "\n")
cat("Estimated time: 6-8 hours\n")
cat("Started:", as.character(Sys.time()), "\n\n")

overall_start <- Sys.time()

# Helper function to log progress
log_progress <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste0("[", timestamp, "] ", msg)
  progress_log <<- c(progress_log, log_entry)
  cat(log_entry, "\n")
  writeLines(progress_log, paste0(base_dir, "/PROGRESS.txt"))
}

# Helper function to log errors
log_error <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  error_entry <- paste0("[", timestamp, "] ❌ ERROR: ", msg)
  error_log <<- c(error_log, error_entry)
  cat(error_entry, "\n")
  writeLines(error_log, paste0(base_dir, "/ERROR_LOG.txt"))
}

# Safe wrapper function
safe_run <- function(func, dataset_info, script_name) {
  script_start <- Sys.time()

  result <- tryCatch({
    log_progress(paste0("Starting ", script_name, " for ", dataset_info$label))
    output <- func(dataset_info)

    script_end <- Sys.time()
    duration <- difftime(script_end, script_start, units = "mins")

    log_progress(paste0("✅ Completed ", script_name, " for ", dataset_info$label,
                       " (", round(duration, 1), " min)"))

    script_times[[paste0(script_name, "_", dataset_info$name)]] <<- as.numeric(duration)

    list(success = TRUE, result = output, duration = duration)

  }, error = function(e) {
    script_end <- Sys.time()
    duration <- difftime(script_end, script_start, units = "mins")

    error_msg <- paste0(script_name, " for ", dataset_info$label, ": ", e$message)
    log_error(error_msg)

    list(success = FALSE, error = e$message, duration = duration)
  })

  return(result)
}

################################################################################
# DEFINE DATASETS
################################################################################

datasets <- list(
  list(name = "Study1", path = "Data/Study1_Only.sav", label = "Study 1 (N=168, CSUS)"),
  list(name = "Study2", path = "Data/Study2_Only.sav", label = "Study 2 (N=269, WSU)"),
  list(name = "Combined", path = "Data/Combined_Study1_Study2.sav", label = "Combined (N=437)")
)

log_progress("Datasets defined: Study 1, Study 2, Combined")

################################################################################
# SCRIPT 02: H1 Comparative Regression - PUBLICATION-READY VERSION
################################################################################

run_script_02 <- function(dataset_info) {
  cat("\n--- Script 02: H1 Comparative Regression Analysis ---\n")
  cat("Dataset:", dataset_info$label, "\n\n")

  # Load and prepare data
  data <- read_sav(dataset_info$path) %>%
    mutate(Meaning = MILjudgements, Depression = DEP_Corrected, Anxiety = GADscore) %>%
    filter(complete.cases(Autonomous, RAI, Meaning, Depression, Anxiety, Age, Sex))

  cat("Sample size: N =", nrow(data), "\n\n")

  # ============================================================================
  # MODEL 0: CONTROLS ONLY (BASELINE)
  # ============================================================================
  cat("MODEL 0: Controls Only (Baseline)\n")
  m_controls <- lm(Meaning ~ Depression + Anxiety + Age + Sex, data = data)
  r2_controls <- summary(m_controls)$r.squared
  cat("  R² =", round(r2_controls, 4), "\n\n")

  # ============================================================================
  # MODEL 1: AUTONOMOUS + CONTROLS
  # ============================================================================
  cat("MODEL 1: Autonomous + Controls\n")
  m1 <- lm(Meaning ~ Autonomous + Depression + Anxiety + Age + Sex, data = data)
  sum1 <- summary(m1)
  r2_auto <- sum1$r.squared
  delta_r2_auto <- r2_auto - r2_controls

  cat("  R² =", round(r2_auto, 4), "\n")
  cat("  ΔR² from controls =", round(delta_r2_auto, 4),
      "(", round(100 * delta_r2_auto, 2), "%)\n")

  # F-change test
  f_change_auto <- anova(m_controls, m1)
  cat("  F-change =", round(f_change_auto$F[2], 2),
      ", p =", format.pval(f_change_auto$`Pr(>F)`[2], digits = 3), "\n")

  auto_b <- coef(m1)["Autonomous"]
  auto_se <- sum1$coefficients["Autonomous", "Std. Error"]
  auto_t <- sum1$coefficients["Autonomous", "t value"]
  auto_p <- sum1$coefficients["Autonomous", "Pr(>|t|)"]
  cat("  Autonomous: b =", round(auto_b, 3), ", SE =", round(auto_se, 3),
      ", t =", round(auto_t, 2), ", p =", format.pval(auto_p, digits = 3), "\n\n")

  # ============================================================================
  # MODEL 2: RAI + CONTROLS
  # ============================================================================
  cat("MODEL 2: RAI + Controls\n")
  m2 <- lm(Meaning ~ RAI + Depression + Anxiety + Age + Sex, data = data)
  sum2 <- summary(m2)
  r2_rai <- sum2$r.squared
  delta_r2_rai <- r2_rai - r2_controls

  cat("  R² =", round(r2_rai, 4), "\n")
  cat("  ΔR² from controls =", round(delta_r2_rai, 4),
      "(", round(100 * delta_r2_rai, 2), "%)\n")

  # F-change test
  f_change_rai <- anova(m_controls, m2)
  cat("  F-change =", round(f_change_rai$F[2], 2),
      ", p =", format.pval(f_change_rai$`Pr(>F)`[2], digits = 3), "\n")

  rai_b <- coef(m2)["RAI"]
  rai_se <- sum2$coefficients["RAI", "Std. Error"]
  rai_t <- sum2$coefficients["RAI", "t value"]
  rai_p <- sum2$coefficients["RAI", "Pr(>|t|)"]
  cat("  RAI: b =", round(rai_b, 3), ", SE =", round(rai_se, 3),
      ", t =", round(rai_t, 2), ", p =", format.pval(rai_p, digits = 3), "\n\n")

  # ============================================================================
  # MODEL 3: BOTH PREDICTORS (CRITICAL TEST)
  # ============================================================================
  cat("MODEL 3: Both Autonomous + RAI + Controls\n")
  m3 <- lm(Meaning ~ Autonomous + RAI + Depression + Anxiety + Age + Sex, data = data)
  sum3 <- summary(m3)
  r2_both <- sum3$r.squared

  cat("  R² =", round(r2_both, 4), "\n")

  auto_b3 <- coef(m3)["Autonomous"]
  auto_se3 <- sum3$coefficients["Autonomous", "Std. Error"]
  auto_t3 <- sum3$coefficients["Autonomous", "t value"]
  auto_p3 <- sum3$coefficients["Autonomous", "Pr(>|t|)"]
  cat("  Autonomous: b =", round(auto_b3, 3), ", SE =", round(auto_se3, 3),
      ", t =", round(auto_t3, 2), ", p =", format.pval(auto_p3, digits = 3))
  if (auto_p3 < 0.05) cat(" ✓ Significant\n") else cat(" ~ Not significant\n")

  rai_b3 <- coef(m3)["RAI"]
  rai_se3 <- sum3$coefficients["RAI", "Std. Error"]
  rai_t3 <- sum3$coefficients["RAI", "t value"]
  rai_p3 <- sum3$coefficients["RAI", "Pr(>|t|)"]
  cat("  RAI: b =", round(rai_b3, 3), ", SE =", round(rai_se3, 3),
      ", t =", round(rai_t3, 2), ", p =", format.pval(rai_p3, digits = 3))
  if (rai_p3 >= 0.05) cat(" ✓ Not significant (H1 prediction)\n") else cat(" ⚠ Significant\n")
  cat("\n")

  # ============================================================================
  # INCREMENTAL VARIANCE TESTS (H1 CRITICAL)
  # ============================================================================
  cat("INCREMENTAL VARIANCE ANALYSIS:\n")

  # Test 1: Does RAI add beyond Autonomous?
  delta_r2_rai_beyond_auto <- r2_both - r2_auto
  f_test_rai_beyond <- anova(m1, m3)
  cat("  1. RAI beyond Autonomous:\n")
  cat("     ΔR² =", round(delta_r2_rai_beyond_auto, 4),
      "(", round(100 * delta_r2_rai_beyond_auto, 2), "%)\n")
  cat("     F =", round(f_test_rai_beyond$F[2], 2),
      ", p =", format.pval(f_test_rai_beyond$`Pr(>F)`[2], digits = 3))
  if (f_test_rai_beyond$`Pr(>F)`[2] >= 0.05) {
    cat(" ✓ RAI adds NO variance (H1c supported)\n\n")
  } else {
    cat(" ⚠ RAI adds variance (challenges H1c)\n\n")
  }

  # Test 2: Does Autonomous add beyond RAI?
  delta_r2_auto_beyond_rai <- r2_both - r2_rai
  f_test_auto_beyond <- anova(m2, m3)
  cat("  2. Autonomous beyond RAI:\n")
  cat("     ΔR² =", round(delta_r2_auto_beyond_rai, 4),
      "(", round(100 * delta_r2_auto_beyond_rai, 2), "%)\n")
  cat("     F =", round(f_test_auto_beyond$F[2], 2),
      ", p =", format.pval(f_test_auto_beyond$`Pr(>F)`[2], digits = 3))
  if (f_test_auto_beyond$`Pr(>F)`[2] < 0.05) {
    cat(" ✓ Autonomous adds variance (H1 supported)\n\n")
  } else {
    cat(" ~ Autonomous does not add variance\n\n")
  }

  # ============================================================================
  # SEMI-PARTIAL CORRELATIONS (UNIQUE VARIANCE)
  # ============================================================================
  cat("SEMI-PARTIAL CORRELATIONS (Unique Variance):\n")

  # Autonomous unique variance
  spcor_auto <- spcor.test(
    x = data$Autonomous,
    y = data$Meaning,
    z = data[, c("RAI", "Depression", "Anxiety", "Age", "Sex")]
  )
  sr2_auto <- spcor_auto$estimate^2
  p_auto_sp <- spcor_auto$p.value

  cat("  Autonomous (controlling RAI + covariates):\n")
  cat("    sr² =", round(sr2_auto, 4),
      "(", round(100 * sr2_auto, 2), "% unique variance)\n")
  cat("    p =", format.pval(p_auto_sp, digits = 3), "\n\n")

  # RAI unique variance
  spcor_rai <- spcor.test(
    x = data$RAI,
    y = data$Meaning,
    z = data[, c("Autonomous", "Depression", "Anxiety", "Age", "Sex")]
  )
  sr2_rai <- spcor_rai$estimate^2
  p_rai_sp <- spcor_rai$p.value

  cat("  RAI (controlling Autonomous + covariates):\n")
  cat("    sr² =", round(sr2_rai, 4),
      "(", round(100 * sr2_rai, 2), "% unique variance)\n")
  cat("    p =", format.pval(p_rai_sp, digits = 3), "\n")

  if (sr2_rai > 0.0001) {
    ratio <- sr2_auto / sr2_rai
    cat("    → Autonomous explains", round(ratio, 1),
        "times more unique variance than RAI\n\n")
  } else {
    cat("    → RAI explains essentially 0% unique variance\n\n")
  }

  # ============================================================================
  # MULTICOLLINEARITY CHECK
  # ============================================================================
  cat("MULTICOLLINEARITY:\n")

  r2_auto_on_others <- summary(lm(
    Autonomous ~ RAI + Depression + Anxiety + Age + Sex, data = data
  ))$r.squared
  vif_auto <- 1 / (1 - r2_auto_on_others)

  r2_rai_on_others <- summary(lm(
    RAI ~ Autonomous + Depression + Anxiety + Age + Sex, data = data
  ))$r.squared
  vif_rai <- 1 / (1 - r2_rai_on_others)

  cat("  Autonomous VIF =", round(vif_auto, 2))
  if (vif_auto < 5) cat(" ✓\n") else if (vif_auto < 10) cat(" ~ Moderate\n") else cat(" ⚠ High\n")

  cat("  RAI VIF =", round(vif_rai, 2))
  if (vif_rai < 5) cat(" ✓\n\n") else if (vif_rai < 10) cat(" ~ Moderate\n\n") else cat(" ⚠ High\n\n")

  # ============================================================================
  # H1 DECISION
  # ============================================================================
  cat("H1 HYPOTHESIS TEST SUMMARY:\n")
  h1a <- auto_p3 < 0.05
  h1b <- rai_p3 >= 0.05
  h1c <- f_test_rai_beyond$`Pr(>F)`[2] >= 0.05

  cat("  H1a (Autonomous remains significant): ", if(h1a) "✓ SUPPORTED\n" else "✗ NOT SUPPORTED\n")
  cat("  H1b (RAI becomes non-significant): ", if(h1b) "✓ SUPPORTED\n" else "✗ NOT SUPPORTED\n")
  cat("  H1c (RAI adds no incremental variance): ", if(h1c) "✓ SUPPORTED\n" else "✗ NOT SUPPORTED\n")

  if (h1a && h1b && h1c) {
    cat("  → OVERALL H1: ✓ FULLY SUPPORTED\n\n")
  } else if (h1a && h1b) {
    cat("  → OVERALL H1: ~ PARTIALLY SUPPORTED (H1a + H1b)\n\n")
  } else {
    cat("  → OVERALL H1: ✗ NOT SUPPORTED\n\n")
  }

  # Save comprehensive results
  results <- list(
    dataset = dataset_info$name,
    n = nrow(data),
    models = list(
      controls = list(model = m_controls, summary = summary(m_controls), r2 = r2_controls),
      m1_autonomous = list(model = m1, summary = sum1, r2 = r2_auto),
      m2_rai = list(model = m2, summary = sum2, r2 = r2_rai),
      m3_both = list(model = m3, summary = sum3, r2 = r2_both)
    ),
    incremental_variance = list(
      delta_r2_auto = delta_r2_auto,
      delta_r2_rai = delta_r2_rai,
      delta_r2_rai_beyond_auto = delta_r2_rai_beyond_auto,
      delta_r2_auto_beyond_rai = delta_r2_auto_beyond_rai,
      f_test_rai_beyond = f_test_rai_beyond,
      f_test_auto_beyond = f_test_auto_beyond
    ),
    semi_partial_correlations = list(
      sr2_auto = sr2_auto,
      sr2_rai = sr2_rai,
      p_auto = p_auto_sp,
      p_rai = p_rai_sp
    ),
    multicollinearity = list(
      vif_auto = vif_auto,
      vif_rai = vif_rai
    ),
    h1_tests = list(
      h1a = h1a,
      h1b = h1b,
      h1c = h1c,
      overall_supported = h1a && h1b && h1c
    )
  )

  saveRDS(results, file = paste0(base_dir, "/Results/Script02_", dataset_info$name, ".rds"))
  cat("Results saved\n")

  return(results)
}

################################################################################
# SCRIPT 03: H1 Path Decomposition Analysis - RAI → Autonomous → Meaning
################################################################################
#
# This Path Decomposition Analysis uses structural equation modeling to isolate
# the predictive pathways. Since RAI is a composite index (Autonomous - Controlled),
# this analysis demonstrates that the association between RAI and Meaning operates
# primarily through the Autonomous component.
#
# The mediation framework allows us to:
# - Decompose RAI's total effect into direct (c') and indirect (a×b) pathways
# - Test whether Autonomous mediates the RAI → Meaning relationship
# - Isolate variance attributable to each component
################################################################################

run_script_03 <- function(dataset_info) {
  data <- read_sav(dataset_info$path) %>%
    mutate(Meaning = MILjudgements, Depression = DEP_Corrected, Anxiety = GADscore) %>%
    filter(complete.cases(RAI, Autonomous, Meaning, Depression, Anxiety, Age, Sex))

  # Path Decomposition Model: Tests WHY RAI drops out (because Autonomous mediates it)
  mediation_model <- '
    # Path a: RAI → Autonomous
    Autonomous ~ a*RAI + Depression + Anxiety + Age + Sex

    # Path b: Autonomous → Meaning (controlling RAI)
    Meaning ~ b*Autonomous + c_prime*RAI + Depression + Anxiety + Age + Sex

    # Indirect effect (a × b)
    indirect := a * b

    # Total effect
    total := c_prime + (a * b)
  '

  fit <- sem(mediation_model, data = data, se = 'bootstrap', bootstrap = 1000)

  results <- list(
    dataset = dataset_info$name, n = nrow(data),
    fit = fit, estimates = parameterEstimates(fit),
    type = "RAI_to_Autonomous_to_Meaning"
  )

  saveRDS(results, file = paste0(base_dir, "/Results/Script03_", dataset_info$name, ".rds"))
  return(results)
}

################################################################################
# SCRIPT 03b: Supplementary Mediation - Autonomous → Depression → Meaning
################################################################################

run_script_03b <- function(dataset_info) {
  data <- read_sav(dataset_info$path) %>%
    mutate(Meaning = MILjudgements, Depression = DEP_Corrected, Anxiety = GADscore) %>%
    filter(complete.cases(Autonomous, Meaning, Depression, Anxiety, Age, Sex))

  # Supplementary: Tests HOW Autonomous works (partly through reducing depression)
  mediation_model <- '
    # Path a: Autonomous → Depression
    Depression ~ a*Autonomous + Age + Sex

    # Path b: Depression → Meaning (controlling Autonomous)
    Meaning ~ b*Depression + c_prime*Autonomous + Anxiety + Age + Sex

    # Indirect effect
    indirect := a * b

    # Total effect
    total := c_prime + (a * b)
  '

  fit <- sem(mediation_model, data = data, se = 'bootstrap', bootstrap = 1000)

  results <- list(
    dataset = dataset_info$name, n = nrow(data),
    fit = fit, estimates = parameterEstimates(fit),
    type = "Autonomous_to_Depression_to_Meaning"
  )

  saveRDS(results, file = paste0(base_dir, "/Results/Script03b_", dataset_info$name, ".rds"))
  return(results)
}

################################################################################
# SCRIPT 04: H1 Bayesian Ordinal - PUBLICATION-READY WITH 3-MODEL COMPARISON
################################################################################

run_script_04 <- function(dataset_info) {
  cat("\n--- Script 04: H1 Bayesian Ordinal (3-Model DAG Comparison) ---\n")
  cat("Dataset:", dataset_info$label, "\n\n")

  # Load and prepare data
  data <- read_sav(dataset_info$path) %>%
    mutate(Meaning = MILjudgements, Depression = DEP_Corrected, Anxiety = GADscore, RAI) %>%
    filter(complete.cases(Autonomous, RAI, Meaning, Depression, Anxiety, Age, Sex)) %>%
    mutate(Meaning_ord = ordered(Meaning))

  cat("N =", nrow(data), "\n")
  cat("Meaning levels:", nlevels(data$Meaning_ord), "\n\n")

  priors <- c(
    prior(student_t(3, 0, 2.5), class = "b"),
    prior(student_t(3, 0, 2.5), class = "Intercept")
  )

  cat("Testing 3 competing ordinal models:\n")
  cat("  Model 1: Autonomous only + controls (DAG 2)\n")
  cat("  Model 2: RAI only + controls (DAG 1)\n")
  cat("  Model 3: Both predictors + controls (DAG 3)\n\n")

  # MODEL 1: AUTONOMOUS ONLY
  cat("MODEL 1: Autonomous Only\n")
  start1 <- Sys.time()
  fit1 <- brm(
    Meaning_ord ~ Autonomous + Depression + Anxiety + Age + Sex,
    data = data, family = cumulative("probit"), prior = priors,
    chains = 4, iter = 3000, warmup = 1500, seed = 42,
    backend = "cmdstanr", cores = 4,
    control = list(adapt_delta = 0.95),
    file = paste0(base_dir, "/Models/Script04_Model1_", dataset_info$name),
    silent = 2, refresh = 0
  )
  elapsed1 <- difftime(Sys.time(), start1, units = "mins")
  cat("  Complete in", round(elapsed1, 1), "min\n")

  rhats1 <- rhat(fit1)
  max_rhat1 <- max(rhats1, na.rm = TRUE)
  cat("  R-hat:", round(max_rhat1, 4))
  if (!is.na(max_rhat1) && max_rhat1 < 1.01) cat(" ✓\n") else cat(" ⚠\n")

  loo1 <- loo(fit1)
  cat("  LOO-IC:", round(loo1$estimates["looic", "Estimate"], 1), "\n\n")

  # MODEL 2: RAI ONLY
  cat("MODEL 2: RAI Only\n")
  start2 <- Sys.time()
  fit2 <- brm(
    Meaning_ord ~ RAI + Depression + Anxiety + Age + Sex,
    data = data, family = cumulative("probit"), prior = priors,
    chains = 4, iter = 3000, warmup = 1500, seed = 42,
    backend = "cmdstanr", cores = 4,
    control = list(adapt_delta = 0.95),
    file = paste0(base_dir, "/Models/Script04_Model2_", dataset_info$name),
    silent = 2, refresh = 0
  )
  elapsed2 <- difftime(Sys.time(), start2, units = "mins")
  cat("  Complete in", round(elapsed2, 1), "min\n")

  rhats2 <- rhat(fit2)
  max_rhat2 <- max(rhats2, na.rm = TRUE)
  cat("  R-hat:", round(max_rhat2, 4))
  if (!is.na(max_rhat2) && max_rhat2 < 1.01) cat(" ✓\n") else cat(" ⚠\n")

  loo2 <- loo(fit2)
  cat("  LOO-IC:", round(loo2$estimates["looic", "Estimate"], 1), "\n\n")

  # MODEL 3: BOTH PREDICTORS
  cat("MODEL 3: Both Autonomous + RAI\n")
  start3 <- Sys.time()
  fit3 <- brm(
    Meaning_ord ~ Autonomous + RAI + Depression + Anxiety + Age + Sex,
    data = data, family = cumulative("probit"), prior = priors,
    chains = 4, iter = 3000, warmup = 1500, seed = 42,
    backend = "cmdstanr", cores = 4,
    control = list(adapt_delta = 0.95),
    file = paste0(base_dir, "/Models/Script04_Model3_", dataset_info$name),
    silent = 2, refresh = 0
  )
  elapsed3 <- difftime(Sys.time(), start3, units = "mins")
  cat("  Complete in", round(elapsed3, 1), "min\n")

  rhats3 <- rhat(fit3)
  max_rhat3 <- max(rhats3, na.rm = TRUE)
  cat("  R-hat:", round(max_rhat3, 4))
  if (!is.na(max_rhat3) && max_rhat3 < 1.01) cat(" ✓\n") else cat(" ⚠\n")

  loo3 <- loo(fit3)
  cat("  LOO-IC:", round(loo3$estimates["looic", "Estimate"], 1), "\n\n")

  # MODEL COMPARISON
  cat("MODEL COMPARISON (LOO-IC):\n")
  loo_comp <- loo_compare(loo1, loo2, loo3)
  print(loo_comp)
  cat("\n")

  best_model <- rownames(loo_comp)[1]
  cat("Best model:", best_model, "\n")
  if (grepl("fit1", best_model)) {
    cat("  → Autonomous-only wins (supports DAG 2 / H1)\n\n")
  } else if (grepl("fit2", best_model)) {
    cat("  → RAI-only wins (supports DAG 1 / traditional SDT)\n\n")
  } else {
    cat("  → Both-predictors wins (supports DAG 3 / dual pathway)\n\n")
  }

  # Extract key parameters
  post1 <- posterior_summary(fit1)
  post3 <- posterior_summary(fit3)

  auto_est1 <- post1["b_Autonomous", "Estimate"]
  auto_lower1 <- post1["b_Autonomous", "Q2.5"]
  auto_upper1 <- post1["b_Autonomous", "Q97.5"]

  auto_est3 <- post3["b_Autonomous", "Estimate"]
  rai_est3 <- post3["b_RAI", "Estimate"]

  cat("KEY FINDINGS:\n")
  cat("  Model 1 - Autonomous:", round(auto_est1, 3),
      "[", round(auto_lower1, 3), ",", round(auto_upper1, 3), "]")
  if (auto_lower1 > 0) cat(" ✓ Credibly positive\n") else cat(" ~ Uncertain\n")

  cat("  Model 3 - Autonomous:", round(auto_est3, 3))
  if (abs(auto_est3 - auto_est1) / auto_est1 < 0.1) {
    cat(" (stable when RAI added)\n")
  } else {
    cat(" (changes when RAI added)\n")
  }

  cat("  Model 3 - RAI:", round(rai_est3, 3))
  rai_lower3 <- post3["b_RAI", "Q2.5"]
  rai_upper3 <- post3["b_RAI", "Q97.5"]
  if (rai_lower3 * rai_upper3 > 0) {
    cat(" (credible)\n")
  } else {
    cat(" (uncertain - crosses zero)\n")
  }
  cat("\n")

  # Posterior predictive checks
  cat("Generating posterior predictive checks...\n")
  pp1 <- pp_check(fit1, ndraws = 100, type = "bars")
  pp_file1 <- paste0(base_dir, "/Results/Script04_Model1_", dataset_info$name, "_pp_check.pdf")
  pdf(pp_file1, width = 10, height = 6)
  print(pp1 + ggtitle(paste("Model 1 (Autonomous):", dataset_info$label)) + theme_minimal())
  dev.off()

  pp3 <- pp_check(fit3, ndraws = 100, type = "bars")
  pp_file3 <- paste0(base_dir, "/Results/Script04_Model3_", dataset_info$name, "_pp_check.pdf")
  pdf(pp_file3, width = 10, height = 6)
  print(pp3 + ggtitle(paste("Model 3 (Both):", dataset_info$label)) + theme_minimal())
  dev.off()
  cat("PP checks saved\n")

  # Save comprehensive results
  results <- list(
    dataset = dataset_info$name,
    n = nrow(data),
    model1 = list(
      fit = fit1,
      post_summary = post1,
      loo = loo1,
      diagnostics = list(max_rhat = max_rhat1, elapsed_mins = as.numeric(elapsed1))
    ),
    model2 = list(
      fit = fit2,
      loo = loo2,
      diagnostics = list(max_rhat = max_rhat2, elapsed_mins = as.numeric(elapsed2))
    ),
    model3 = list(
      fit = fit3,
      post_summary = post3,
      loo = loo3,
      diagnostics = list(max_rhat = max_rhat3, elapsed_mins = as.numeric(elapsed3))
    ),
    comparison = list(
      loo_compare = loo_comp,
      best_model = best_model
    ),
    key_findings = list(
      model1_autonomous = c(est = auto_est1, lower = auto_lower1, upper = auto_upper1),
      model3_autonomous = auto_est3,
      model3_rai = rai_est3
    )
  )

  saveRDS(results, file = paste0(base_dir, "/Results/Script04_", dataset_info$name, ".rds"))
  cat("Results saved\n\n")

  return(results)
}

################################################################################
# SCRIPT 05: H2 Controlled Motivation - PUBLICATION-READY VERSION
################################################################################

run_script_05 <- function(dataset_info) {
  cat("\n--- Script 05: H2 Controlled Motivation Analysis ---\n")
  cat("Dataset:", dataset_info$label, "\n\n")

  # Load and prepare data
  data <- read_sav(dataset_info$path) %>%
    mutate(Meaning = MILjudgements, Depression = DEP_Corrected, Anxiety = GADscore) %>%
    filter(complete.cases(Autonomous, Controlled, Meaning, Depression, Anxiety, Age, Sex))

  cat("Sample size: N =", nrow(data), "\n\n")
  cat("H2: Testing if Controlled motivation is irrelevant once Autonomous is accounted for\n\n")

  # ============================================================================
  # MODEL 0: CONTROLS ONLY
  # ============================================================================
  m_controls <- lm(Meaning ~ Depression + Anxiety + Age + Sex, data = data)
  r2_controls <- summary(m_controls)$r.squared

  # ============================================================================
  # MODEL 1: CONTROLLED ALONE
  # ============================================================================
  cat("MODEL 1: Controlled + Controls\n")
  m_controlled <- lm(Meaning ~ Controlled + Depression + Anxiety + Age + Sex, data = data)
  sum_controlled <- summary(m_controlled)
  r2_controlled <- sum_controlled$r.squared
  delta_r2_controlled <- r2_controlled - r2_controls

  cat("  R² =", round(r2_controlled, 4), "\n")
  cat("  ΔR² from controls =", round(delta_r2_controlled, 4),
      "(", round(100 * delta_r2_controlled, 2), "%)\n")

  f_change_controlled <- anova(m_controls, m_controlled)
  cat("  F-change =", round(f_change_controlled$F[2], 2),
      ", p =", format.pval(f_change_controlled$`Pr(>F)`[2], digits = 3), "\n")

  ctrl_b <- coef(m_controlled)["Controlled"]
  ctrl_se <- sum_controlled$coefficients["Controlled", "Std. Error"]
  ctrl_t <- sum_controlled$coefficients["Controlled", "t value"]
  ctrl_p <- sum_controlled$coefficients["Controlled", "Pr(>|t|)"]
  cat("  Controlled: b =", round(ctrl_b, 3), ", SE =", round(ctrl_se, 3),
      ", t =", round(ctrl_t, 2), ", p =", format.pval(ctrl_p, digits = 3), "\n\n")

  # ============================================================================
  # MODEL 2: AUTONOMOUS ALONE
  # ============================================================================
  cat("MODEL 2: Autonomous + Controls\n")
  m_autonomous <- lm(Meaning ~ Autonomous + Depression + Anxiety + Age + Sex, data = data)
  sum_autonomous <- summary(m_autonomous)
  r2_autonomous <- sum_autonomous$r.squared
  delta_r2_autonomous <- r2_autonomous - r2_controls

  cat("  R² =", round(r2_autonomous, 4), "\n")
  cat("  ΔR² from controls =", round(delta_r2_autonomous, 4),
      "(", round(100 * delta_r2_autonomous, 2), "%)\n\n")

  # ============================================================================
  # MODEL 3: BOTH PREDICTORS (H2 CRITICAL TEST)
  # ============================================================================
  cat("MODEL 3: Both Autonomous + Controlled + Controls\n")
  m_both <- lm(Meaning ~ Autonomous + Controlled + Depression + Anxiety + Age + Sex, data = data)
  sum_both <- summary(m_both)
  r2_both <- sum_both$r.squared

  cat("  R² =", round(r2_both, 4), "\n")

  auto_b3 <- coef(m_both)["Autonomous"]
  auto_p3 <- sum_both$coefficients["Autonomous", "Pr(>|t|)"]
  cat("  Autonomous: b =", round(auto_b3, 3), ", p =", format.pval(auto_p3, digits = 3))
  if (auto_p3 < 0.05) cat(" ✓ Significant (H2a)\n") else cat(" ~ Not significant\n")

  ctrl_b3 <- coef(m_both)["Controlled"]
  ctrl_p3 <- sum_both$coefficients["Controlled", "Pr(>|t|)"]
  cat("  Controlled: b =", round(ctrl_b3, 3), ", p =", format.pval(ctrl_p3, digits = 3))
  if (ctrl_p3 >= 0.05) cat(" ✓ Not significant (H2b)\n") else cat(" ⚠ Significant\n")
  cat("\n")

  # ============================================================================
  # INCREMENTAL VARIANCE TESTS (H2c)
  # ============================================================================
  cat("INCREMENTAL VARIANCE ANALYSIS:\n")

  # Test: Does Controlled add beyond Autonomous?
  delta_r2_ctrl_beyond_auto <- r2_both - r2_autonomous
  f_test_ctrl_beyond <- anova(m_autonomous, m_both)

  cat("  Controlled beyond Autonomous:\n")
  cat("    ΔR² =", round(delta_r2_ctrl_beyond_auto, 4),
      "(", round(100 * delta_r2_ctrl_beyond_auto, 2), "%)\n")
  cat("    F =", round(f_test_ctrl_beyond$F[2], 2),
      ", p =", format.pval(f_test_ctrl_beyond$`Pr(>F)`[2], digits = 3))
  if (f_test_ctrl_beyond$`Pr(>F)`[2] >= 0.05) {
    cat(" ✓ Controlled adds NO variance (H2c supported)\n\n")
  } else {
    cat(" ⚠ Controlled adds variance (challenges H2c)\n\n")
  }

  # ============================================================================
  # SEMI-PARTIAL CORRELATIONS
  # ============================================================================
  cat("SEMI-PARTIAL CORRELATIONS:\n")

  # Autonomous unique variance
  spcor_auto <- spcor.test(
    x = data$Autonomous,
    y = data$Meaning,
    z = data[, c("Controlled", "Depression", "Anxiety", "Age", "Sex")]
  )
  sr2_auto <- spcor_auto$estimate^2

  cat("  Autonomous (controlling Controlled + covariates):\n")
  cat("    sr² =", round(sr2_auto, 4),
      "(", round(100 * sr2_auto, 2), "% unique variance)\n\n")

  # Controlled unique variance
  spcor_ctrl <- spcor.test(
    x = data$Controlled,
    y = data$Meaning,
    z = data[, c("Autonomous", "Depression", "Anxiety", "Age", "Sex")]
  )
  sr2_ctrl <- spcor_ctrl$estimate^2

  cat("  Controlled (controlling Autonomous + covariates):\n")
  cat("    sr² =", round(sr2_ctrl, 4),
      "(", round(100 * sr2_ctrl, 2), "% unique variance)\n\n")

  # ============================================================================
  # H2 DECISION
  # ============================================================================
  cat("H2 HYPOTHESIS TEST SUMMARY:\n")
  h2a <- auto_p3 < 0.05
  h2b <- ctrl_p3 >= 0.05
  h2c <- f_test_ctrl_beyond$`Pr(>F)`[2] >= 0.05

  cat("  H2a (Autonomous remains significant): ", if(h2a) "✓ SUPPORTED\n" else "✗ NOT SUPPORTED\n")
  cat("  H2b (Controlled becomes non-significant): ", if(h2b) "✓ SUPPORTED\n" else "✗ NOT SUPPORTED\n")
  cat("  H2c (Controlled adds no incremental variance): ", if(h2c) "✓ SUPPORTED\n" else "✗ NOT SUPPORTED\n")

  if (h2a && h2b && h2c) {
    cat("  → OVERALL H2: ✓ FULLY SUPPORTED\n\n")
  } else if (h2b && h2c) {
    cat("  → OVERALL H2: ~ PARTIALLY SUPPORTED (H2b + H2c)\n\n")
  } else {
    cat("  → OVERALL H2: ✗ NOT SUPPORTED\n\n")
  }

  # Save comprehensive results
  results <- list(
    dataset = dataset_info$name,
    n = nrow(data),
    models = list(
      controls = list(model = m_controls, r2 = r2_controls),
      m1_controlled = list(model = m_controlled, summary = sum_controlled, r2 = r2_controlled),
      m2_autonomous = list(model = m_autonomous, summary = sum_autonomous, r2 = r2_autonomous),
      m3_both = list(model = m_both, summary = sum_both, r2 = r2_both)
    ),
    incremental_variance = list(
      delta_r2_controlled = delta_r2_controlled,
      delta_r2_autonomous = delta_r2_autonomous,
      delta_r2_ctrl_beyond_auto = delta_r2_ctrl_beyond_auto,
      f_test_ctrl_beyond = f_test_ctrl_beyond
    ),
    semi_partial_correlations = list(
      sr2_auto = sr2_auto,
      sr2_ctrl = sr2_ctrl
    ),
    h2_tests = list(
      h2a = h2a,
      h2b = h2b,
      h2c = h2c,
      overall_supported = h2a && h2b && h2c
    )
  )

  saveRDS(results, file = paste0(base_dir, "/Results/Script05_", dataset_info$name, ".rds"))
  cat("Results saved\n")

  return(results)
}

################################################################################
# SCRIPT 06: H3 SEM DAG Comparison (Proper DAG modeling)
################################################################################

run_script_06 <- function(dataset_info) {
  data_raw <- read_sav(dataset_info$path) %>%
    mutate(Meaning = MILjudgements, Depression = DEP_Corrected, Anxiety = GADscore) %>%
    filter(complete.cases(Autonomous, Controlled, RAI, Meaning, Depression, Anxiety, Age, Sex))

  # Standardize all variables for proper DAG comparison
  data <- data_raw %>%
    mutate(
      Meaning_z = scale(Meaning)[,1],
      Autonomous_z = scale(Autonomous)[,1],
      Controlled_z = scale(Controlled)[,1],
      RAI_z = scale(RAI)[,1],
      Depression_z = scale(Depression)[,1],
      Anxiety_z = scale(Anxiety)[,1],
      Age_z = scale(Age)[,1],
      Sex_z = as.numeric(Sex) - mean(as.numeric(Sex))
    )

  # DAG 1: RAI is causally primary (traditional SDT)
  dag1 <- '
    # Regressions
    Meaning_z ~ b_rai*RAI_z + b_dep*Depression_z + b_anx*Anxiety_z + b_age*Age_z + b_sex*Sex_z
    RAI_z ~ Depression_z + Anxiety_z + Age_z + Sex_z

    # Covariances among exogenous variables
    Depression_z ~~ Anxiety_z + Age_z + Sex_z
    Anxiety_z ~~ Age_z + Sex_z
    Age_z ~~ Sex_z
  '

  # DAG 2: Autonomous primary, Controlled = 0 (YOUR HYPOTHESIS)
  dag2 <- '
    # Regressions
    Meaning_z ~ b_aut*Autonomous_z + 0*Controlled_z + b_dep*Depression_z + b_anx*Anxiety_z + b_age*Age_z + b_sex*Sex_z
    Autonomous_z ~ Depression_z + Anxiety_z + Age_z + Sex_z
    Controlled_z ~ Depression_z + Anxiety_z + Age_z + Sex_z
    RAI_z ~ Autonomous_z + Controlled_z

    # Covariances
    Depression_z ~~ Anxiety_z + Age_z + Sex_z
    Anxiety_z ~~ Age_z + Sex_z
    Age_z ~~ Sex_z
    Autonomous_z ~~ Controlled_z

    # Residual covariance
    Meaning_z ~~ RAI_z
  '

  # DAG 3: Dual-pathway (both Autonomous and Controlled free)
  dag3 <- '
    # Regressions
    Meaning_z ~ b_aut*Autonomous_z + b_ctrl*Controlled_z + b_dep*Depression_z + b_anx*Anxiety_z + b_age*Age_z + b_sex*Sex_z
    Autonomous_z ~ Depression_z + Anxiety_z + Age_z + Sex_z
    Controlled_z ~ Depression_z + Anxiety_z + Age_z + Sex_z
    RAI_z ~ Autonomous_z + Controlled_z

    # Covariances
    Depression_z ~~ Anxiety_z + Age_z + Sex_z
    Anxiety_z ~~ Age_z + Sex_z
    Age_z ~~ Sex_z
    Autonomous_z ~~ Controlled_z

    # Residual covariance
    Meaning_z ~~ RAI_z
  '

  # Fit models with proper estimator
  fit1 <- sem(dag1, data = data, estimator = "MLR", se = "robust")
  fit2 <- sem(dag2, data = data, estimator = "MLR", se = "robust")
  fit3 <- sem(dag3, data = data, estimator = "MLR", se = "robust")

  results <- list(
    dataset = dataset_info$name, n = nrow(data),
    dag1 = fit1, dag2 = fit2, dag3 = fit3,
    fit_indices = list(
      dag1 = fitMeasures(fit1, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic")),
      dag2 = fitMeasures(fit2, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic")),
      dag3 = fitMeasures(fit3, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic"))
    )
  )

  saveRDS(results, file = paste0(base_dir, "/Results/Script06_", dataset_info$name, ".rds"))
  return(results)
}

################################################################################
# SCRIPT 07: H4 Robustness Outliers - PUBLICATION-READY VERSION
################################################################################

run_script_07 <- function(dataset_info) {
  cat("\n--- Script 07: H4 Robustness Analysis (Influential Observations) ---\n")
  cat("Dataset:", dataset_info$label, "\n\n")

  # Load and prepare data
  data <- read_sav(dataset_info$path) %>%
    mutate(Meaning = MILjudgements, Depression = DEP_Corrected, Anxiety = GADscore) %>%
    filter(complete.cases(Autonomous, RAI, Meaning, Depression, Anxiety, Age, Sex))

  cat("Original sample size: N =", nrow(data), "\n\n")

  # ============================================================================
  # BASELINE MODEL (MODEL 3 FROM H1)
  # ============================================================================
  cat("BASELINE MODEL: Both Autonomous + RAI\n")
  m_baseline <- lm(Meaning ~ Autonomous + RAI + Depression + Anxiety + Age + Sex, data = data)
  sum_baseline <- summary(m_baseline)

  auto_b_base <- coef(m_baseline)["Autonomous"]
  auto_se_base <- sum_baseline$coefficients["Autonomous", "Std. Error"]
  auto_t_base <- sum_baseline$coefficients["Autonomous", "t value"]
  auto_p_base <- sum_baseline$coefficients["Autonomous", "Pr(>|t|)"]

  rai_b_base <- coef(m_baseline)["RAI"]
  rai_se_base <- sum_baseline$coefficients["RAI", "Std. Error"]
  rai_t_base <- sum_baseline$coefficients["RAI", "t value"]
  rai_p_base <- sum_baseline$coefficients["RAI", "Pr(>|t|)"]

  cat("  Autonomous: b =", round(auto_b_base, 3), ", SE =", round(auto_se_base, 3),
      ", t =", round(auto_t_base, 2), ", p =", format.pval(auto_p_base, digits = 3), "\n")
  cat("  RAI: b =", round(rai_b_base, 3), ", SE =", round(rai_se_base, 3),
      ", t =", round(rai_t_base, 2), ", p =", format.pval(rai_p_base, digits = 3), "\n")
  cat("  R² =", round(sum_baseline$r.squared, 4), "\n\n")

  # ============================================================================
  # IDENTIFY INFLUENTIAL CASES (COOK'S DISTANCE)
  # ============================================================================
  cat("INFLUENTIAL CASES (Cook's Distance):\n")
  cooks_d <- cooks.distance(m_baseline)
  threshold <- 4 / nrow(data)

  cat("  Threshold (4/n):", round(threshold, 4), "\n")

  influential_idx <- which(cooks_d > threshold)
  n_influential <- length(influential_idx)

  cat("  N influential cases:", n_influential,
      "(", round(100 * n_influential / nrow(data), 1), "%)\n\n")

  if (n_influential > 0) {
    # Create robust dataset
    data_robust <- data[-influential_idx, ]
    cat("Robust sample size: N =", nrow(data_robust),
        "(excluded", n_influential, "cases)\n\n")

    # Re-fit model without influential cases
    cat("ROBUST MODEL (Influential cases excluded):\n")
    m_robust <- lm(Meaning ~ Autonomous + RAI + Depression + Anxiety + Age + Sex, data = data_robust)
    sum_robust <- summary(m_robust)

    auto_b_robust <- coef(m_robust)["Autonomous"]
    auto_se_robust <- sum_robust$coefficients["Autonomous", "Std. Error"]
    auto_t_robust <- sum_robust$coefficients["Autonomous", "t value"]
    auto_p_robust <- sum_robust$coefficients["Autonomous", "Pr(>|t|)"]

    rai_b_robust <- coef(m_robust)["RAI"]
    rai_se_robust <- sum_robust$coefficients["RAI", "Std. Error"]
    rai_t_robust <- sum_robust$coefficients["RAI", "t value"]
    rai_p_robust <- sum_robust$coefficients["RAI", "Pr(>|t|)"]

    cat("  Autonomous: b =", round(auto_b_robust, 3), ", SE =", round(auto_se_robust, 3),
        ", t =", round(auto_t_robust, 2), ", p =", format.pval(auto_p_robust, digits = 3), "\n")
    cat("  RAI: b =", round(rai_b_robust, 3), ", SE =", round(rai_se_robust, 3),
        ", t =", round(rai_t_robust, 2), ", p =", format.pval(rai_p_robust, digits = 3), "\n")
    cat("  R² =", round(sum_robust$r.squared, 4), "\n\n")

    # ============================================================================
    # COMPARE BASELINE VS ROBUST
    # ============================================================================
    cat("COMPARISON (Baseline vs Robust):\n")

    # Autonomous comparison
    auto_b_change <- ((auto_b_robust - auto_b_base) / auto_b_base) * 100
    auto_se_change <- ((auto_se_robust - auto_se_base) / auto_se_base) * 100

    cat("  Autonomous:\n")
    cat("    Baseline: b =", round(auto_b_base, 3), ", SE =", round(auto_se_base, 3), "\n")
    cat("    Robust:   b =", round(auto_b_robust, 3), ", SE =", round(auto_se_robust, 3), "\n")
    cat("    Change: b", sprintf("%+.1f%%", auto_b_change), ", SE", sprintf("%+.1f%%", auto_se_change), "\n")

    if (auto_p_base < 0.05 && auto_p_robust < 0.05) {
      cat("    ✓ Significant in both (robust finding)\n")
      auto_robust <- TRUE
    } else if (auto_p_base >= 0.05 && auto_p_robust >= 0.05) {
      cat("    ✓ Non-significant in both (consistent)\n")
      auto_robust <- TRUE
    } else {
      cat("    ⚠ Significance differs (not robust)\n")
      auto_robust <- FALSE
    }
    cat("\n")

    # RAI comparison
    rai_b_change <- ((rai_b_robust - rai_b_base) / abs(rai_b_base)) * 100
    rai_se_change <- ((rai_se_robust - rai_se_base) / rai_se_base) * 100

    cat("  RAI:\n")
    cat("    Baseline: b =", round(rai_b_base, 3), ", SE =", round(rai_se_base, 3), "\n")
    cat("    Robust:   b =", round(rai_b_robust, 3), ", SE =", round(rai_se_robust, 3), "\n")
    cat("    Change: b", sprintf("%+.1f%%", rai_b_change), ", SE", sprintf("%+.1f%%", rai_se_change), "\n")

    if (rai_p_base >= 0.05 && rai_p_robust >= 0.05) {
      cat("    ✓ Non-significant in both (H1 robust)\n")
      rai_robust <- TRUE
    } else if (rai_p_base < 0.05 && rai_p_robust < 0.05) {
      cat("    ~ Significant in both\n")
      rai_robust <- TRUE
    } else {
      cat("    ⚠ Significance differs\n")
      rai_robust <- FALSE
    }
    cat("\n")

    # Overall robustness
    cat("H4 ROBUSTNESS CHECK:\n")
    if (auto_robust && rai_robust) {
      cat("  ✓ ROBUST: Results consistent after excluding influential cases\n\n")
    } else {
      cat("  ⚠ NOT FULLY ROBUST: Some changes in significance\n\n")
    }

  } else {
    cat("No influential cases detected - results robust by default\n\n")
    m_robust <- m_baseline
    sum_robust <- sum_baseline
    auto_b_robust <- auto_b_base
    auto_p_robust <- auto_p_base
    rai_b_robust <- rai_b_base
    rai_p_robust <- rai_p_base
    auto_robust <- TRUE
    rai_robust <- TRUE
    auto_b_change <- 0
    rai_b_change <- 0
  }

  # Save comprehensive results
  results <- list(
    dataset = dataset_info$name,
    n_original = nrow(data),
    n_robust = if(n_influential > 0) nrow(data_robust) else nrow(data),
    n_influential = n_influential,
    threshold = threshold,
    influential_idx = influential_idx,
    model_baseline = list(model = m_baseline, summary = sum_baseline),
    model_robust = list(model = m_robust, summary = sum_robust),
    comparison = list(
      auto_baseline = c(b = auto_b_base, p = auto_p_base),
      auto_robust = c(b = auto_b_robust, p = auto_p_robust),
      auto_change_pct = auto_b_change,
      rai_baseline = c(b = rai_b_base, p = rai_p_base),
      rai_robust = c(b = rai_b_robust, p = rai_p_robust),
      rai_change_pct = rai_b_change,
      auto_robust = auto_robust,
      rai_robust = rai_robust,
      overall_robust = auto_robust && rai_robust
    )
  )

  saveRDS(results, file = paste0(base_dir, "/Results/Script07_", dataset_info$name, ".rds"))
  cat("Results saved\n")

  return(results)
}

################################################################################
# SCRIPT 08: Depression Ordinal Model - PUBLICATION-READY VERSION
################################################################################

run_script_08 <- function(dataset_info) {
  cat("\n--- Script 08: Depression Ordinal Model ---\n")
  cat("Dataset:", dataset_info$label, "\n\n")

  # Load and prepare data
  data <- read_sav(dataset_info$path) %>%
    mutate(Depression = DEP_Corrected) %>%
    filter(complete.cases(Autonomous, Controlled, Depression, Age, Sex)) %>%
    mutate(Depression_ord = ordered(Depression))

  cat("N =", nrow(data), "\n")
  cat("Depression levels:", nlevels(data$Depression_ord), "\n\n")

  # Define priors
  priors <- c(
    prior(student_t(3, 0, 2.5), class = "b"),
    prior(student_t(3, 0, 2.5), class = "Intercept")
  )

  # Fit model
  cat("Fitting Bayesian cumulative ordinal regression...\n")
  cat("Model: Depression ~ Autonomous + Controlled + Age + Sex\n")
  start_time <- Sys.time()

  fit <- brm(
    formula = Depression_ord ~ Autonomous + Controlled + Age + Sex,
    data = data, family = cumulative("probit"), prior = priors,
    chains = 4, iter = 3000, warmup = 1500, seed = 42,
    backend = "cmdstanr", cores = 4,
    control = list(adapt_delta = 0.95),
    file = paste0(base_dir, "/Models/Script08_", dataset_info$name),
    silent = 2, refresh = 0
  )

  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat("Model complete in", round(elapsed, 1), "minutes\n\n")

  # Convergence diagnostics
  cat("Convergence Diagnostics:\n")
  rhats <- rhat(fit)
  max_rhat <- max(rhats, na.rm = TRUE)
  cat("  Max R-hat:", round(max_rhat, 4))
  if (!is.na(max_rhat) && max_rhat < 1.01) {
    cat(" ✓ Excellent\n")
  } else {
    cat(" ⚠ Check convergence\n")
  }

  ess_ratios <- neff_ratio(fit)
  min_ess <- min(ess_ratios, na.rm = TRUE)
  cat("  Min ESS ratio:", round(min_ess, 3))
  if (!is.na(min_ess) && min_ess > 0.1) {
    cat(" ✓ Adequate\n")
  } else {
    cat(" ⚠ Low ESS\n")
  }

  np <- nuts_params(fit)
  divs <- sum(subset(np, Parameter == "divergent__")$Value)
  cat("  Divergences:", divs)
  if (divs == 0) {
    cat(" ✓ None\n\n")
  } else {
    cat(" ⚠ Present\n\n")
  }

  # Extract results
  post_sum <- posterior_summary(fit)
  fixed_effects <- c("b_Autonomous", "b_Controlled", "b_Age", "b_Sex")

  cat("Fixed Effects:\n")
  for (fx in fixed_effects) {
    if (fx %in% rownames(post_sum)) {
      est <- post_sum[fx, "Estimate"]
      lower <- post_sum[fx, "Q2.5"]
      upper <- post_sum[fx, "Q97.5"]
      pred_name <- sub("b_", "", fx)
      credible <- if (lower > 0) "+" else if (upper < 0) "-" else "~"
      cat(sprintf("  %-12s: %7.3f [%6.3f, %6.3f] %s\n",
                  pred_name, est, lower, upper, credible))
    }
  }
  cat("\n")

  # Key findings
  auto_est <- post_sum["b_Autonomous", "Estimate"]
  auto_lower <- post_sum["b_Autonomous", "Q2.5"]
  auto_upper <- post_sum["b_Autonomous", "Q97.5"]
  ctrl_est <- post_sum["b_Controlled", "Estimate"]
  ctrl_lower <- post_sum["b_Controlled", "Q2.5"]
  ctrl_upper <- post_sum["b_Controlled", "Q97.5"]

  cat("Key Findings:\n")
  cat("  Autonomous → Depression:", round(auto_est, 3),
      "[", round(auto_lower, 3), ",", round(auto_upper, 3), "]")
  if (auto_upper < 0) {
    cat(" ✓ Credibly negative (protective)\n")
  } else if (auto_lower > 0) {
    cat(" ⚠ Credibly positive (risk factor)\n")
  } else {
    cat(" ~ Uncertain\n")
  }

  cat("  Controlled → Depression:", round(ctrl_est, 3),
      "[", round(ctrl_lower, 3), ",", round(ctrl_upper, 3), "]")
  if (ctrl_lower > 0) {
    cat(" ✓ Credibly positive (risk factor)\n")
  } else if (ctrl_upper < 0) {
    cat(" ⚠ Credibly negative (protective)\n")
  } else {
    cat(" ~ Uncertain\n")
  }
  cat("\n")

  # Posterior predictive check
  cat("Generating posterior predictive check...\n")
  pp <- pp_check(fit, ndraws = 100, type = "bars")

  pp_file <- paste0(base_dir, "/Results/Script08_", dataset_info$name, "_pp_check.pdf")
  pdf(pp_file, width = 10, height = 6)
  print(pp +
        ggtitle(paste("Depression Ordinal Model:", dataset_info$label)) +
        theme_minimal())
  dev.off()
  cat("PP check saved:", pp_file, "\n")

  # Save summary results
  results <- list(
    dataset = dataset_info$name,
    n = nrow(data),
    fit = fit,
    post_summary = post_sum,
    fixed_effects = post_sum[fixed_effects, ],
    key_findings = list(
      autonomous = c(est = auto_est, lower = auto_lower, upper = auto_upper),
      controlled = c(est = ctrl_est, lower = ctrl_lower, upper = ctrl_upper)
    ),
    diagnostics = list(
      max_rhat = max_rhat,
      min_ess = min_ess,
      divergences = divs,
      elapsed_mins = as.numeric(elapsed)
    )
  )

  saveRDS(results, file = paste0(base_dir, "/Results/Script08_", dataset_info$name, ".rds"))
  cat("Results saved\n\n")

  return(results)
}

################################################################################
# SCRIPT 09: Item-level Depression (CESD) - PUBLICATION-READY VERSION
################################################################################

run_script_09 <- function(dataset_info) {
  cat("\n--- Script 09: Item-Level Depression Analysis ---\n")
  cat("Dataset:", dataset_info$label, "\n\n")

  # Load and prepare data
  data_raw <- read_sav(dataset_info$path)
  cesd_items <- paste0("CESD_", 1:10)

  # Handle different variable names across datasets:
  # Study 1 has SCReasons/nonSCReasons, Study 2 & Combined have Autonomous/Controlled
  if ("SCReasons" %in% names(data_raw)) {
    data_wide <- data_raw %>%
      dplyr::select(any_of(c("id", cesd_items, "SCReasons", "nonSCReasons", "Age", "Sex"))) %>%
      rename(Autonomous = SCReasons, Controlled = nonSCReasons) %>%
      filter(complete.cases(.))
  } else {
    data_wide <- data_raw %>%
      dplyr::select(any_of(c("id", cesd_items, "Autonomous", "Controlled", "Age", "Sex"))) %>%
      filter(complete.cases(.))
  }

  cat("N (people):", nrow(data_wide), "\n")
  cat("N (items):", length(cesd_items), "\n")
  cat("N (observations):", nrow(data_wide) * length(cesd_items), "\n\n")

  # Reshape to long format
  # Create Person ID using row_number() to avoid issues with missing/empty IDs
  data_long <- data_wide %>%
    mutate(Person = row_number()) %>%
    pivot_longer(cols = all_of(cesd_items), names_to = "Item", values_to = "Rating") %>%
    mutate(Rating = ordered(Rating), Person = as.factor(Person), Item = as.factor(Item))

  # Define priors
  priors <- c(
    prior(student_t(3, 0, 2.5), class = "b"),
    prior(student_t(3, 0, 2.5), class = "sd"),
    prior(student_t(3, 0, 2.5), class = "Intercept")
  )

  # Fit model
  cat("Fitting Bayesian cumulative ordinal regression...\n")
  start_time <- Sys.time()

  fit <- brm(
    formula = Rating ~ Autonomous + Controlled + Age + Sex + (1|Person) + (1|Item),
    data = data_long, family = cumulative("probit"), prior = priors,
    chains = 4, iter = 3000, warmup = 1500, seed = 42,
    backend = "cmdstanr", cores = 4,
    control = list(adapt_delta = 0.95),
    file = paste0(base_dir, "/Models/Script09_", dataset_info$name),
    silent = 2, refresh = 0
  )

  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat("Model complete in", round(elapsed, 1), "minutes\n\n")

  # Convergence diagnostics
  cat("Convergence Diagnostics:\n")
  rhats <- rhat(fit)
  max_rhat <- max(rhats, na.rm = TRUE)
  cat("  Max R-hat:", round(max_rhat, 4))
  if (!is.na(max_rhat) && max_rhat < 1.01) {
    cat(" ✓ Excellent\n")
  } else {
    cat(" ⚠ Check convergence\n")
  }

  ess_ratios <- neff_ratio(fit)
  min_ess <- min(ess_ratios, na.rm = TRUE)
  cat("  Min ESS ratio:", round(min_ess, 3))
  if (!is.na(min_ess) && min_ess > 0.1) {
    cat(" ✓ Adequate\n")
  } else {
    cat(" ⚠ Low ESS\n")
  }

  np <- nuts_params(fit)
  divs <- sum(subset(np, Parameter == "divergent__")$Value)
  cat("  Divergences:", divs)
  if (divs == 0) {
    cat(" ✓ None\n\n")
  } else {
    cat(" ⚠ Present\n\n")
  }

  # Extract results
  post_sum <- posterior_summary(fit)
  fixed_effects <- c("b_Autonomous", "b_Controlled", "b_Age", "b_Sex")

  cat("Fixed Effects:\n")
  for (fx in fixed_effects) {
    if (fx %in% rownames(post_sum)) {
      est <- post_sum[fx, "Estimate"]
      lower <- post_sum[fx, "Q2.5"]
      upper <- post_sum[fx, "Q97.5"]
      pred_name <- sub("b_", "", fx)
      credible <- if (lower > 0) "+" else if (upper < 0) "-" else "~"
      cat(sprintf("  %-12s: %7.3f [%6.3f, %6.3f] %s\n",
                  pred_name, est, lower, upper, credible))
    }
  }
  cat("\n")

  # Random effects
  re_person <- post_sum[grep("sd_Person__Intercept", rownames(post_sum)), "Estimate"]
  re_item <- post_sum[grep("sd_Item__Intercept", rownames(post_sum)), "Estimate"]

  cat("Random Effects:\n")
  cat("  Person SD:", round(re_person, 3), "\n")
  cat("  Item SD:", round(re_item, 3), "\n\n")

  # Key interpretation
  auto_est <- post_sum["b_Autonomous", "Estimate"]
  auto_lower <- post_sum["b_Autonomous", "Q2.5"]
  auto_upper <- post_sum["b_Autonomous", "Q97.5"]
  ctrl_est <- post_sum["b_Controlled", "Estimate"]
  ctrl_lower <- post_sum["b_Controlled", "Q2.5"]
  ctrl_upper <- post_sum["b_Controlled", "Q97.5"]

  cat("Key Findings:\n")
  cat("  Autonomous → Depression:", round(auto_est, 3),
      "[", round(auto_lower, 3), ",", round(auto_upper, 3), "]")
  if (auto_upper < 0) {
    cat(" ✓ Credibly negative\n")
  } else if (auto_lower > 0) {
    cat(" ⚠ Credibly positive\n")
  } else {
    cat(" ~ Uncertain\n")
  }

  cat("  Controlled → Depression:", round(ctrl_est, 3),
      "[", round(ctrl_lower, 3), ",", round(ctrl_upper, 3), "]")
  if (ctrl_lower > 0) {
    cat(" ✓ Credibly positive\n")
  } else if (ctrl_upper < 0) {
    cat(" ⚠ Credibly negative\n")
  } else {
    cat(" ~ Uncertain\n")
  }
  cat("\n")

  # Posterior predictive check
  cat("Generating posterior predictive check...\n")
  pp <- pp_check(fit, ndraws = 100, type = "bars")

  pp_file <- paste0(base_dir, "/Results/Script09_", dataset_info$name, "_pp_check.pdf")
  pdf(pp_file, width = 10, height = 6)
  print(pp +
        ggtitle(paste("CESD Item-Level Model:", dataset_info$label)) +
        theme_minimal())
  dev.off()
  cat("PP check saved:", pp_file, "\n")

  # Save summary results
  results <- list(
    fit = fit,
    data_long = data_long,
    post_summary = post_sum,
    fixed_effects = post_sum[fixed_effects, ],
    random_effects = list(person_sd = re_person, item_sd = re_item),
    key_findings = list(
      autonomous = c(est = auto_est, lower = auto_lower, upper = auto_upper),
      controlled = c(est = ctrl_est, lower = ctrl_lower, upper = ctrl_upper)
    ),
    diagnostics = list(
      max_rhat = max_rhat,
      min_ess = min_ess,
      divergences = divs,
      elapsed_mins = as.numeric(elapsed)
    )
  )

  saveRDS(results, file = paste0(base_dir, "/Results/Script09_", dataset_info$name, ".rds"))
  cat("Results saved\n")

  return(results)
}

################################################################################
# SCRIPT 10: Item-level MMLS Judgements - PUBLICATION-READY VERSION WITH DAG COMPARISON
################################################################################

run_script_10 <- function(dataset_info) {
  cat("\n--- Script 10: Item-Level MMLS General Meaning Judgements (3-Model DAG Comparison) ---\n")
  cat("Dataset:", dataset_info$label, "\n\n")

  # Load and prepare data
  data_raw <- read_sav(dataset_info$path)
  mil_items <- paste0("MIL", 1:4)  # MMLS General Meaning Judgements subscale (4 items)

  # Handle different variable names across datasets:
  # Study 1 has SCReasons/nonSCReasons, Study 2 & Combined have Autonomous/Controlled
  if ("SCReasons" %in% names(data_raw)) {
    data_wide <- data_raw %>%
      dplyr::select(any_of(c("id", mil_items, "SCReasons", "nonSCReasons",
                             "RAI", "DEP_Corrected", "GADscore", "Age", "Sex"))) %>%
      rename(Autonomous = SCReasons, Controlled = nonSCReasons,
             Depression = DEP_Corrected, Anxiety = GADscore) %>%
      filter(complete.cases(.))
  } else {
    data_wide <- data_raw %>%
      dplyr::select(any_of(c("id", mil_items, "Autonomous", "Controlled",
                             "RAI", "DEP_Corrected", "GADscore", "Age", "Sex"))) %>%
      rename(Depression = DEP_Corrected, Anxiety = GADscore) %>%
      filter(complete.cases(.))
  }

  cat("N (people):", nrow(data_wide), "\n")
  cat("N (items):", length(mil_items), "\n")
  cat("N (observations):", nrow(data_wide) * length(mil_items), "\n\n")

  # Reshape to long format
  # Create Person ID using row_number() to avoid issues with missing/empty IDs
  data_long <- data_wide %>%
    mutate(Person = row_number()) %>%
    pivot_longer(cols = all_of(mil_items), names_to = "Item", values_to = "Rating") %>%
    mutate(Rating = ordered(Rating), Person = as.factor(Person), Item = as.factor(Item))

  # Define priors
  priors <- c(
    prior(student_t(3, 0, 2.5), class = "b"),
    prior(student_t(3, 0, 2.5), class = "sd"),
    prior(student_t(3, 0, 2.5), class = "Intercept")
  )

  cat("Testing 3 competing models (item-level DAG comparison):\n")
  cat("  Model 1: Autonomous only + controls (DAG 2)\n")
  cat("  Model 2: RAI only + controls (DAG 1)\n")
  cat("  Model 3: Both Autonomous + Controlled + controls (DAG 3)\n\n")

  # ============================================================================
  # MODEL 1: AUTONOMOUS ONLY (DAG 2 ITEM-LEVEL)
  # ============================================================================
  cat("MODEL 1: Autonomous Only\n")
  cat("Formula: Rating ~ Autonomous + Depression + Anxiety + Age + Sex + (1|Person) + (1|Item)\n")
  start_time1 <- Sys.time()

  fit1 <- brm(
    Rating ~ Autonomous + Depression + Anxiety + Age + Sex + (1|Person) + (1|Item),
    data = data_long, family = cumulative("probit"), prior = priors,
    chains = 4, iter = 2000, warmup = 1000, seed = 42,
    backend = "cmdstanr", cores = 4,
    control = list(adapt_delta = 0.95),
    file = paste0(base_dir, "/Models/Script10_Model1_", dataset_info$name),
    silent = 2, refresh = 0
  )

  elapsed1 <- difftime(Sys.time(), start_time1, units = "mins")
  cat("Complete in", round(elapsed1, 1), "minutes\n")

  # Diagnostics Model 1
  rhats1 <- rhat(fit1)
  max_rhat1 <- max(rhats1, na.rm = TRUE)
  ess1 <- neff_ratio(fit1)
  min_ess1 <- min(ess1, na.rm = TRUE)
  np1 <- nuts_params(fit1)
  divs1 <- sum(subset(np1, Parameter == "divergent__")$Value)

  cat("  R-hat:", round(max_rhat1, 4))
  if (!is.na(max_rhat1) && max_rhat1 < 1.01) cat(" ✓\n") else cat(" ⚠\n")
  cat("  ESS:", round(min_ess1, 3))
  if (!is.na(min_ess1) && min_ess1 > 0.1) cat(" ✓\n") else cat(" ⚠\n")
  cat("  Divergences:", divs1)
  if (divs1 == 0) cat(" ✓\n") else cat(" ⚠\n")

  # LOO-IC Model 1
  loo1 <- loo(fit1)
  cat("  LOO-IC:", round(loo1$estimates["looic", "Estimate"], 1), "\n\n")

  # ============================================================================
  # MODEL 2: RAI ONLY (DAG 1 ITEM-LEVEL)
  # ============================================================================
  cat("MODEL 2: RAI Only\n")
  cat("Formula: Rating ~ RAI + Depression + Anxiety + Age + Sex + (1|Person) + (1|Item)\n")
  start_time2 <- Sys.time()

  fit2 <- brm(
    Rating ~ RAI + Depression + Anxiety + Age + Sex + (1|Person) + (1|Item),
    data = data_long, family = cumulative("probit"), prior = priors,
    chains = 4, iter = 2000, warmup = 1000, seed = 42,
    backend = "cmdstanr", cores = 4,
    control = list(adapt_delta = 0.95),
    file = paste0(base_dir, "/Models/Script10_Model2_", dataset_info$name),
    silent = 2, refresh = 0
  )

  elapsed2 <- difftime(Sys.time(), start_time2, units = "mins")
  cat("Complete in", round(elapsed2, 1), "minutes\n")

  # Diagnostics Model 2
  rhats2 <- rhat(fit2)
  max_rhat2 <- max(rhats2, na.rm = TRUE)
  ess2 <- neff_ratio(fit2)
  min_ess2 <- min(ess2, na.rm = TRUE)
  np2 <- nuts_params(fit2)
  divs2 <- sum(subset(np2, Parameter == "divergent__")$Value)

  cat("  R-hat:", round(max_rhat2, 4))
  if (!is.na(max_rhat2) && max_rhat2 < 1.01) cat(" ✓\n") else cat(" ⚠\n")
  cat("  ESS:", round(min_ess2, 3))
  if (!is.na(min_ess2) && min_ess2 > 0.1) cat(" ✓\n") else cat(" ⚠\n")
  cat("  Divergences:", divs2)
  if (divs2 == 0) cat(" ✓\n") else cat(" ⚠\n")

  # LOO-IC Model 2
  loo2 <- loo(fit2)
  cat("  LOO-IC:", round(loo2$estimates["looic", "Estimate"], 1), "\n\n")

  # ============================================================================
  # MODEL 3: BOTH PREDICTORS (DAG 3 ITEM-LEVEL)
  # ============================================================================
  cat("MODEL 3: Both Autonomous + Controlled\n")
  cat("Formula: Rating ~ Autonomous + Controlled + Depression + Anxiety + Age + Sex + (1|Person) + (1|Item)\n")
  start_time3 <- Sys.time()

  fit3 <- brm(
    Rating ~ Autonomous + Controlled + Depression + Anxiety + Age + Sex + (1|Person) + (1|Item),
    data = data_long, family = cumulative("probit"), prior = priors,
    chains = 4, iter = 2000, warmup = 1000, seed = 42,
    backend = "cmdstanr", cores = 4,
    control = list(adapt_delta = 0.95),
    file = paste0(base_dir, "/Models/Script10_Model3_", dataset_info$name),
    silent = 2, refresh = 0
  )

  elapsed3 <- difftime(Sys.time(), start_time3, units = "mins")
  cat("Complete in", round(elapsed3, 1), "minutes\n")

  # Diagnostics Model 3
  rhats3 <- rhat(fit3)
  max_rhat3 <- max(rhats3, na.rm = TRUE)
  ess3 <- neff_ratio(fit3)
  min_ess3 <- min(ess3, na.rm = TRUE)
  np3 <- nuts_params(fit3)
  divs3 <- sum(subset(np3, Parameter == "divergent__")$Value)

  cat("  R-hat:", round(max_rhat3, 4))
  if (!is.na(max_rhat3) && max_rhat3 < 1.01) cat(" ✓\n") else cat(" ⚠\n")
  cat("  ESS:", round(min_ess3, 3))
  if (!is.na(min_ess3) && min_ess3 > 0.1) cat(" ✓\n") else cat(" ⚠\n")
  cat("  Divergences:", divs3)
  if (divs3 == 0) cat(" ✓\n") else cat(" ⚠\n")

  # LOO-IC Model 3
  loo3 <- loo(fit3)
  cat("  LOO-IC:", round(loo3$estimates["looic", "Estimate"], 1), "\n\n")

  # ============================================================================
  # MODEL COMPARISON
  # ============================================================================
  cat("MODEL COMPARISON (LOO-IC):\n")

  loo_comp <- loo_compare(loo1, loo2, loo3)
  print(loo_comp)
  cat("\n")

  best_model <- rownames(loo_comp)[1]
  cat("Best model:", best_model, "\n")

  if (grepl("fit1", best_model)) {
    cat("  → Autonomous-only model wins (supports DAG 2 / H1)\n\n")
  } else if (grepl("fit2", best_model)) {
    cat("  → RAI-only model wins (supports DAG 1 / traditional SDT)\n\n")
  } else {
    cat("  → Both-predictors model wins (supports DAG 3 / dual pathway)\n\n")
  }

  # Extract key parameters from best model
  cat("KEY FINDINGS (Model 1 - Autonomous Only):\n")
  post_sum1 <- posterior_summary(fit1)

  auto_est <- post_sum1["b_Autonomous", "Estimate"]
  auto_lower <- post_sum1["b_Autonomous", "Q2.5"]
  auto_upper <- post_sum1["b_Autonomous", "Q97.5"]

  cat("  Autonomous → MMLS Judgements:", round(auto_est, 3),
      "[", round(auto_lower, 3), ",", round(auto_upper, 3), "]")
  if (auto_lower > 0) {
    cat(" ✓ Credibly positive\n")
  } else if (auto_upper < 0) {
    cat(" ⚠ Credibly negative\n")
  } else {
    cat(" ~ Uncertain\n")
  }

  cat("\nKEY FINDINGS (Model 3 - Both Predictors):\n")
  post_sum3 <- posterior_summary(fit3)

  auto_est3 <- post_sum3["b_Autonomous", "Estimate"]
  auto_lower3 <- post_sum3["b_Autonomous", "Q2.5"]
  auto_upper3 <- post_sum3["b_Autonomous", "Q97.5"]
  ctrl_est3 <- post_sum3["b_Controlled", "Estimate"]
  ctrl_lower3 <- post_sum3["b_Controlled", "Q2.5"]
  ctrl_upper3 <- post_sum3["b_Controlled", "Q97.5"]

  cat("  Autonomous → MMLS Judgements:", round(auto_est3, 3),
      "[", round(auto_lower3, 3), ",", round(auto_upper3, 3), "]")
  if (auto_lower3 > 0) {
    cat(" ✓ Credibly positive\n")
  } else {
    cat(" ~ Uncertain\n")
  }

  cat("  Controlled → MMLS Judgements:", round(ctrl_est3, 3),
      "[", round(ctrl_lower3, 3), ",", round(ctrl_upper3, 3), "]")
  if (ctrl_upper3 < 0) {
    cat(" ✓ Credibly negative\n")
  } else if (ctrl_lower3 > 0) {
    cat(" ⚠ Credibly positive\n")
  } else {
    cat(" ~ Uncertain\n")
  }
  cat("\n")

  # Posterior predictive checks
  cat("Generating posterior predictive checks...\n")
  pp1 <- pp_check(fit1, ndraws = 100, type = "bars")
  pp3 <- pp_check(fit3, ndraws = 100, type = "bars")

  pp_file1 <- paste0(base_dir, "/Results/Script10_Model1_", dataset_info$name, "_pp_check.pdf")
  pdf(pp_file1, width = 10, height = 6)
  print(pp1 +
        ggtitle(paste("MMLS Judgements Item-Level Model 1 (Autonomous):", dataset_info$label)) +
        theme_minimal())
  dev.off()

  pp_file3 <- paste0(base_dir, "/Results/Script10_Model3_", dataset_info$name, "_pp_check.pdf")
  pdf(pp_file3, width = 10, height = 6)
  print(pp3 +
        ggtitle(paste("MMLS Judgements Item-Level Model 3 (Both):", dataset_info$label)) +
        theme_minimal())
  dev.off()
  cat("PP checks saved\n")

  # Save comprehensive results
  results <- list(
    model1 = list(
      fit = fit1,
      post_summary = post_sum1,
      loo = loo1,
      diagnostics = list(max_rhat = max_rhat1, min_ess = min_ess1,
                         divergences = divs1, elapsed_mins = as.numeric(elapsed1))
    ),
    model2 = list(
      fit = fit2,
      loo = loo2,
      diagnostics = list(max_rhat = max_rhat2, min_ess = min_ess2,
                         divergences = divs2, elapsed_mins = as.numeric(elapsed2))
    ),
    model3 = list(
      fit = fit3,
      post_summary = post_sum3,
      loo = loo3,
      diagnostics = list(max_rhat = max_rhat3, min_ess = min_ess3,
                         divergences = divs3, elapsed_mins = as.numeric(elapsed3))
    ),
    comparison = list(
      loo_compare = loo_comp,
      best_model = best_model
    ),
    data_long = data_long,
    key_findings = list(
      model1_autonomous = c(est = auto_est, lower = auto_lower, upper = auto_upper),
      model3_autonomous = c(est = auto_est3, lower = auto_lower3, upper = auto_upper3),
      model3_controlled = c(est = ctrl_est3, lower = ctrl_lower3, upper = ctrl_upper3)
    )
  )

  saveRDS(results, file = paste0(base_dir, "/Results/Script10_", dataset_info$name, ".rds"))
  cat("Results saved\n\n")

  return(results)
}

################################################################################
# MAIN EXECUTION LOOP
################################################################################

all_results <- list()
completion_summary <- list()

for (ds in datasets) {
  cat("\n################################################################\n")
  cat("PROCESSING:", ds$label, "\n")
  cat("################################################################\n\n")

  dataset_results <- list()
  dataset_results$info <- ds

  # Run all 10 scripts (skip 01 - data prep)
  scripts <- list(
    list(num = "02", func = run_script_02, name = "H1 Comparative Regression"),
    list(num = "03", func = run_script_03, name = "H1 Mediation Analysis (RAI→Auto→Meaning)"),
    list(num = "03b", func = run_script_03b, name = "Supplementary Mediation (Auto→Dep→Meaning)"),
    list(num = "04", func = run_script_04, name = "H1 Bayesian Ordinal"),
    list(num = "05", func = run_script_05, name = "H2 Controlled Motivation"),
    list(num = "06", func = run_script_06, name = "H3 SEM DAG Comparison"),
    list(num = "07", func = run_script_07, name = "H4 Robustness Outliers"),
    list(num = "08", func = run_script_08, name = "Depression Ordinal Model"),
    list(num = "09", func = run_script_09, name = "Item-level Depression"),
    list(num = "10", func = run_script_10, name = "Item-level MMLS Judgements")
  )

  for (script in scripts) {
    script_key <- paste0("script", script$num)
    dataset_results[[script_key]] <- safe_run(script$func, ds, paste0("Script ", script$num))

    # Track completion
    if (dataset_results[[script_key]]$success) {
      completion_summary[[paste0(ds$name, "_", script$num)]] <- "✅"
    } else {
      completion_summary[[paste0(ds$name, "_", script$num)]] <- "❌"
    }
  }

  all_results[[ds$name]] <- dataset_results

  log_progress(paste0("Completed all scripts for ", ds$label))
}

################################################################################
# SAVE FINAL RESULTS
################################################################################

saveRDS(all_results, file = paste0(base_dir, "/ALL_RESULTS_MASTER.rds"))

overall_end <- Sys.time()
total_duration <- difftime(overall_end, overall_start, units = "hours")

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n\n################################################################\n")
cat("################################################################\n")
cat("ALL ANALYSES COMPLETED!\n")
cat("################################################################\n")
cat("################################################################\n\n")

cat("Total execution time:", round(total_duration, 2), "hours\n\n")

cat("COMPLETION SUMMARY:\n")
cat("==================\n\n")

# Print completion matrix
cat("         | Study1 | Study2 | Combined\n")
cat("---------+--------+--------+---------\n")
for (i in 2:10) {
  script_num <- sprintf("%02d", i)
  cat(sprintf("Script %s | %s      | %s      | %s\n",
              script_num,
              completion_summary[[paste0("Study1_", script_num)]],
              completion_summary[[paste0("Study2_", script_num)]],
              completion_summary[[paste0("Combined_", script_num)]]))
}

cat("\n\nERROR SUMMARY:\n")
cat("==============\n")
if (length(error_log) == 0) {
  cat("✅ NO ERRORS - All 27 analyses completed successfully!\n")
} else {
  cat("⚠️ ", length(error_log), " error(s) encountered:\n\n")
  for (err in error_log) {
    cat(err, "\n")
  }
  cat("\nSee ERROR_LOG.txt for details\n")
}

cat("\n\nOUTPUT LOCATIONS:\n")
cat("=================\n")
cat("Main folder:", base_dir, "/\n")
cat("  - output.log          (this log)\n")
cat("  - PROGRESS.txt        (progress tracking)\n")
cat("  - ERROR_LOG.txt       (error details)\n")
cat("  - Results/            (27 .rds files - one per script per dataset)\n")
cat("  - Models/             (12 Bayesian model files)\n")
cat("  - ALL_RESULTS_MASTER.rds (master results file)\n\n")

log_progress("MASTER SCRIPT COMPLETED")

cat("✅ DONE!\n")
