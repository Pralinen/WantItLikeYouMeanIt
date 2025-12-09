################################################################################
# SCRIPT 11b: BAYESIAN LATENT CONSTRUCT MODELING (blavaan)
################################################################################
#
# Purpose: Bayesian version of Script 11 using blavaan
#          Full Bayesian SEM with latent variables and posterior inference
#
# This addresses the "real Bayesian" requirement by:
#   - Using MCMC sampling for latent factor models
#   - Providing posterior distributions (not just point estimates)
#   - Allowing for Bayesian model comparison
#
# Structure:
#   PART 1: Full Measurement Model CFA (all latent constructs)
#   PART 2: DAG 1 - Traditional SDT (RAI → Meaning)
#   PART 3: DAG 2 - Your Hypothesis (Autonomous → Meaning, Controlled = 0)
#   PART 4: DAG 3 - Dual Pathway (Both AUTO and CTRL)
#   PART 5: Bayesian Path Decomposition (Mediation Analysis)
#   PART 6: Model Comparison
#
# Uses the same goal-level parcelling approach as Script 11 (proven good fit)
#
################################################################################

# library(haven) # Not needed for CSV
library(dplyr)
library(lavaan)
library(blavaan)

# Helper function to extract parameter estimates (handles blavaan column names)
get_param_est <- function(params, col = "est") {
  if (col == "est") {
    if ("est" %in% names(params)) return(params$est)
    if ("Estimate" %in% names(params)) return(params$Estimate)
    return(NA)
  } else if (col == "ci.lower") {
    if ("ci.lower" %in% names(params)) return(params$ci.lower)
    if ("pi.lower" %in% names(params)) return(params$pi.lower)
    return(NA)
  } else if (col == "ci.upper") {
    if ("ci.upper" %in% names(params)) return(params$ci.upper)
    if ("pi.upper" %in% names(params)) return(params$pi.upper)
    return(NA)
  }
  return(NA)
}

cat("================================================================================\n")
cat("SCRIPT 11b: BAYESIAN LATENT CONSTRUCT MODELING (blavaan)\n")
cat("================================================================================\n\n")

# Set paths
base_dir <- "."
data_dir <- "."

# Create output directory
if (!dir.exists("Results/Latent")) {
  dir.create("Results/Latent", recursive = TRUE)
}

# MCMC settings (kept light for reasonable runtime)
n_chains <- 3
n_burnin <- 1000
n_sample <- 1000

cat("MCMC Settings: Chains =", n_chains, ", Burnin =", n_burnin, ", Samples =", n_sample, "\n\n")

################################################################################
# LOAD AND PREPARE DATA
################################################################################

cat("Loading Combined dataset...\n")
combined <- read.csv("dataset.csv")
cat("  N =", nrow(combined), "\n\n")

################################################################################
# CREATE GOAL-LEVEL PARCELS (same as Script 11)
################################################################################

cat("Creating goal-level parcels for motivation...\n")

for (i in 1:6) {
  # Autonomous parcel: mean of Identified (_2) and Intrinsic (_4)
  combined[[paste0("Goal_", i, "_Auto")]] <- rowMeans(
    combined[, c(paste0("GSC", i, "_2"), paste0("GSC", i, "_4"))], na.rm = TRUE
  )
  # Controlled parcel: mean of External (_1) and Introjected (_3)
  combined[[paste0("Goal_", i, "_Ctrl")]] <- rowMeans(
    combined[, c(paste0("GSC", i, "_1"), paste0("GSC", i, "_3"))], na.rm = TRUE
  )
}

cat("  Created 12 goal-level parcels\n\n")

# Define CESD items (use reverse-coded item 8)
cesd_items <- c("CESD_1", "CESD_2", "CESD_3", "CESD_4", "CESD_5",
                "CESD_6", "CESD_7", "CESD_8R", "CESD_9", "CESD_10")
cesd_items <- cesd_items[cesd_items %in% names(combined)]

################################################################################
# PART 1: FULL MEASUREMENT MODEL CFA (All Latent Constructs)
################################################################################

cat("================================================================================\n")
cat("PART 1: BAYESIAN CFA - FULL MEASUREMENT MODEL\n")
cat("================================================================================\n\n")

full_cfa_model <- paste0('
  # Meaning in Life
  MIL =~ MIL1 + MIL2 + MIL3 + MIL4

  # Autonomous Motivation (6 goal-level parcels)
  AUTO =~ Goal_1_Auto + Goal_2_Auto + Goal_3_Auto + Goal_4_Auto + Goal_5_Auto + Goal_6_Auto

  # Controlled Motivation (6 goal-level parcels)
  CTRL =~ Goal_1_Ctrl + Goal_2_Ctrl + Goal_3_Ctrl + Goal_4_Ctrl + Goal_5_Ctrl + Goal_6_Ctrl

  # Depression
  DEP =~ ', paste(cesd_items, collapse = " + "), '
')

cat("Fitting Full Measurement Model CFA...\n")
cat("  4 latent factors: MIL, AUTO, CTRL, DEP\n\n")

bcfa_full <- tryCatch({
  bcfa(full_cfa_model, data = combined,
       n.chains = n_chains, burnin = n_burnin, sample = n_sample,
       std.lv = TRUE, seed = 42)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(bcfa_full)) {
  cat("\n--- Bayesian Full CFA Results ---\n\n")

  # Summary
  print(summary(bcfa_full))

  # Fit indices
  cat("\nBayesian Fit Indices:\n")
  fit_full <- blavFitIndices(bcfa_full)
  print(fit_full)

  # Extract factor correlations
  cat("\nFactor Correlations (Posterior Estimates):\n")
  params <- parameterEstimates(bcfa_full)
  cors <- params[params$op == "~~" & params$lhs != params$rhs, ]
  if (nrow(cors) > 0) {
    for (i in 1:nrow(cors)) {
      est_val <- if("est" %in% names(cors)) cors$est[i] else cors$Estimate[i]
      ci_lo <- if("ci.lower" %in% names(cors)) cors$ci.lower[i] else cors$Post.SD[i]
      ci_hi <- if("ci.upper" %in% names(cors)) cors$ci.upper[i] else NA
      cat("  ", cors$lhs[i], "<->", cors$rhs[i], ":",
          round(est_val, 3), "\n")
    }
  }
  cat("\n")
}

################################################################################
# PART 2: DAG 1 - Traditional SDT (RAI → Meaning)
################################################################################

cat("================================================================================\n")
cat("PART 2: BAYESIAN SEM - DAG 1 (Traditional SDT: RAI → Meaning)\n")
cat("================================================================================\n\n")

dag1_model <- paste0('
  # Measurement Model
  MIL =~ MIL1 + MIL2 + MIL3 + MIL4
  DEP =~ ', paste(cesd_items, collapse = " + "), '

  # Structural Model - RAI only (balance score)
  MIL ~ RAI + DEP + Age + Sex
')

cat("Fitting DAG 1: RAI (balance) predicts Meaning...\n\n")

bsem_dag1 <- tryCatch({
  bsem(dag1_model, data = combined,
       n.chains = n_chains, burnin = n_burnin, sample = n_sample,
       std.lv = TRUE, seed = 42)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(bsem_dag1)) {
  cat("\n--- Bayesian SEM Results: DAG 1 ---\n\n")
  print(summary(bsem_dag1))

  # Extract RAI effect
  dag1_params <- parameterEstimates(bsem_dag1)
  rai_effect <- dag1_params[dag1_params$lhs == "MIL" & dag1_params$rhs == "RAI", ]
  rai_est <- get_param_est(rai_effect, "est")
  rai_lo <- get_param_est(rai_effect, "ci.lower")
  rai_hi <- get_param_est(rai_effect, "ci.upper")

  cat("\nKEY RESULT - RAI → Meaning:\n")
  cat("  Posterior Mean:", round(rai_est, 3), "\n")
  cat("  95% Credible Interval: [", round(rai_lo, 3), ",", round(rai_hi, 3), "]\n")
  if (!is.na(rai_lo) && rai_lo > 0) {
    cat("  → CREDIBLY POSITIVE\n")
  } else if (!is.na(rai_hi) && rai_hi < 0) {
    cat("  → CREDIBLY NEGATIVE\n")
  } else {
    cat("  → INCLUDES 0\n")
  }

  cat("\nBayesian Fit Indices:\n")
  print(blavFitIndices(bsem_dag1))
  cat("\n")
}

################################################################################
# PART 3: DAG 2 - Your Hypothesis (Autonomous → Meaning, Controlled = 0)
################################################################################

cat("================================================================================\n")
cat("PART 3: BAYESIAN SEM - DAG 2 (YOUR HYPOTHESIS: Autonomous → Meaning)\n")
cat("================================================================================\n\n")
cat("Controlled motivation constrained to 0 (not in model)\n\n")

dag2_model <- paste0('
  # Measurement Model
  MIL =~ MIL1 + MIL2 + MIL3 + MIL4
  AUTO =~ Goal_1_Auto + Goal_2_Auto + Goal_3_Auto + Goal_4_Auto + Goal_5_Auto + Goal_6_Auto
  DEP =~ ', paste(cesd_items, collapse = " + "), '

  # Structural Model - Autonomous only (Controlled OMITTED = constrained to 0)
  MIL ~ AUTO + DEP + Age + Sex
')

cat("Fitting DAG 2: Autonomous predicts Meaning (Controlled = 0)...\n\n")

bsem_dag2 <- tryCatch({
  bsem(dag2_model, data = combined,
       n.chains = n_chains, burnin = n_burnin, sample = n_sample,
       std.lv = TRUE, seed = 42)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(bsem_dag2)) {
  cat("\n--- Bayesian SEM Results: DAG 2 ---\n\n")
  print(summary(bsem_dag2))

  # Extract AUTO effect
  dag2_params <- parameterEstimates(bsem_dag2)
  auto_effect <- dag2_params[dag2_params$lhs == "MIL" & dag2_params$rhs == "AUTO", ]
  auto_est <- get_param_est(auto_effect, "est")
  auto_lo <- get_param_est(auto_effect, "ci.lower")
  auto_hi <- get_param_est(auto_effect, "ci.upper")

  cat("\nKEY RESULT - Autonomous → Meaning:\n")
  cat("  Posterior Mean:", round(auto_est, 3), "\n")
  cat("  95% Credible Interval: [", round(auto_lo, 3), ",", round(auto_hi, 3), "]\n")
  if (!is.na(auto_lo) && auto_lo > 0) {
    cat("  → CREDIBLY POSITIVE\n")
  } else {
    cat("  → INCLUDES 0\n")
  }

  cat("\nBayesian Fit Indices:\n")
  print(blavFitIndices(bsem_dag2))
  cat("\n")
}

################################################################################
# PART 4: DAG 3 - Dual Pathway (Both AUTO and CTRL)
################################################################################

cat("================================================================================\n")
cat("PART 4: BAYESIAN SEM - DAG 3 (Dual Pathway: Both AUTO and CTRL)\n")
cat("================================================================================\n\n")
cat("Testing if Controlled adds anything beyond Autonomous\n\n")

dag3_model <- paste0('
  # Measurement Model
  MIL =~ MIL1 + MIL2 + MIL3 + MIL4
  AUTO =~ Goal_1_Auto + Goal_2_Auto + Goal_3_Auto + Goal_4_Auto + Goal_5_Auto + Goal_6_Auto
  CTRL =~ Goal_1_Ctrl + Goal_2_Ctrl + Goal_3_Ctrl + Goal_4_Ctrl + Goal_5_Ctrl + Goal_6_Ctrl
  DEP =~ ', paste(cesd_items, collapse = " + "), '

  # Structural Model - Both pathways
  MIL ~ AUTO + CTRL + DEP + Age + Sex
')

cat("Fitting DAG 3: Both Autonomous and Controlled predict Meaning...\n\n")

bsem_dag3 <- tryCatch({
  bsem(dag3_model, data = combined,
       n.chains = n_chains, burnin = n_burnin, sample = n_sample,
       std.lv = TRUE, seed = 42)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(bsem_dag3)) {
  cat("\n--- Bayesian SEM Results: DAG 3 ---\n\n")
  print(summary(bsem_dag3))

  # Extract effects
  dag3_params <- parameterEstimates(bsem_dag3)
  auto_effect <- dag3_params[dag3_params$lhs == "MIL" & dag3_params$rhs == "AUTO", ]
  ctrl_effect <- dag3_params[dag3_params$lhs == "MIL" & dag3_params$rhs == "CTRL", ]

  auto_est <- get_param_est(auto_effect, "est")
  auto_lo <- get_param_est(auto_effect, "ci.lower")
  auto_hi <- get_param_est(auto_effect, "ci.upper")
  ctrl_est <- get_param_est(ctrl_effect, "est")
  ctrl_lo <- get_param_est(ctrl_effect, "ci.lower")
  ctrl_hi <- get_param_est(ctrl_effect, "ci.upper")

  cat("\nKEY RESULTS:\n")
  cat("\nAutonomous → Meaning:\n")
  cat("  Posterior Mean:", round(auto_est, 3), "\n")
  cat("  95% CI: [", round(auto_lo, 3), ",", round(auto_hi, 3), "]\n")
  if (!is.na(auto_lo) && auto_lo > 0) {
    cat("  → CREDIBLY POSITIVE\n")
  } else {
    cat("  → Includes 0\n")
  }

  cat("\nControlled → Meaning:\n")
  cat("  Posterior Mean:", round(ctrl_est, 3), "\n")
  cat("  95% CI: [", round(ctrl_lo, 3), ",", round(ctrl_hi, 3), "]\n")
  if (!is.na(ctrl_lo) && !is.na(ctrl_hi) && (ctrl_lo > 0 || ctrl_hi < 0)) {
    cat("  → CREDIBLY DIFFERENT FROM 0\n")
  } else {
    cat("  → INCLUDES 0 (not credible) - supports constraining to 0!\n")
  }

  cat("\nBayesian Fit Indices:\n")
  print(blavFitIndices(bsem_dag3))
  cat("\n")
}

################################################################################
# PART 5: BAYESIAN PATH DECOMPOSITION (Mediation Analysis)
################################################################################

cat("================================================================================\n")
cat("PART 5: BAYESIAN PATH DECOMPOSITION (Mediation Analysis)\n")
cat("================================================================================\n\n")
cat("Testing: Does Autonomous mediate the RAI → Meaning relationship?\n\n")

path_model <- paste0('
  # Measurement Model
  MIL =~ MIL1 + MIL2 + MIL3 + MIL4
  AUTO =~ Goal_1_Auto + Goal_2_Auto + Goal_3_Auto + Goal_4_Auto + Goal_5_Auto + Goal_6_Auto
  DEP =~ ', paste(cesd_items, collapse = " + "), '

  # Structural Model (Causal Paths)
  # Path a: RAI → Autonomous
  AUTO ~ a*RAI + DEP + Age + Sex

  # Path b: Autonomous → Meaning (controlling for RAI)
  # Path c_prime: Direct effect of RAI on Meaning
  MIL ~ b*AUTO + c_prime*RAI + DEP + Age + Sex

  # Effect Decomposition
  indirect := a * b
  direct := c_prime
  total := c_prime + (a * b)
')

cat("Fitting Bayesian Path Decomposition Model...\n\n")

bsem_path <- tryCatch({
  bsem(path_model, data = combined,
       n.chains = n_chains, burnin = n_burnin, sample = n_sample,
       std.lv = TRUE, seed = 42)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(bsem_path)) {
  cat("\n--- Bayesian Path Decomposition Results ---\n\n")
  print(summary(bsem_path))

  # Extract path coefficients
  path_params <- parameterEstimates(bsem_path)

  a_path <- path_params[path_params$label == "a", ]
  b_path <- path_params[path_params$label == "b", ]
  indirect <- path_params[path_params$label == "indirect", ]
  direct <- path_params[path_params$label == "direct", ]
  total <- path_params[path_params$label == "total", ]

  # Get values using helper
  a_est <- get_param_est(a_path, "est"); a_lo <- get_param_est(a_path, "ci.lower"); a_hi <- get_param_est(a_path, "ci.upper")
  b_est <- get_param_est(b_path, "est"); b_lo <- get_param_est(b_path, "ci.lower"); b_hi <- get_param_est(b_path, "ci.upper")
  ind_est <- get_param_est(indirect, "est"); ind_lo <- get_param_est(indirect, "ci.lower"); ind_hi <- get_param_est(indirect, "ci.upper")
  dir_est <- get_param_est(direct, "est"); dir_lo <- get_param_est(direct, "ci.lower"); dir_hi <- get_param_est(direct, "ci.upper")
  tot_est <- get_param_est(total, "est"); tot_lo <- get_param_est(total, "ci.lower"); tot_hi <- get_param_est(total, "ci.upper")

  cat("\n========================================\n")
  cat("CAUSAL PATH ESTIMATES (Posterior)\n")
  cat("========================================\n\n")

  cat("Path a (RAI → AUTO):\n")
  cat("  Posterior Mean:", round(a_est, 3), "\n")
  cat("  95% CI: [", round(a_lo, 3), ",", round(a_hi, 3), "]\n\n")

  cat("Path b (AUTO → MIL):\n")
  cat("  Posterior Mean:", round(b_est, 3), "\n")
  cat("  95% CI: [", round(b_lo, 3), ",", round(b_hi, 3), "]\n\n")

  cat("========================================\n")
  cat("EFFECT DECOMPOSITION\n")
  cat("========================================\n\n")

  cat("Indirect Effect (a × b):\n")
  cat("  Posterior Mean:", round(ind_est, 3), "\n")
  cat("  95% CI: [", round(ind_lo, 3), ",", round(ind_hi, 3), "]\n")
  ind_credible <- !is.na(ind_lo) && !is.na(ind_hi) && (ind_lo > 0 || ind_hi < 0)
  if (ind_credible) {
    cat("  → CREDIBLY DIFFERENT FROM 0 (significant mediation)\n\n")
  } else {
    cat("  → Includes 0\n\n")
  }

  cat("Direct Effect (c'):\n")
  cat("  Posterior Mean:", round(dir_est, 3), "\n")
  cat("  95% CI: [", round(dir_lo, 3), ",", round(dir_hi, 3), "]\n")
  dir_credible <- !is.na(dir_lo) && !is.na(dir_hi) && (dir_lo > 0 || dir_hi < 0)
  if (dir_credible) {
    cat("  → CREDIBLY DIFFERENT FROM 0 (partial mediation)\n\n")
  } else {
    cat("  → INCLUDES 0 (full mediation - supports your hypothesis!)\n\n")
  }

  cat("Total Effect:\n")
  cat("  Posterior Mean:", round(tot_est, 3), "\n")
  cat("  95% CI: [", round(tot_lo, 3), ",", round(tot_hi, 3), "]\n\n")

  # Proportion mediated
  if (!is.na(tot_est) && tot_est != 0) {
    prop_med <- ind_est / tot_est
    cat("Proportion Mediated:", round(prop_med * 100, 1), "%\n\n")
  }

  cat("========================================\n")
  cat("INTERPRETATION\n")
  cat("========================================\n")
  if (ind_credible && !dir_credible) {
    cat("FULL MEDIATION SUPPORTED!\n")
    cat("→ RAI's effect on Meaning operates entirely through Autonomous motivation\n")
    cat("→ This supports your DAG 2 hypothesis\n")
  } else if (ind_credible) {
    cat("PARTIAL MEDIATION\n")
    cat("→ Autonomous partially mediates RAI → Meaning\n")
    cat("→ But RAI retains some direct effect\n")
  } else {
    cat("NO MEDIATION\n")
    cat("→ Indirect effect not credibly different from 0\n")
  }

  cat("\nBayesian Fit Indices:\n")
  print(blavFitIndices(bsem_path))
  cat("\n")
}

################################################################################
# PART 6: MODEL COMPARISON SUMMARY
################################################################################

cat("================================================================================\n")
cat("PART 6: MODEL COMPARISON SUMMARY\n")
cat("================================================================================\n\n")

cat("SUMMARY OF ALL DAG MODELS:\n")
cat("--------------------------\n\n")

# Collect results if available
if (!is.null(bsem_dag1)) {
  dag1_params <- parameterEstimates(bsem_dag1)
  rai_eff <- dag1_params[dag1_params$lhs == "MIL" & dag1_params$rhs == "RAI", ]
  rai_e <- get_param_est(rai_eff, "est"); rai_l <- get_param_est(rai_eff, "ci.lower"); rai_h <- get_param_est(rai_eff, "ci.upper")
  cat("DAG 1 (Traditional SDT - RAI only):\n")
  cat("  RAI → MIL: β =", round(rai_e, 3), "[", round(rai_l, 3), ",", round(rai_h, 3), "]\n\n")
}

if (!is.null(bsem_dag2)) {
  dag2_params <- parameterEstimates(bsem_dag2)
  auto_eff <- dag2_params[dag2_params$lhs == "MIL" & dag2_params$rhs == "AUTO", ]
  auto_e <- get_param_est(auto_eff, "est"); auto_l <- get_param_est(auto_eff, "ci.lower"); auto_h <- get_param_est(auto_eff, "ci.upper")
  cat("DAG 2 (YOUR HYPOTHESIS - Autonomous only):\n")
  cat("  AUTO → MIL: β =", round(auto_e, 3), "[", round(auto_l, 3), ",", round(auto_h, 3), "]\n")
  cat("  CTRL → MIL: constrained to 0\n\n")
}

if (!is.null(bsem_dag3)) {
  dag3_params <- parameterEstimates(bsem_dag3)
  auto_eff <- dag3_params[dag3_params$lhs == "MIL" & dag3_params$rhs == "AUTO", ]
  ctrl_eff <- dag3_params[dag3_params$lhs == "MIL" & dag3_params$rhs == "CTRL", ]
  auto_e <- get_param_est(auto_eff, "est"); auto_l <- get_param_est(auto_eff, "ci.lower"); auto_h <- get_param_est(auto_eff, "ci.upper")
  ctrl_e <- get_param_est(ctrl_eff, "est"); ctrl_l <- get_param_est(ctrl_eff, "ci.lower"); ctrl_h <- get_param_est(ctrl_eff, "ci.upper")
  cat("DAG 3 (Dual Pathway - Both):\n")
  cat("  AUTO → MIL: β =", round(auto_e, 3), "[", round(auto_l, 3), ",", round(auto_h, 3), "]\n")
  cat("  CTRL → MIL: β =", round(ctrl_e, 3), "[", round(ctrl_l, 3), ",", round(ctrl_h, 3), "]")
  if (!is.na(ctrl_l) && !is.na(ctrl_h) && !(ctrl_l > 0 || ctrl_h < 0)) {
    cat(" ← INCLUDES 0!\n")
  } else {
    cat("\n")
  }
  cat("\n")
}

cat("CONCLUSION:\n")
cat("-----------\n")
cat("If CTRL → MIL includes 0 in DAG 3, the constraint in DAG 2 is justified.\n")
cat("If indirect effect is credible and direct is not, full mediation is supported.\n")
cat("Both findings support YOUR HYPOTHESIS (DAG 2).\n")

################################################################################
# SAVE RESULTS
################################################################################

cat("\n================================================================================\n")
cat("SAVING RESULTS\n")
cat("================================================================================\n\n")

results <- list(
  bcfa_full = bcfa_full,
  bsem_dag1 = bsem_dag1,
  bsem_dag2 = bsem_dag2,
  bsem_dag3 = bsem_dag3,
  bsem_path = bsem_path
)

saveRDS(results, file = "Results/Latent/Script11b_Bayesian_Results.rds")
cat("Results saved to: Results/Latent/Script11b_Bayesian_Results.rds\n")

cat("\n================================================================================\n")
cat("SCRIPT 11b COMPLETE - BAYESIAN LATENT SEM\n")
cat("================================================================================\n")
