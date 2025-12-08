################################################################################
# SCRIPT 13: IEA FRAMEWORK - DISTRESS MEDIATION MODEL
################################################################################
#
# Purpose: Test the Intention-Experience Alignment (IEA) Framework
#          The "missing link" - explicitly modeling CTRL → DEP → MIL
#
# Theoretical Framework:
#   - Meaning = experience tracking intended trajectory
#   - Depression = signal that experience is NOT tracking intended trajectory
#   - Autonomous motivation = goals that ARE your trajectory → meaning directly
#   - Controlled motivation = goals outside trajectory → misalignment → distress → meaning loss
#
# The Key Hypothesis:
#   - AUTO → MIL (direct, positive) = Alignment creates meaning
#   - CTRL → MIL (direct) = ZERO (no direct effect)
#   - CTRL → DEP → MIL (indirect, negative) = Misalignment creates distress signal
#
# This is the definitive test of the IEA framework.
#
################################################################################

library(haven)
library(dplyr)
library(lavaan)
library(blavaan)

cat("================================================================================\n")
cat("SCRIPT 13: IEA FRAMEWORK - DISTRESS MEDIATION MODEL\n")
cat("================================================================================\n")
cat("Testing: Controlled motivation destroys meaning ONLY via the distress signal\n")
cat("================================================================================\n\n")

# Set paths
base_dir <- getwd()
data_dir <- file.path(dirname(base_dir), "Data")

# Create output directory
if (!dir.exists("Results/Latent")) {
  dir.create("Results/Latent", recursive = TRUE)
}

# MCMC settings
n_chains <- 3
n_burnin <- 1000
n_sample <- 1000

cat("MCMC Settings: Chains =", n_chains, ", Burnin =", n_burnin, ", Samples =", n_sample, "\n\n")

################################################################################
# LOAD AND PREPARE DATA
################################################################################

cat("Loading Combined dataset...\n")
combined <- read_sav(file.path(data_dir, "Combined_Study1_Study2_FULL.sav"))
cat("  N =", nrow(combined), "\n\n")

################################################################################
# CREATE GOAL-LEVEL PARCELS
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

# Define CESD items
cesd_items <- c("CESD_1", "CESD_2", "CESD_3", "CESD_4", "CESD_5",
                "CESD_6", "CESD_7", "CESD_8R", "CESD_9", "CESD_10")
cesd_items <- cesd_items[cesd_items %in% names(combined)]

################################################################################
# THE IEA MODEL: Full Distress Mediation
################################################################################

cat("================================================================================\n")
cat("THE IEA MODEL: INTENTION-EXPERIENCE ALIGNMENT FRAMEWORK\n")
cat("================================================================================\n\n")

cat("Model Structure:\n")
cat("  MIL ~ a*AUTO + b*DEP + c_prime*CTRL    (Meaning predicted by all three)\n")
cat("  DEP ~ d*CTRL + e*AUTO                   (Depression predicted by both motivations)\n")
cat("\n")
cat("Key Paths to Test:\n")
cat("  - AUTO → MIL (direct): Should be POSITIVE (alignment creates meaning)\n")
cat("  - CTRL → MIL (direct): Should be ZERO (no direct effect)\n")
cat("  - CTRL → DEP (d): Should be POSITIVE (misalignment creates distress)\n")
cat("  - AUTO → DEP (e): Should be NEGATIVE (alignment protects from distress)\n")
cat("  - distress_path (d*b): CTRL → DEP → MIL - the indirect harm pathway\n")
cat("\n")

iea_model <- paste0('
  # ============================================
  # MEASUREMENT MODEL
  # ============================================

  # Meaning in Life (4 items)
  MIL =~ MIL1 + MIL2 + MIL3 + MIL4

  # Autonomous Motivation (6 goal-level parcels)
  AUTO =~ Goal_1_Auto + Goal_2_Auto + Goal_3_Auto + Goal_4_Auto + Goal_5_Auto + Goal_6_Auto

  # Controlled Motivation (6 goal-level parcels)
  CTRL =~ Goal_1_Ctrl + Goal_2_Ctrl + Goal_3_Ctrl + Goal_4_Ctrl + Goal_5_Ctrl + Goal_6_Ctrl

  # Depression/Distress (CESD items)
  DEP =~ ', paste(cesd_items, collapse = " + "), '

  # ============================================
  # STRUCTURAL MODEL (IEA Framework)
  # ============================================

  # Meaning predicted by: AUTO (direct), CTRL (direct), DEP
  MIL ~ a*AUTO + c_prime*CTRL + b*DEP

  # Depression predicted by: CTRL and AUTO
  DEP ~ d*CTRL + e*AUTO

  # ============================================
  # EFFECT DECOMPOSITION
  # ============================================

  # The distress pathway: CTRL → DEP → MIL
  distress_path := d * b

  # AUTO indirect pathway: AUTO → DEP → MIL
  auto_indirect := e * b

  # Total effects
  ctrl_total := c_prime + (d * b)
  auto_total := a + (e * b)
')

cat("Fitting IEA Model (Bayesian SEM)...\n")
cat("This may take a few minutes...\n\n")

bsem_iea <- tryCatch({
  bsem(iea_model, data = combined,
       n.chains = n_chains, burnin = n_burnin, sample = n_sample,
       std.lv = TRUE, seed = 42)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

################################################################################
# EXTRACT AND REPORT RESULTS
################################################################################

if (!is.null(bsem_iea)) {

  cat("\n================================================================================\n")
  cat("IEA MODEL RESULTS\n")
  cat("================================================================================\n\n")

  # Get parameter estimates
  params <- parameterEstimates(bsem_iea)

  # Determine column names (blavaan uses different names)
  est_col <- ifelse("Estimate" %in% names(params), "Estimate", "est")
  ci_lo_col <- ifelse("pi.lower" %in% names(params), "pi.lower", "ci.lower")
  ci_hi_col <- ifelse("pi.upper" %in% names(params), "pi.upper", "ci.upper")

  # Helper to extract and format
  get_result <- function(label_or_filter, is_label = TRUE) {
    if (is_label) {
      row <- params[params$label == label_or_filter, ]
    } else {
      row <- label_or_filter
    }
    if (nrow(row) == 0) return(list(est = NA, lo = NA, hi = NA, credible = NA))

    est <- row[[est_col]][1]
    lo <- row[[ci_lo_col]][1]
    hi <- row[[ci_hi_col]][1]
    credible <- !is.na(lo) && !is.na(hi) && (lo > 0 || hi < 0)

    return(list(est = est, lo = lo, hi = hi, credible = credible))
  }

  # Extract all key paths
  a <- get_result("a")           # AUTO → MIL (direct)
  b <- get_result("b")           # DEP → MIL
  c_prime <- get_result("c_prime")  # CTRL → MIL (direct)
  d <- get_result("d")           # CTRL → DEP
  e <- get_result("e")           # AUTO → DEP
  distress <- get_result("distress_path")
  auto_ind <- get_result("auto_indirect")
  ctrl_tot <- get_result("ctrl_total")
  auto_tot <- get_result("auto_total")

  # ============================================
  # PRINT RESULTS
  # ============================================

  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("STRUCTURAL PATHS\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

  cat("PATHS TO MEANING (MIL):\n")
  cat("------------------------\n")
  cat(sprintf("  AUTO → MIL (direct):  β = %6.3f [%6.3f, %6.3f]  %s\n",
              a$est, a$lo, a$hi,
              ifelse(a$credible, "✓ CREDIBLE (Alignment → Meaning)", "✗ Not credible")))
  cat(sprintf("  CTRL → MIL (direct):  β = %6.3f [%6.3f, %6.3f]  %s\n",
              c_prime$est, c_prime$lo, c_prime$hi,
              ifelse(!c_prime$credible, "✓ NOT CREDIBLE (No direct harm)", "✗ Credible")))
  cat(sprintf("  DEP → MIL:            β = %6.3f [%6.3f, %6.3f]  %s\n",
              b$est, b$lo, b$hi,
              ifelse(b$credible, "✓ CREDIBLE", "Not credible")))

  cat("\nPATHS TO DEPRESSION (DEP):\n")
  cat("---------------------------\n")
  cat(sprintf("  CTRL → DEP:           β = %6.3f [%6.3f, %6.3f]  %s\n",
              d$est, d$lo, d$hi,
              ifelse(d$credible & d$est > 0, "✓ CREDIBLE (Misalignment → Distress)",
                     ifelse(d$credible, "✓ CREDIBLE", "✗ Not credible"))))
  cat(sprintf("  AUTO → DEP:           β = %6.3f [%6.3f, %6.3f]  %s\n",
              e$est, e$lo, e$hi,
              ifelse(e$credible & e$est < 0, "✓ CREDIBLE (Alignment protects)",
                     ifelse(e$credible, "✓ CREDIBLE", "✗ Not credible"))))

  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("INDIRECT & TOTAL EFFECTS\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

  cat("THE DISTRESS PATHWAY (CTRL → DEP → MIL):\n")
  cat("-----------------------------------------\n")
  cat(sprintf("  Indirect effect:      β = %6.3f [%6.3f, %6.3f]  %s\n",
              distress$est, distress$lo, distress$hi,
              ifelse(distress$credible, "✓ CREDIBLE", "✗ Not credible")))

  cat("\nAUTONOMOUS INDIRECT (AUTO → DEP → MIL):\n")
  cat("----------------------------------------\n")
  cat(sprintf("  Indirect effect:      β = %6.3f [%6.3f, %6.3f]  %s\n",
              auto_ind$est, auto_ind$lo, auto_ind$hi,
              ifelse(auto_ind$credible, "✓ CREDIBLE (bonus via less distress)", "✗ Not credible")))

  cat("\nTOTAL EFFECTS ON MEANING:\n")
  cat("--------------------------\n")
  cat(sprintf("  AUTO total:           β = %6.3f [%6.3f, %6.3f]\n",
              auto_tot$est, auto_tot$lo, auto_tot$hi))
  cat(sprintf("  CTRL total:           β = %6.3f [%6.3f, %6.3f]\n",
              ctrl_tot$est, ctrl_tot$lo, ctrl_tot$hi))

  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("IEA FRAMEWORK EVALUATION\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

  # Evaluate framework support
  auto_direct_pos <- a$credible && a$est > 0
  ctrl_direct_zero <- !c_prime$credible
  ctrl_dep_pos <- d$credible && d$est > 0
  auto_dep_neg <- e$credible && e$est < 0
  distress_path_neg <- distress$credible && distress$est < 0

  tests_passed <- sum(c(auto_direct_pos, ctrl_direct_zero, ctrl_dep_pos, auto_dep_neg, distress_path_neg))

  cat("HYPOTHESIS TESTS:\n")
  cat("-----------------\n")
  cat(sprintf("  [%s] AUTO → MIL is positive (alignment creates meaning)\n",
              ifelse(auto_direct_pos, "✓", "✗")))
  cat(sprintf("  [%s] CTRL → MIL is zero (no direct harm)\n",
              ifelse(ctrl_direct_zero, "✓", "✗")))
  cat(sprintf("  [%s] CTRL → DEP is positive (misalignment creates distress)\n",
              ifelse(ctrl_dep_pos, "✓", "✗")))
  cat(sprintf("  [%s] AUTO → DEP is negative (alignment protects from distress)\n",
              ifelse(auto_dep_neg, "✓", "✗")))
  cat(sprintf("  [%s] CTRL → DEP → MIL is negative (the distress pathway)\n",
              ifelse(distress_path_neg, "✓", "✗")))

  cat(sprintf("\nTESTS PASSED: %d/5\n", tests_passed))

  if (tests_passed == 5) {
    cat("\n╔══════════════════════════════════════════════════════════════════════════════╗\n")
    cat("║                    IEA FRAMEWORK FULLY SUPPORTED                             ║\n")
    cat("╠══════════════════════════════════════════════════════════════════════════════╣\n")
    cat("║  • Autonomous motivation creates meaning DIRECTLY (intention-experience     ║\n")
    cat("║    alignment)                                                               ║\n")
    cat("║  • Controlled motivation has NO direct effect on meaning                    ║\n")
    cat("║  • Controlled motivation harms meaning ONLY through the distress signal     ║\n")
    cat("║  • Depression is the felt consequence of experiencing life out of           ║\n")
    cat("║    alignment with your overall meta-trajectory                              ║\n")
    cat("╚══════════════════════════════════════════════════════════════════════════════╝\n")
  } else if (tests_passed >= 3) {
    cat("\n*** IEA FRAMEWORK PARTIALLY SUPPORTED ***\n")
  } else {
    cat("\n*** IEA FRAMEWORK NOT SUPPORTED ***\n")
  }

  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("COMPARISON: AUTO vs CTRL EFFECTS\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

  cat("AUTONOMOUS MOTIVATION:\n")
  cat(sprintf("  Direct effect:   β = %+.3f (%5.1f%% of total)\n",
              a$est, abs(a$est / auto_tot$est) * 100))
  cat(sprintf("  Indirect effect: β = %+.3f (%5.1f%% of total)\n",
              auto_ind$est, abs(auto_ind$est / auto_tot$est) * 100))
  cat(sprintf("  TOTAL:           β = %+.3f\n\n", auto_tot$est))

  cat("CONTROLLED MOTIVATION:\n")
  cat(sprintf("  Direct effect:   β = %+.3f (%5.1f%% of total)\n",
              c_prime$est, ifelse(ctrl_tot$est != 0, abs(c_prime$est / ctrl_tot$est) * 100, 0)))
  cat(sprintf("  Indirect effect: β = %+.3f (%5.1f%% of total)  ← THE DISTRESS PATHWAY\n",
              distress$est, ifelse(ctrl_tot$est != 0, abs(distress$est / ctrl_tot$est) * 100, 100)))
  cat(sprintf("  TOTAL:           β = %+.3f\n\n", ctrl_tot$est))

  cat(sprintf("AUTO total effect is %.1fx stronger than CTRL total effect\n",
              abs(auto_tot$est) / abs(ctrl_tot$est)))

  # Fit indices
  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("MODEL FIT\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
  fit <- blavFitIndices(bsem_iea)
  print(fit)

  # Save results
  cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("SAVING RESULTS\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

  results <- list(
    model = bsem_iea,
    paths = list(
      auto_direct = a,
      ctrl_direct = c_prime,
      dep_to_mil = b,
      ctrl_to_dep = d,
      auto_to_dep = e,
      distress_path = distress,
      auto_indirect = auto_ind,
      ctrl_total = ctrl_tot,
      auto_total = auto_tot
    ),
    framework_supported = tests_passed == 5
  )

  saveRDS(results, file = "Results/Latent/Script13_IEA_Results.rds")
  cat("Results saved to: Results/Latent/Script13_IEA_Results.rds\n")

} else {
  cat("\nModel failed to converge. Check error messages above.\n")
}

cat("\n================================================================================\n")
cat("SCRIPT 13 COMPLETE - IEA FRAMEWORK TEST\n")
cat("================================================================================\n")
