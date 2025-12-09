################################################################################
# SCRIPT 11: LATENT CONSTRUCT MODELING
################################################################################
#
# Purpose: Supplement the observed-variable analyses (Scripts 02-10) with
#          latent variable models that account for measurement error
#
# Structure:
#   PART 1: Confirmatory Factor Analysis (CFA) - Measurement Models
#   PART 2: Latent Structural Models - DAG Comparison
#   PART 3: Latent Path Decomposition Analysis
#
# Theoretical Framework:
#   - Judea Pearl's causal inference (do-calculus)
#   - Self-Determination Theory (SDT) motivation continuum
#   - Your DAG 2 hypothesis: Autonomous → Meaning (Controlled = 0)
#
# GSC Item Structure (verified from data):
#   GSC#_1 = External regulation (CONTROLLED) - "somebody else wants you to"
#   GSC#_2 = Identified regulation (AUTONOMOUS) - "you really believe it's important"
#   GSC#_3 = Introjected regulation (CONTROLLED) - "feel ashamed, guilty, anxious"
#   GSC#_4 = Intrinsic regulation (AUTONOMOUS) - "enjoyment or stimulation"
#
################################################################################

# library(haven) # Not needed for CSV
library(dplyr)
library(lavaan)
# library(semPlot)  # Optional - for path diagrams

cat("================================================================================\n")
cat("SCRIPT 11: LATENT CONSTRUCT MODELING\n")
cat("================================================================================\n\n")

# Set paths
base_dir <- "."
data_dir <- "."

# Create output directory for latent models
if (!dir.exists("Results/Latent")) {
  dir.create("Results/Latent", recursive = TRUE)
}

################################################################################
# LOAD AND PREPARE DATA
################################################################################

cat("Loading datasets...\n")

study1 <- read.csv("dataset.csv")
study2 <- read.csv("dataset.csv")
combined <- read.csv("dataset.csv")  # With GSC items

cat("  Study 1: N =", nrow(study1), "\n")
cat("  Study 2: N =", nrow(study2), "\n")
cat("  Combined: N =", nrow(combined), "\n\n")

# Define item names
mil_items <- c("MIL1", "MIL2", "MIL3", "MIL4")
cesd_items <- paste0("CESD_", 1:10)
gad_items <- paste0("GAD_", 1:7)

# Autonomous items (Identified + Intrinsic for each of 6 goals)
auto_items <- c(paste0("GSC", 1:6, "_2"), paste0("GSC", 1:6, "_4"))  # _2 and _4

# Controlled items (External + Introjected for each of 6 goals)
ctrl_items <- c(paste0("GSC", 1:6, "_1"), paste0("GSC", 1:6, "_3"))  # _1 and _3

################################################################################
# PRE-PROCESSING: CREATE GOAL-LEVEL PARCELS
################################################################################
# Item parcelling approach - standard in SDT literature to improve model fit
# by reducing item-level noise and method effects across 6 distinct goals
#
# For each goal i (1-6):
#   Goal_i_Auto = Mean(GSCi_2, GSCi_4)  # Identified + Intrinsic
#   Goal_i_Ctrl = Mean(GSCi_1, GSCi_3)  # External + Introjected
################################################################################

cat("Creating goal-level parcels for motivation...\n\n")

create_parcels <- function(data) {
  for (i in 1:6) {
    # Autonomous parcel: mean of Identified (_2) and Intrinsic (_4)
    auto_col <- paste0("Goal_", i, "_Auto")
    data[[auto_col]] <- rowMeans(data[, c(paste0("GSC", i, "_2"), paste0("GSC", i, "_4"))], na.rm = TRUE)

    # Controlled parcel: mean of External (_1) and Introjected (_3)
    ctrl_col <- paste0("Goal_", i, "_Ctrl")
    data[[ctrl_col]] <- rowMeans(data[, c(paste0("GSC", i, "_1"), paste0("GSC", i, "_3"))], na.rm = TRUE)
  }
  return(data)
}

# Apply parcelling to all datasets
study1 <- create_parcels(study1)
study2 <- create_parcels(study2)
combined <- create_parcels(combined)

# Define parcel names for latent models
auto_parcels <- paste0("Goal_", 1:6, "_Auto")
ctrl_parcels <- paste0("Goal_", 1:6, "_Ctrl")

cat("  Created 6 Autonomous parcels:", paste(auto_parcels, collapse = ", "), "\n")
cat("  Created 6 Controlled parcels:", paste(ctrl_parcels, collapse = ", "), "\n\n")

################################################################################
# PART 1: CONFIRMATORY FACTOR ANALYSIS (CFA)
################################################################################

cat("================================================================================\n")
cat("PART 1: CONFIRMATORY FACTOR ANALYSIS - MEASUREMENT MODELS\n")
cat("================================================================================\n\n")

run_cfa <- function(data, dataset_name) {

  cat("--- CFA for", dataset_name, "---\n\n")

  # Check which variables exist in this dataset
  has_gsc_items <- all(c("GSC1_1", "GSC1_2", "GSC1_3", "GSC1_4") %in% names(data))

  if (!has_gsc_items) {
    cat("  GSC items not available in", dataset_name, "- using composite scores\n\n")
    return(NULL)
  }

  # ============================================================================
  # MODEL 1: Single-factor MIL
  # ============================================================================
  cat("1. Meaning in Life (MIL) - Single Factor\n")

  mil_model <- '
    MIL =~ MIL1 + MIL2 + MIL3 + MIL4
  '

  mil_fit <- cfa(mil_model, data = data, estimator = "MLR", std.lv = TRUE)
  mil_summary <- fitMeasures(mil_fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

  cat("   Chi-square:", round(mil_summary["chisq"], 2), "(df =", mil_summary["df"], ", p =", round(mil_summary["pvalue"], 3), ")\n")
  cat("   CFI:", round(mil_summary["cfi"], 3), "\n")
  cat("   TLI:", round(mil_summary["tli"], 3), "\n")
  cat("   RMSEA:", round(mil_summary["rmsea"], 3), "\n")
  cat("   SRMR:", round(mil_summary["srmr"], 3), "\n")

  # Check fit
  if (mil_summary["cfi"] >= 0.95 & mil_summary["rmsea"] <= 0.08) {
    cat("   Fit: GOOD\n\n")
  } else if (mil_summary["cfi"] >= 0.90) {
    cat("   Fit: ACCEPTABLE\n\n")
  } else {
    cat("   Fit: POOR\n\n")
  }

  # ============================================================================
  # MODEL 2: Two-factor Motivation using GOAL-LEVEL PARCELS
  # ============================================================================
  cat("2. Motivation - Two Factors using Goal-Level Parcels\n")
  cat("   (6 parcels per factor instead of 12 raw items)\n")

  # Use parcels instead of raw items
  motivation_model <- '
    # Autonomous Motivation (6 goal-level parcels)
    AUTO =~ Goal_1_Auto + Goal_2_Auto + Goal_3_Auto + Goal_4_Auto + Goal_5_Auto + Goal_6_Auto

    # Controlled Motivation (6 goal-level parcels)
    CTRL =~ Goal_1_Ctrl + Goal_2_Ctrl + Goal_3_Ctrl + Goal_4_Ctrl + Goal_5_Ctrl + Goal_6_Ctrl
  '

  motivation_fit <- tryCatch({
    cfa(motivation_model, data = data, estimator = "MLR", std.lv = TRUE)
  }, error = function(e) {
    cat("   Error fitting motivation CFA:", e$message, "\n\n")
    return(NULL)
  })

  if (!is.null(motivation_fit)) {
    mot_summary <- fitMeasures(motivation_fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

    cat("   Chi-square:", round(mot_summary["chisq"], 2), "(df =", mot_summary["df"], ")\n")
    cat("   CFI:", round(mot_summary["cfi"], 3), "\n")
    cat("   TLI:", round(mot_summary["tli"], 3), "\n")
    cat("   RMSEA:", round(mot_summary["rmsea"], 3), "\n")
    cat("   SRMR:", round(mot_summary["srmr"], 3), "\n")

    # Check fit
    if (mot_summary["cfi"] >= 0.95 & mot_summary["rmsea"] <= 0.08) {
      cat("   Fit: GOOD\n")
    } else if (mot_summary["cfi"] >= 0.90) {
      cat("   Fit: ACCEPTABLE\n")
    } else {
      cat("   Fit: POOR\n")
    }

    # Correlation between factors
    mot_params <- parameterEstimates(motivation_fit, standardized = TRUE)
    auto_ctrl_cor <- mot_params[mot_params$lhs == "AUTO" & mot_params$rhs == "CTRL", "std.all"]
    if (length(auto_ctrl_cor) > 0) {
      cat("   AUTO-CTRL correlation:", round(auto_ctrl_cor, 3), "\n")
    }
    cat("\n")
  }

  # ============================================================================
  # MODEL 3: Depression (CESD-10)
  # ============================================================================
  cat("3. Depression (CESD-10) - Single Factor\n")

  cesd_items_present <- cesd_items[cesd_items %in% names(data)]

  # Use CESD_8R if available (reverse-coded)
  if ("CESD_8R" %in% names(data)) {
    cesd_items_present <- c(cesd_items_present[cesd_items_present != "CESD_8"], "CESD_8R")
  }

  dep_model <- paste0('
    DEP =~ ', paste(cesd_items_present, collapse = " + "), '
  ')

  dep_fit <- tryCatch({
    cfa(dep_model, data = data, estimator = "MLR", std.lv = TRUE)
  }, error = function(e) {
    cat("   Error fitting depression CFA:", e$message, "\n\n")
    return(NULL)
  })

  if (!is.null(dep_fit)) {
    dep_summary <- fitMeasures(dep_fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

    cat("   Chi-square:", round(dep_summary["chisq"], 2), "(df =", dep_summary["df"], ")\n")
    cat("   CFI:", round(dep_summary["cfi"], 3), "\n")
    cat("   TLI:", round(dep_summary["tli"], 3), "\n")
    cat("   RMSEA:", round(dep_summary["rmsea"], 3), "\n")
    cat("   SRMR:", round(dep_summary["srmr"], 3), "\n\n")
  }

  # ============================================================================
  # MODEL 4: Full Measurement Model using PARCELS
  # ============================================================================
  cat("4. Full Measurement Model (with Goal-Level Parcels)\n")

  full_measurement_model <- paste0('
    # Meaning in Life
    MIL =~ MIL1 + MIL2 + MIL3 + MIL4

    # Autonomous Motivation (using parcels)
    AUTO =~ Goal_1_Auto + Goal_2_Auto + Goal_3_Auto + Goal_4_Auto + Goal_5_Auto + Goal_6_Auto

    # Controlled Motivation (using parcels)
    CTRL =~ Goal_1_Ctrl + Goal_2_Ctrl + Goal_3_Ctrl + Goal_4_Ctrl + Goal_5_Ctrl + Goal_6_Ctrl

    # Depression
    DEP =~ ', paste(cesd_items_present, collapse = " + "), '
  ')

  full_fit <- tryCatch({
    cfa(full_measurement_model, data = data, estimator = "MLR", std.lv = TRUE)
  }, error = function(e) {
    cat("   Error fitting full measurement model:", e$message, "\n\n")
    return(NULL)
  })

  if (!is.null(full_fit)) {
    full_summary <- fitMeasures(full_fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

    cat("   Chi-square:", round(full_summary["chisq"], 2), "(df =", full_summary["df"], ")\n")
    cat("   CFI:", round(full_summary["cfi"], 3), "\n")
    cat("   TLI:", round(full_summary["tli"], 3), "\n")
    cat("   RMSEA:", round(full_summary["rmsea"], 3), "\n")
    cat("   SRMR:", round(full_summary["srmr"], 3), "\n")

    # Check fit
    if (full_summary["cfi"] >= 0.95 & full_summary["rmsea"] <= 0.08) {
      cat("   Fit: GOOD\n\n")
    } else if (full_summary["cfi"] >= 0.90) {
      cat("   Fit: ACCEPTABLE\n\n")
    } else {
      cat("   Fit: POOR\n\n")
    }

    # Factor correlations
    cat("   Factor Correlations:\n")
    full_params <- parameterEstimates(full_fit, standardized = TRUE)
    cors <- full_params[full_params$op == "~~" & full_params$lhs != full_params$rhs,
                        c("lhs", "rhs", "std.all")]
    for (i in 1:nrow(cors)) {
      cat("     ", cors$lhs[i], "<->", cors$rhs[i], ":", round(cors$std.all[i], 3), "\n")
    }
    cat("\n")
  }

  # Return results
  return(list(
    dataset = dataset_name,
    mil_fit = mil_fit,
    motivation_fit = motivation_fit,
    depression_fit = dep_fit,
    full_fit = full_fit
  ))
}

# Run CFA for Combined dataset
cfa_combined <- run_cfa(combined, "Combined")


################################################################################
# PART 2: LATENT STRUCTURAL MODELS - DAG COMPARISON
################################################################################

cat("\n================================================================================\n")
cat("PART 2: LATENT STRUCTURAL MODELS - DAG COMPARISON\n")
cat("================================================================================\n\n")

run_latent_dag <- function(data, dataset_name) {

  cat("--- Latent DAG Comparison for", dataset_name, "---\n\n")

  # Check for required variables
  has_gsc_items <- all(c("GSC1_1", "GSC1_2", "GSC1_3", "GSC1_4") %in% names(data))

  if (!has_gsc_items) {
    cat("  GSC items not available - skipping latent DAG analysis\n\n")
    return(NULL)
  }

  # Get available items
  auto_items_present <- auto_items[auto_items %in% names(data)]
  ctrl_items_present <- ctrl_items[ctrl_items %in% names(data)]
  cesd_items_present <- cesd_items[cesd_items %in% names(data)]
  if ("CESD_8R" %in% names(data)) {
    cesd_items_present <- c(cesd_items_present[cesd_items_present != "CESD_8"], "CESD_8R")
  }

  # ============================================================================
  # DAG 1: Traditional SDT (RAI → Meaning)
  # Uses observed RAI since it's a computed balance score
  # ============================================================================
  cat("DAG 1: Traditional SDT - RAI (Balance) → Meaning\n")

  dag1_model <- paste0('
    # Measurement Model
    MIL =~ MIL1 + MIL2 + MIL3 + MIL4
    DEP =~ ', paste(cesd_items_present, collapse = " + "), '

    # Structural Model
    MIL ~ RAI + DEP + Age + Sex
  ')

  dag1_fit <- tryCatch({
    sem(dag1_model, data = data, estimator = "MLR")
  }, error = function(e) {
    cat("   Error:", e$message, "\n")
    return(NULL)
  })

  if (!is.null(dag1_fit)) {
    dag1_fitmeasures <- fitMeasures(dag1_fit, c("chisq", "df", "cfi", "rmsea", "aic", "bic"))
    dag1_params <- parameterEstimates(dag1_fit, standardized = TRUE)

    rai_effect <- dag1_params[dag1_params$lhs == "MIL" & dag1_params$rhs == "RAI", ]

    cat("   RAI → MIL: β =", round(rai_effect$std.all, 3),
        "[", round(rai_effect$ci.lower, 3), ",", round(rai_effect$ci.upper, 3), "]\n")
    cat("   AIC:", round(dag1_fitmeasures["aic"], 1), "\n")
    cat("   BIC:", round(dag1_fitmeasures["bic"], 1), "\n\n")
  }

  # ============================================================================
  # DAG 2: Your Hypothesis - Autonomous → Meaning (Controlled = 0)
  # Using GOAL-LEVEL PARCELS
  # ============================================================================
  cat("DAG 2: Autonomous-Primary (YOUR HYPOTHESIS)\n")
  cat("   Controlled motivation constrained to 0 (not in model)\n")
  cat("   Using goal-level parcels for motivation\n")

  dag2_model <- paste0('
    # Measurement Model (with parcels)
    MIL =~ MIL1 + MIL2 + MIL3 + MIL4
    AUTO =~ Goal_1_Auto + Goal_2_Auto + Goal_3_Auto + Goal_4_Auto + Goal_5_Auto + Goal_6_Auto
    DEP =~ ', paste(cesd_items_present, collapse = " + "), '

    # Structural Model - Controlled OMITTED (constrained to 0)
    MIL ~ AUTO + DEP + Age + Sex
  ')

  dag2_fit <- tryCatch({
    sem(dag2_model, data = data, estimator = "MLR")
  }, error = function(e) {
    cat("   Error:", e$message, "\n")
    return(NULL)
  })

  if (!is.null(dag2_fit)) {
    dag2_fitmeasures <- fitMeasures(dag2_fit, c("chisq", "df", "cfi", "rmsea", "aic", "bic"))
    dag2_params <- parameterEstimates(dag2_fit, standardized = TRUE)

    auto_effect <- dag2_params[dag2_params$lhs == "MIL" & dag2_params$rhs == "AUTO", ]

    cat("   AUTO → MIL: β =", round(auto_effect$std.all, 3),
        "[", round(auto_effect$ci.lower, 3), ",", round(auto_effect$ci.upper, 3), "]\n")
    cat("   CFI:", round(dag2_fitmeasures["cfi"], 3), "\n")
    cat("   RMSEA:", round(dag2_fitmeasures["rmsea"], 3), "\n")
    cat("   AIC:", round(dag2_fitmeasures["aic"], 1), "\n")
    cat("   BIC:", round(dag2_fitmeasures["bic"], 1), "\n\n")
  }

  # ============================================================================
  # DAG 3: Both Pathways - Autonomous + Controlled → Meaning
  # Using GOAL-LEVEL PARCELS
  # ============================================================================
  cat("DAG 3: Dual Pathway (Both AUTO and CTRL)\n")
  cat("   Using goal-level parcels for motivation\n")

  dag3_model <- paste0('
    # Measurement Model (with parcels)
    MIL =~ MIL1 + MIL2 + MIL3 + MIL4
    AUTO =~ Goal_1_Auto + Goal_2_Auto + Goal_3_Auto + Goal_4_Auto + Goal_5_Auto + Goal_6_Auto
    CTRL =~ Goal_1_Ctrl + Goal_2_Ctrl + Goal_3_Ctrl + Goal_4_Ctrl + Goal_5_Ctrl + Goal_6_Ctrl
    DEP =~ ', paste(cesd_items_present, collapse = " + "), '

    # Structural Model - Both pathways estimated
    MIL ~ AUTO + CTRL + DEP + Age + Sex
  ')

  dag3_fit <- tryCatch({
    sem(dag3_model, data = data, estimator = "MLR")
  }, error = function(e) {
    cat("   Error:", e$message, "\n")
    return(NULL)
  })

  if (!is.null(dag3_fit)) {
    dag3_fitmeasures <- fitMeasures(dag3_fit, c("chisq", "df", "cfi", "rmsea", "aic", "bic"))
    dag3_params <- parameterEstimates(dag3_fit, standardized = TRUE)

    auto_effect <- dag3_params[dag3_params$lhs == "MIL" & dag3_params$rhs == "AUTO", ]
    ctrl_effect <- dag3_params[dag3_params$lhs == "MIL" & dag3_params$rhs == "CTRL", ]

    cat("   AUTO → MIL: β =", round(auto_effect$std.all, 3),
        "[", round(auto_effect$ci.lower, 3), ",", round(auto_effect$ci.upper, 3), "]\n")
    cat("   CTRL → MIL: β =", round(ctrl_effect$std.all, 3),
        "[", round(ctrl_effect$ci.lower, 3), ",", round(ctrl_effect$ci.upper, 3), "]\n")
    cat("   CFI:", round(dag3_fitmeasures["cfi"], 3), "\n")
    cat("   RMSEA:", round(dag3_fitmeasures["rmsea"], 3), "\n")
    cat("   AIC:", round(dag3_fitmeasures["aic"], 1), "\n")
    cat("   BIC:", round(dag3_fitmeasures["bic"], 1), "\n\n")
  }

  # ============================================================================
  # MODEL COMPARISON
  # ============================================================================
  cat("MODEL COMPARISON:\n")
  cat("-----------------\n")

  if (!is.null(dag1_fit) && !is.null(dag2_fit) && !is.null(dag3_fit)) {
    aics <- c(DAG1 = dag1_fitmeasures["aic"],
              DAG2 = dag2_fitmeasures["aic"],
              DAG3 = dag3_fitmeasures["aic"])
    bics <- c(DAG1 = dag1_fitmeasures["bic"],
              DAG2 = dag2_fitmeasures["bic"],
              DAG3 = dag3_fitmeasures["bic"])

    cat("   AIC: DAG1 =", round(aics[1], 1),
        ", DAG2 =", round(aics[2], 1),
        ", DAG3 =", round(aics[3], 1), "\n")
    cat("   BIC: DAG1 =", round(bics[1], 1),
        ", DAG2 =", round(bics[2], 1),
        ", DAG3 =", round(bics[3], 1), "\n")

    best_aic <- names(which.min(aics))
    best_bic <- names(which.min(bics))

    cat("\n   Best by AIC:", best_aic, "\n")
    cat("   Best by BIC:", best_bic, "\n")

    if (best_aic == "DAG2" || best_bic == "DAG2") {
      cat("\n   → YOUR HYPOTHESIS (DAG 2) SUPPORTED!\n")
    }
    cat("\n")
  }

  return(list(
    dataset = dataset_name,
    dag1 = dag1_fit,
    dag2 = dag2_fit,
    dag3 = dag3_fit
  ))
}

# Run latent DAG comparison for Combined dataset
dag_combined <- run_latent_dag(combined, "Combined")


################################################################################
# PART 3: LATENT PATH DECOMPOSITION ANALYSIS
################################################################################

cat("\n================================================================================\n")
cat("PART 3: LATENT PATH DECOMPOSITION ANALYSIS\n")
cat("================================================================================\n")
cat("Testing: Does Autonomous mediate the RAI → Meaning relationship?\n\n")

run_latent_path_decomposition <- function(data, dataset_name) {

  cat("--- Path Decomposition for", dataset_name, "---\n\n")

  # Check for required variables
  has_gsc_items <- all(c("GSC1_1", "GSC1_2", "GSC1_3", "GSC1_4") %in% names(data))

  if (!has_gsc_items) {
    cat("  GSC items not available - skipping\n\n")
    return(NULL)
  }

  # Get available items for Depression
  cesd_items_present <- cesd_items[cesd_items %in% names(data)]
  if ("CESD_8R" %in% names(data)) {
    cesd_items_present <- c(cesd_items_present[cesd_items_present != "CESD_8"], "CESD_8R")
  }

  cat("Using goal-level parcels for Autonomous motivation\n\n")

  path_model <- paste0('
    # ================================================
    # MEASUREMENT MODEL (with Goal-Level Parcels)
    # ================================================

    # Meaning in Life (latent)
    MIL =~ MIL1 + MIL2 + MIL3 + MIL4

    # Autonomous Motivation (latent, using parcels)
    AUTO =~ Goal_1_Auto + Goal_2_Auto + Goal_3_Auto + Goal_4_Auto + Goal_5_Auto + Goal_6_Auto

    # Depression (latent)
    DEP =~ ', paste(cesd_items_present, collapse = " + "), '

    # ================================================
    # STRUCTURAL MODEL (Causal Paths)
    # ================================================

    # Path a: RAI → Autonomous
    AUTO ~ a*RAI + DEP + Age + Sex

    # Path b: Autonomous → Meaning (controlling for RAI)
    # Path c_prime: Direct effect of RAI on Meaning
    MIL ~ b*AUTO + c_prime*RAI + DEP + Age + Sex

    # ================================================
    # CAUSAL EFFECT DECOMPOSITION
    # ================================================

    # Indirect effect (mediated through Autonomous)
    indirect := a * b

    # Direct effect (RAI → Meaning, not through Autonomous)
    direct := c_prime

    # Total effect
    total := c_prime + (a * b)

    # Proportion mediated
    prop_mediated := (a * b) / (c_prime + (a * b))
  ')

  path_fit <- tryCatch({
    sem(path_model, data = data, estimator = "MLR", se = "bootstrap", bootstrap = 1000)
  }, warning = function(w) {
    cat("   Warning:", w$message, "\n")
    sem(path_model, data = data, estimator = "MLR")
  }, error = function(e) {
    cat("   Error:", e$message, "\n")
    cat("   Trying without bootstrap...\n")
    tryCatch({
      sem(path_model, data = data, estimator = "MLR")
    }, error = function(e2) {
      cat("   Still failed:", e2$message, "\n")
      return(NULL)
    })
  })

  if (!is.null(path_fit)) {
    path_params <- parameterEstimates(path_fit, standardized = TRUE)

    # Extract key effects
    a_path <- path_params[path_params$label == "a", ]
    b_path <- path_params[path_params$label == "b", ]
    indirect <- path_params[path_params$label == "indirect", ]
    direct <- path_params[path_params$label == "direct", ]
    total <- path_params[path_params$label == "total", ]
    prop_med <- path_params[path_params$label == "prop_mediated", ]

    cat("CAUSAL PATH ESTIMATES:\n")
    cat("----------------------\n")
    cat("Path a (RAI → AUTO):     ", round(a_path$est, 3),
        " [", round(a_path$ci.lower, 3), ", ", round(a_path$ci.upper, 3), "]\n", sep = "")
    cat("Path b (AUTO → MIL):     ", round(b_path$est, 3),
        " [", round(b_path$ci.lower, 3), ", ", round(b_path$ci.upper, 3), "]\n", sep = "")
    cat("\n")

    cat("EFFECT DECOMPOSITION:\n")
    cat("---------------------\n")
    cat("Indirect (a × b):        ", round(indirect$est, 3),
        " [", round(indirect$ci.lower, 3), ", ", round(indirect$ci.upper, 3), "]\n", sep = "")
    cat("Direct (c'):             ", round(direct$est, 3),
        " [", round(direct$ci.lower, 3), ", ", round(direct$ci.upper, 3), "]\n", sep = "")
    cat("Total:                   ", round(total$est, 3),
        " [", round(total$ci.lower, 3), ", ", round(total$ci.upper, 3), "]\n", sep = "")
    cat("\n")

    cat("Proportion Mediated:     ", round(prop_med$est * 100, 1), "%\n", sep = "")

    # Interpretation
    cat("\nINTERPRETATION:\n")
    cat("---------------\n")

    if (indirect$ci.lower > 0 || indirect$ci.upper < 0) {
      cat("Indirect effect is SIGNIFICANT\n")
      cat("→ Autonomous DOES mediate the RAI → Meaning relationship\n")
    } else {
      cat("Indirect effect includes 0 (not significant)\n")
    }

    if (direct$ci.lower <= 0 && direct$ci.upper >= 0) {
      cat("→ Direct effect of RAI is NOT significant when controlling for Autonomous\n")
      cat("→ This supports FULL MEDIATION (and your DAG 2 hypothesis)\n")
    } else {
      cat("→ Direct effect of RAI remains significant\n")
      cat("→ This suggests PARTIAL MEDIATION\n")
    }
    cat("\n")

    # Model fit
    fit_measures <- fitMeasures(path_fit, c("chisq", "df", "cfi", "tli", "rmsea", "srmr"))
    cat("MODEL FIT:\n")
    cat("----------\n")
    cat("CFI:", round(fit_measures["cfi"], 3), "\n")
    cat("TLI:", round(fit_measures["tli"], 3), "\n")
    cat("RMSEA:", round(fit_measures["rmsea"], 3), "\n")
    cat("SRMR:", round(fit_measures["srmr"], 3), "\n\n")
  }

  return(list(
    dataset = dataset_name,
    fit = path_fit
  ))
}

# Run path decomposition for Combined dataset
path_combined <- run_latent_path_decomposition(combined, "Combined")


################################################################################
# SAVE ALL RESULTS
################################################################################

cat("\n================================================================================\n")
cat("SAVING RESULTS\n")
cat("================================================================================\n\n")

results <- list(
  cfa = list(combined = cfa_combined),
  dag = list(combined = dag_combined),
  path = list(combined = path_combined)
)

saveRDS(results, file = "Results/Latent/Script11_Results.rds")
cat("Results saved to: Results/Latent/Script11_Results.rds\n")

cat("\n================================================================================\n")
cat("SCRIPT 11 COMPLETE\n")
cat("================================================================================\n")
