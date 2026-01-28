# =============================================================================
# Script 15: Interaction Analysis (Autonomous x Controlled)
# =============================================================================


cat("\n")
cat("=============================================================================\n")
cat("SCRIPT 15: INTERACTION ANALYSIS (AUTONOMOUS x CONTROLLED)\n")
cat("=============================================================================\n\n")

# Load combined dataset
data <- read.csv("dataset.csv")

# Prepare analysis data
model_data <- data[complete.cases(data$Autonomous, data$MILjudgements,
                                   data$DEP_Corrected, data$Age, data$Sex,
                                   data$Controlled), ]

# Rename variables for clarity
model_data$Meaning <- model_data$MILjudgements
model_data$Depression <- model_data$DEP_Corrected

cat("Sample size: N =", nrow(model_data), "\n\n")

# =============================================================================
# 2. INTERACTION ANALYSIS
# =============================================================================

cat("=============================================================================\n")
cat("INTERACTION: AUTONOMOUS x CONTROLLED PREDICTING MEANING\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Model 1: Without Depression (covariates: Age, Sex)
# -----------------------------------------------------------------------------

cat("MODEL 1: Without Depression as Covariate\n")
cat("Formula: Meaning ~ Autonomous * Controlled + Age + Sex\n")
cat("-----------------------------------------------------------------------------\n")

m1 <- lm(Meaning ~ Autonomous * Controlled + Age + Sex, data = model_data)
s1 <- summary(m1)

# Print full model summary
cat("\nFull Model Summary:\n")
print(s1)

# Extract interaction term
int1_b <- coef(m1)["Autonomous:Controlled"]
int1_se <- s1$coefficients["Autonomous:Controlled", "Std. Error"]
int1_t <- s1$coefficients["Autonomous:Controlled", "t value"]
int1_p <- s1$coefficients["Autonomous:Controlled", "Pr(>|t|)"]

cat("\n>>> INTERACTION TERM (Autonomous:Controlled):\n")
cat("    b =", round(int1_b, 4), "\n")
cat("    SE =", round(int1_se, 4), "\n")
cat("    t =", round(int1_t, 3), "\n")
cat("    p =", round(int1_p, 3), "\n")
cat("    Significant:", ifelse(int1_p < .05, "YES", "NO"), "\n\n")

# -----------------------------------------------------------------------------
# Model 2: With Depression (covariates: Depression, Age, Sex)
# -----------------------------------------------------------------------------

cat("\nMODEL 2: With Depression as Covariate\n")
cat("Formula: Meaning ~ Autonomous * Controlled + Depression + Age + Sex\n")
cat("-----------------------------------------------------------------------------\n")

m2 <- lm(Meaning ~ Autonomous * Controlled + Depression + Age + Sex, data = model_data)
s2 <- summary(m2)

# Print full model summary
cat("\nFull Model Summary:\n")
print(s2)

# Extract interaction term
int2_b <- coef(m2)["Autonomous:Controlled"]
int2_se <- s2$coefficients["Autonomous:Controlled", "Std. Error"]
int2_t <- s2$coefficients["Autonomous:Controlled", "t value"]
int2_p <- s2$coefficients["Autonomous:Controlled", "Pr(>|t|)"]

cat("\n>>> INTERACTION TERM (Autonomous:Controlled):\n")
cat("    b =", round(int2_b, 4), "\n")
cat("    SE =", round(int2_se, 4), "\n")
cat("    t =", round(int2_t, 3), "\n")
cat("    p =", round(int2_p, 3), "\n")
cat("    Significant:", ifelse(int2_p < .05, "YES", "NO"), "\n\n")

# =============================================================================
# 3. SUMMARY
# =============================================================================

cat("=============================================================================\n")
cat("SUMMARY OF INTERACTION RESULTS\n")
cat("=============================================================================\n\n")

cat("Model                          | b       | SE      | t      | p      | Sig?\n")
cat("-------------------------------|---------|---------|--------|--------|------\n")
cat(sprintf("Without Depression             | %.4f  | %.4f  | %.3f  | %.3f  | %s\n",
            int1_b, int1_se, int1_t, int1_p, ifelse(int1_p < .05, "YES", "NO")))
cat(sprintf("With Depression                | %.4f  | %.4f  | %.3f  | %.3f  | %s\n",
            int2_b, int2_se, int2_t, int2_p, ifelse(int2_p < .05, "YES", "NO")))

cat("\n=============================================================================\n")
cat("CONCLUSION\n")
cat("=============================================================================\n\n")

if (int1_p > .05 && int2_p > .05) {
  cat("The interaction between Autonomous and Controlled motivation is\n")
  cat("NON-SIGNIFICANT in both models (p > .05).\n\n")
  cat("INTERPRETATION: The positive predictive value of autonomous motivation\n")
  cat("is ROBUST and does not vary as a function of controlled motivation.\n")
  cat("External pressure does not undermine the relationship between\n")
  cat("autonomous motivation and meaning in life.\n\n")
  cat("H3 SUPPORTED: Autonomous motivation predicts meaning regardless of\n")
  cat("the level of controlled motivation (no undermining effect).\n")
} else {
  cat("At least one interaction term is significant (p < .05).\n")
  cat("Further investigation needed.\n")
}

cat("\n=============================================================================\n")
cat("VALUES FOR THESIS\n")
cat("=============================================================================\n\n")
cat(sprintf("'The interaction term was non-significant both with depression\n"))
cat(sprintf("included (b = %.3f, p = %.3f) and without it (b = %.3f, p = %.3f).'\n",
            int2_b, int2_p, int1_b, int1_p))

cat("\n=============================================================================\n")
cat("Script 15 completed.\n")
cat("=============================================================================\n")
