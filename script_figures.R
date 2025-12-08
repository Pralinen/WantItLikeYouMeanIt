################################################################################
# THESIS FIGURES: Figure 3 (Path Model) and Figure 4 (Forest Plot)
################################################################################

library(ggplot2)
library(brms)
library(tidyr)
library(dplyr)

# Create Figures directory if needed
if (!dir.exists("Figures")) {
  dir.create("Figures", recursive = TRUE)
}

################################################################################
# FIGURE 3: Bayesian Path Analysis Results (The Structural Model)
################################################################################

cat("Creating Figure 3: Bayesian Path Analysis (Structural Model)...\n")

# The coefficients from Script 13 (IEA Model)
# AUTO → MIL: 0.41
# CTRL → MIL: -0.04 (n.s.)
# CTRL → DEP: 0.37
# DEP → MIL: -0.81

# Using ggplot2 to draw a path diagram manually
# This gives more control than semPlot for publication-quality figures

# Define node positions
nodes <- data.frame(
  name = c("AUTO", "CTRL", "DEP", "MIL"),
  x = c(0, 0, 1.5, 3),
  y = c(2, 0, 0, 1),
  label = c("Autonomous\nMotivation", "Controlled\nMotivation", "Depression", "Meaning\nin Life")
)

# Define edges with coefficients
edges <- data.frame(
  from_x = c(0, 0, 0, 1.5),
  from_y = c(2, 0, 0, 0),
  to_x = c(3, 3, 1.5, 3),
  to_y = c(1, 1, 0, 1),
  label = c("β = 0.41***", "β = -0.04", "β = 0.37***", "β = -0.81***"),
  significant = c(TRUE, FALSE, TRUE, TRUE),
  curve = c(0.2, -0.2, 0, 0)  # curvature for paths
)

# Create the path diagram
fig3 <- ggplot() +
  # Draw edges (arrows)
  # AUTO → MIL (curved up)
  geom_curve(aes(x = 0.3, y = 2, xend = 2.7, yend = 1.2),
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             curvature = 0.2, linewidth = 1.2, color = "black") +
  # CTRL → MIL (curved down, dashed for n.s.)
  geom_curve(aes(x = 0.3, y = 0.1, xend = 2.7, yend = 0.85),
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             curvature = -0.15, linewidth = 0.8, linetype = "dashed", color = "gray50") +
  # CTRL → DEP
  geom_segment(aes(x = 0.3, y = 0, xend = 1.2, yend = 0),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               linewidth = 1.2, color = "black") +
  # DEP → MIL
  geom_segment(aes(x = 1.8, y = 0.1, xend = 2.7, yend = 0.9),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               linewidth = 1.2, color = "black") +
  # Draw nodes (rectangles)
  geom_rect(aes(xmin = -0.6, xmax = 0.6, ymin = 1.5, ymax = 2.5),
            fill = "white", color = "black", linewidth = 1) +
  geom_rect(aes(xmin = -0.6, xmax = 0.6, ymin = -0.5, ymax = 0.5),
            fill = "white", color = "black", linewidth = 1) +
  geom_rect(aes(xmin = 0.9, xmax = 2.1, ymin = -0.5, ymax = 0.5),
            fill = "gray90", color = "black", linewidth = 1) +
  geom_rect(aes(xmin = 2.4, xmax = 3.6, ymin = 0.5, ymax = 1.5),
            fill = "white", color = "black", linewidth = 1) +
  # Node labels
  geom_text(aes(x = 0, y = 2, label = "Autonomous\nMotivation"),
            size = 3.5, fontface = "bold") +
  geom_text(aes(x = 0, y = 0, label = "Controlled\nMotivation"),
            size = 3.5, fontface = "bold") +
  geom_text(aes(x = 1.5, y = 0, label = "Depression"),
            size = 3.5, fontface = "bold") +
  geom_text(aes(x = 3, y = 1, label = "Meaning\nin Life"),
            size = 3.5, fontface = "bold") +
  # Path labels
  geom_label(aes(x = 1.5, y = 1.85, label = "b = 0.41***"),
             size = 3.2, fill = "white", linewidth = 0) +
  geom_label(aes(x = 1.8, y = 0.35, label = "b = -0.04"),
             size = 3.2, fill = "white", linewidth = 0, color = "gray50") +
  geom_label(aes(x = 0.75, y = -0.25, label = "b = 0.37***"),
             size = 3.2, fill = "white", linewidth = 0) +
  geom_label(aes(x = 2.5, y = 0.35, label = "b = -0.81***"),
             size = 3.2, fill = "white", linewidth = 0) +
  # Add note about indirect effect
  annotate("text", x = 1.5, y = -0.9,
           label = "Indirect effect (CTRL -> DEP -> MIL): b = -0.30***",
           size = 3, fontface = "italic") +
  # Styling
  coord_fixed(ratio = 0.8, xlim = c(-1, 4), ylim = c(-1.2, 2.8)) +
  theme_void() +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9, face = "italic")
  ) +
  labs(
    title = "IEA Framework: Dual Pathway Model",
    caption = "Note: ***p < .001. Solid lines = credible effects; Dashed line = non-credible (95% CI includes 0).\nBayesian SEM with N = 437. Depression shown in gray as mediator."
  )

# Save Figure 3
ggsave("Figures/Figure3_Path_Model.png", fig3, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("Figures/Figure3_Path_Model.pdf", fig3, width = 8, height = 6, bg = "white")

cat("  Figure 3 saved to Figures/Figure3_Path_Model.png and .pdf\n\n")

################################################################################
# FIGURE 4: Bayesian Posterior Distributions (Forest Plot)
################################################################################

cat("Creating Figure 4: Forest Plot of Fixed Effects...\n")

# Load Script 10 results
results10 <- readRDS("Results/Script10_Combined.rds")
brm_fit <- results10$model3$fit

# Extract posterior samples for key predictors
post <- as_draws_df(brm_fit)

# Get the fixed effects we care about
key_effects <- data.frame(
  Parameter = c("Autonomous", "Controlled", "Depression"),
  Estimate = c(
    mean(post$b_Autonomous),
    mean(post$b_Controlled),
    mean(post$b_Depression)
  ),
  Q2.5 = c(
    quantile(post$b_Autonomous, 0.025),
    quantile(post$b_Controlled, 0.025),
    quantile(post$b_Depression, 0.025)
  ),
  Q97.5 = c(
    quantile(post$b_Autonomous, 0.975),
    quantile(post$b_Controlled, 0.975),
    quantile(post$b_Depression, 0.975)
  )
)

# Reorder for display (Depression at bottom, then Controlled, then Autonomous at top)
key_effects$Parameter <- factor(key_effects$Parameter,
                                 levels = c("Depression", "Controlled", "Autonomous"))

# Create forest plot
fig4 <- ggplot(key_effects, aes(x = Estimate, y = Parameter)) +
  # Add vertical line at 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  # Add credible intervals
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5),
                 height = 0.2, linewidth = 0.8, color = "black") +
  # Add point estimates
  geom_point(size = 4, shape = 21, fill = "black", color = "black") +
  # Add estimate labels
  geom_text(aes(label = sprintf("%.3f", Estimate)),
            hjust = -0.3, vjust = -0.8, size = 3.5) +
  # Styling
  scale_x_continuous(
    name = "Standardized Effect Size (Probit Scale)",
    limits = c(-0.25, 0.12),
    breaks = seq(-0.25, 0.1, 0.05)
  ) +
  scale_y_discrete(name = "Predictors") +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9, face = "italic")
  ) +
  labs(
    title = "Fixed Effects on Meaning in Life (Item-Level)",
    caption = "Note: Points show posterior means; error bars show 95% Bayesian credible intervals.\nEffects crossing zero are not credibly different from zero."
  )

# Save Figure 4
ggsave("Figures/Figure4_Forest_Plot.png", fig4, width = 8, height = 5, dpi = 300, bg = "white")
ggsave("Figures/Figure4_Forest_Plot.pdf", fig4, width = 8, height = 5, bg = "white")

cat("  Figure 4 saved to Figures/Figure4_Forest_Plot.png and .pdf\n\n")

################################################################################
# BONUS: Combined version showing AUTO vs CTRL contrast more clearly
################################################################################

cat("Creating Bonus Figure: AUTO vs CTRL Contrast...\n")

# Just AUTO and CTRL for cleaner comparison
contrast_effects <- key_effects[key_effects$Parameter %in% c("Autonomous", "Controlled"), ]
contrast_effects$Parameter <- factor(contrast_effects$Parameter,
                                      levels = c("Controlled", "Autonomous"))

fig4b <- ggplot(contrast_effects, aes(x = Estimate, y = Parameter)) +
  # Shaded region for "non-credible" zone around 0
  annotate("rect", xmin = -0.02, xmax = 0.02, ymin = 0, ymax = 3,
           fill = "gray90", alpha = 0.5) +
  # Vertical line at 0
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  # Credible intervals
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5),
                 height = 0.15, linewidth = 1, color = "black") +
  # Point estimates
  geom_point(size = 5, shape = 21,
             aes(fill = ifelse(Q2.5 > 0 | Q97.5 < 0, "Credible", "Not Credible")),
             color = "black") +
  scale_fill_manual(values = c("Credible" = "black", "Not Credible" = "white"),
                    guide = "none") +
  # Labels
  geom_text(aes(label = sprintf("b = %.3f [%.3f, %.3f]", Estimate, Q2.5, Q97.5)),
            hjust = -0.1, vjust = -1.2, size = 3.5) +
  # Styling
  scale_x_continuous(
    name = "Effect on Meaning in Life (Probit Scale)",
    limits = c(-0.03, 0.11),
    breaks = seq(-0.02, 0.10, 0.02)
  ) +
  scale_y_discrete(name = "") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9, face = "italic")
  ) +
  labs(
    title = "Autonomous vs Controlled Motivation Effects",
    caption = "Note: Filled points indicate credible effects (95% CI excludes 0).\nAutonomous motivation credibly predicts meaning; Controlled does not."
  )

ggsave("Figures/Figure4b_AUTO_vs_CTRL.png", fig4b, width = 8, height = 4, dpi = 300, bg = "white")
ggsave("Figures/Figure4b_AUTO_vs_CTRL.pdf", fig4b, width = 8, height = 4, bg = "white")

cat("  Bonus Figure saved to Figures/Figure4b_AUTO_vs_CTRL.png and .pdf\n\n")

################################################################################
# Summary
################################################################################

cat("================================================================================\n")
cat("ALL FIGURES CREATED SUCCESSFULLY\n")
cat("================================================================================\n\n")

cat("Files saved in Figures/:\n")
cat("  - Figure3_Path_Model.png / .pdf    (Structural Model Diagram)\n")
cat("  - Figure4_Forest_Plot.png / .pdf   (Forest Plot with Depression)\n")
cat("  - Figure4b_AUTO_vs_CTRL.png / .pdf (Clean AUTO vs CTRL comparison)\n")
