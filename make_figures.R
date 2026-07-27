################################################################################
# Generate publication-ready figures (600 DPI, Times New Roman)
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(showtext)
  library(sysfonts)
  library(grid)
  library(gridExtra)
})

# Setup Times New Roman
font_add("Times New Roman",
         regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
         bold = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",
         italic = "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf",
         bolditalic = "/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf")
showtext_auto()
showtext_opts(dpi = 600)

if (!dir.exists("Figures")) dir.create("Figures", recursive = TRUE)

# Text sizes (showtext uses mm-like units; these are calibrated for 600 DPI)
node_size  <- 3.5
coef_size  <- 3.0
label_size <- 3.8

################################################################################
# FIGURE 1: Three Competing DAGs
################################################################################
cat("Creating Figure 1: Competing DAGs...\n")

make_dag <- function(panel_label, pred_labels, pred_fills, pred_borders,
                     pred_lty, arr_lwd, arr_lty, arr_col,
                     coef_labels, coef_colors) {

  p <- ggplot() +
    coord_fixed(ratio = 1, xlim = c(-1.3, 3.3), ylim = c(-0.8, 3.5)) +
    theme_void() +
    theme(text = element_text(family = "Times New Roman"),
          plot.margin = margin(2, 2, 2, 2))

  # Panel label
  p <- p + annotate("text", x = 1, y = 3.3, label = panel_label,
                     size = label_size, family = "Times New Roman", fontface = "bold")

  if (length(pred_labels) == 1) {
    p <- p +
      geom_rect(aes(xmin = 0, xmax = 2, ymin = 1.6, ymax = 2.4),
                fill = pred_fills[1], color = pred_borders[1],
                linewidth = 0.7, linetype = pred_lty[1]) +
      annotate("text", x = 1, y = 2, label = pred_labels[1],
               size = node_size, family = "Times New Roman", fontface = "bold") +
      geom_rect(aes(xmin = 0, xmax = 2, ymin = -0.4, ymax = 0.4),
                fill = "#D6E4F0", color = "black", linewidth = 0.7) +
      annotate("text", x = 1, y = 0, label = "Meaning",
               size = node_size, family = "Times New Roman", fontface = "bold") +
      geom_segment(aes(x = 1, y = 1.6, xend = 1, yend = 0.4),
                   arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                   linewidth = arr_lwd[1], linetype = arr_lty[1], color = arr_col[1]) +
      annotate("text", x = 1.4, y = 1, label = coef_labels[1],
               size = coef_size, family = "Times New Roman", fontface = "italic",
               color = coef_colors[1], hjust = 0)
  } else {
    p <- p +
      geom_rect(aes(xmin = -0.7, xmax = 0.7, ymin = 1.6, ymax = 2.4),
                fill = pred_fills[1], color = pred_borders[1],
                linewidth = 0.7, linetype = pred_lty[1]) +
      annotate("text", x = 0, y = 2, label = pred_labels[1],
               size = node_size, family = "Times New Roman", fontface = "bold",
               color = ifelse(pred_lty[1] == "dashed", "gray50", "black")) +
      geom_rect(aes(xmin = 1.3, xmax = 2.7, ymin = 1.6, ymax = 2.4),
                fill = pred_fills[2], color = pred_borders[2],
                linewidth = 0.7, linetype = pred_lty[2]) +
      annotate("text", x = 2, y = 2, label = pred_labels[2],
               size = node_size, family = "Times New Roman", fontface = "bold",
               color = ifelse(pred_lty[2] == "dashed", "gray50", "black")) +
      geom_rect(aes(xmin = 0.3, xmax = 1.7, ymin = -0.4, ymax = 0.4),
                fill = "#D6E4F0", color = "black", linewidth = 0.7) +
      annotate("text", x = 1, y = 0, label = "Meaning",
               size = node_size, family = "Times New Roman", fontface = "bold") +
      geom_segment(aes(x = 0, y = 1.6, xend = 0.8, yend = 0.4),
                   arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                   linewidth = arr_lwd[1], linetype = arr_lty[1], color = arr_col[1]) +
      geom_segment(aes(x = 2, y = 1.6, xend = 1.2, yend = 0.4),
                   arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
                   linewidth = arr_lwd[2], linetype = arr_lty[2], color = arr_col[2]) +
      annotate("text", x = -0.25, y = 0.95, label = coef_labels[1],
               size = coef_size, family = "Times New Roman", fontface = "italic",
               color = coef_colors[1], hjust = 1) +
      annotate("text", x = 2.25, y = 0.95, label = coef_labels[2],
               size = coef_size, family = "Times New Roman", fontface = "italic",
               color = coef_colors[2], hjust = 0)
  }
  return(p)
}

dag1 <- make_dag("A. GSC Only", "GSC", "#E8D5F0", "black", "solid",
                 0.6, "solid", "black", "b", "black")

dag2 <- make_dag("B. Autonomous Only\n(Hypothesized)",
                 c("Auto", "Ctrl"), c("#D5E8D4", "#F8D7DA"), c("black", "gray60"),
                 c("solid", "dashed"), c(0.6, 0.4), c("solid", "solid"), c("black", "gray60"),
                 c("b***", "b = 0"), c("black", "gray50"))

dag3 <- make_dag("C. Both Predictors",
                 c("Auto", "Ctrl"), c("#D5E8D4", "#F8D7DA"), c("black", "black"),
                 c("solid", "solid"), c(0.6, 0.4), c("solid", "dashed"), c("black", "gray50"),
                 c("b***", "b NS"), c("black", "gray50"))

fig1 <- gridExtra::arrangeGrob(dag1, dag2, dag3, ncol = 3)

ggsave("Figures/Figure1_DAGs.png", fig1, width = 7.5, height = 3, dpi = 600, bg = "white")
ggsave("Figures/Figure1_DAGs.pdf", fig1, width = 7.5, height = 3, bg = "white", device = cairo_pdf)
cat("  Figure 1 saved.\n")

################################################################################
# FIGURE 2: Forest Plot (hardcoded from verified Script 10 Model 3 output)
################################################################################
cat("Creating Figure 2: Forest Plot...\n")

# Values from Script 10 Model 3 Combined (verified earlier)
key_effects <- data.frame(
  Parameter = factor(c("Autonomous", "Controlled", "Depression"),
                     levels = c("Depression", "Controlled", "Autonomous")),
  Estimate = c(0.0699, -0.0056, -0.1831),
  Q2.5     = c(0.0508, -0.0190, -0.2304),
  Q97.5    = c(0.0888,  0.0077, -0.1384)
)

fig2 <- ggplot(key_effects, aes(x = Estimate, y = Parameter)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5),
                 height = 0.2, linewidth = 0.6, color = "black") +
  geom_point(size = 3, shape = 16, color = "black") +
  geom_text(aes(label = sprintf("%.3f", Estimate)),
            hjust = -0.3, vjust = -0.8, size = coef_size, family = "Times New Roman") +
  scale_x_continuous(
    name = "Standardized Effect Size (Probit Scale)",
    limits = c(-0.26, 0.12),
    breaks = seq(-0.25, 0.10, 0.05)
  ) +
  scale_y_discrete(name = "Predictors") +
  theme_classic(base_size = 10) +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 9),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.2),
    plot.margin = margin(10, 15, 10, 10)
  )

ggsave("Figures/Figure2_Forest_Plot.png", fig2, width = 5.5, height = 3, dpi = 600, bg = "white")
ggsave("Figures/Figure2_Forest_Plot.pdf", fig2, width = 5.5, height = 3, bg = "white", device = cairo_pdf)
cat("  Figure 2 saved.\n")

################################################################################
# FIGURE 3: IEA Path Model (from Script 13 RDS)
################################################################################
cat("Creating Figure 3: IEA Path Model...\n")

r13 <- readRDS("Results/Latent/Script13_IEA_Results.rds")
pp <- r13$paths

ax <- 0;   ay <- 2.5   # Autonomous
cx <- 0;   cy <- 0     # Controlled
dx <- 2.5; dy <- 0     # Depression
mx <- 4;   my <- 1.5   # Meaning
bw <- 0.9; bh <- 0.45

fig3 <- ggplot() +
  coord_fixed(ratio = 0.75, xlim = c(-1.5, 5.2), ylim = c(-1.3, 3.2)) +
  theme_void() +
  theme(text = element_text(family = "Times New Roman"),
        plot.margin = margin(5, 10, 10, 10)) +

  # Nodes
  geom_rect(aes(xmin = ax-bw, xmax = ax+bw, ymin = ay-bh, ymax = ay+bh),
            fill = "white", color = "black", linewidth = 0.8) +
  annotate("text", x = ax, y = ay, label = "Autonomous\nMotivation",
           size = node_size, family = "Times New Roman", fontface = "bold", lineheight = 0.85) +

  geom_rect(aes(xmin = cx-bw, xmax = cx+bw, ymin = cy-bh, ymax = cy+bh),
            fill = "white", color = "black", linewidth = 0.8) +
  annotate("text", x = cx, y = cy, label = "Controlled\nMotivation",
           size = node_size, family = "Times New Roman", fontface = "bold", lineheight = 0.85) +

  geom_rect(aes(xmin = dx-bw, xmax = dx+bw, ymin = dy-bh, ymax = dy+bh),
            fill = "gray90", color = "black", linewidth = 0.8) +
  annotate("text", x = dx, y = dy, label = "Depression",
           size = node_size, family = "Times New Roman", fontface = "bold") +

  geom_rect(aes(xmin = mx-bw, xmax = mx+bw, ymin = my-bh, ymax = my+bh),
            fill = "white", color = "black", linewidth = 0.8) +
  annotate("text", x = mx, y = my, label = "Meaning\nin Life",
           size = node_size, family = "Times New Roman", fontface = "bold", lineheight = 0.85) +

  # AUTO -> MIL
  geom_curve(aes(x = ax + bw, y = ay, xend = mx - bw, yend = my + 0.15),
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
             curvature = 0.15, linewidth = 0.8, color = "black") +
  annotate("text", x = 2.2, y = 2.45,
           label = sprintf("b = %.2f***", round(pp$auto_direct$est, 2)),
           size = coef_size, family = "Times New Roman") +

  # CTRL -> MIL (dashed, n.s.)
  geom_curve(aes(x = cx + bw, y = cy + 0.3, xend = mx - bw, yend = my - 0.2),
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             curvature = -0.15, linewidth = 0.4, linetype = "dashed", color = "gray50") +
  annotate("text", x = 2.3, y = 1.0,
           label = sprintf("b = %.2f", round(pp$ctrl_direct$est, 2)),
           size = coef_size, family = "Times New Roman", color = "gray50") +

  # CTRL -> DEP
  geom_segment(aes(x = cx + bw, y = cy, xend = dx - bw, yend = dy),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth = 0.8, color = "black") +
  annotate("text", x = 1.25, y = -0.35,
           label = sprintf("b = %.2f***", round(pp$ctrl_to_dep$est, 2)),
           size = coef_size, family = "Times New Roman") +

  # DEP -> MIL
  geom_segment(aes(x = dx + bw, y = dy + 0.15, xend = mx - bw, yend = my - 0.3),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               linewidth = 0.8, color = "black") +
  annotate("text", x = 3.6, y = 0.5,
           label = sprintf("b = %.2f***", round(pp$dep_to_mil$est, 2)),
           size = coef_size, family = "Times New Roman") +

  # Indirect effect
  annotate("text", x = 2, y = -1.1,
           label = sprintf("Indirect effect (CTRL \u2192 DEP \u2192 MIL): b = %.2f***",
                          round(pp$distress_path$est, 2)),
           size = coef_size, family = "Times New Roman", fontface = "italic")

ggsave("Figures/Figure3_IEA_Model.png", fig3, width = 7, height = 4.5, dpi = 600, bg = "white")
ggsave("Figures/Figure3_IEA_Model.pdf", fig3, width = 7, height = 4.5, bg = "white", device = cairo_pdf)
cat("  Figure 3 saved.\n")

cat("\nAll figures saved to Figures/ at 600 DPI with Times New Roman.\n")
