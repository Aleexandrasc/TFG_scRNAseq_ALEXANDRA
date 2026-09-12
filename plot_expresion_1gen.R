###############################################################################
#PLOTS DE EXPRESIÓN DE UN GEN CONCRETO 
###############################################################################

library(ggplot2)
library(patchwork)
library(glmmTMB)

source("06_functions.R")

dir.create("results/paradigm_panels", showWarnings = FALSE, recursive = TRUE)

build_panel_from_rds <- function(path_sim_data, 
                                 path_resultados, 
                                 output_filename, 
                                 regime_title, 
                                 parameter_subtitle, 
                                 target_gene = "G1") {
  
  if (!file.exists(path_sim_data) || !file.exists(path_resultados)) {
    cat("Warning: Files not found at:\n", path_sim_data, "\n")
    return(NULL)
  }
  
  cat("Generating:", regime_title, "\n")
  
  sim_data   <- readRDS(path_sim_data)
  resultados <- readRDS(path_resultados)
  
  # ---------------------------------------------------------------------------
  # Panel A: Single-cell expression per experimental unit
  # ---------------------------------------------------------------------------
  expr_vals <- as.numeric(sim_data$matrix_expresion_genica[target_gene, ])
  df_a <- sim_data$cell_metadata
  df_a$expression <- expr_vals
  
  ue_info <- sim_data$ue_metadata[, c("UE", "condition")]
  df_a <- merge(df_a, ue_info, by = "UE")
  
  n_ue <- length(unique(df_a$UE))
  df_a$UE <- factor(df_a$UE, levels = paste0("U", 1:n_ue))
  df_a$condition <- factor(df_a$condition, levels = c(0, 1), labels = c("Control", "Treated"))
  
  p_expr <- ggplot(df_a, aes(x = UE, y = expression, fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.65, colour = "black", linewidth = 0.35) +
    geom_boxplot(width = 0.12, fill = "white", colour = "black", outlier.shape = NA, alpha = 0.85) +
    scale_fill_manual(values = c("Control" = "#4E79A7", "Treated" = "#F28E2B")) +
    labs(
      title = paste0("A. Single-cell expression profiles (", target_gene, ")"),
      subtitle = paste0(regime_title, " (", parameter_subtitle, ")"),
      x = "Experimental Unit (UE)",
      y = "Gene expression level",
      fill = "Experimental Condition:"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(face = "italic", size = 9.5, hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom"
    )
  
  # ---------------------------------------------------------------------------
  # Panel B: Theoretical Beta distributions of winning probabilities
  # ---------------------------------------------------------------------------
  fit_gen <- resultados$fits[[target_gene]]
  
  X_data <- fit_gen$obj$env$data$X
  p_dim  <- ncol(X_data)
  
  beta_hat  <- fit_gen$fit$par[1:p_dim]
  gamma_hat <- fit_gen$fit$par[(p_dim + 1):(2 * p_dim)]
  
  # Within-condition comparison (Baseline: Delta = 0)
  mu_base   <- 0.5
  phi_base  <- 20
  alpha_base <- mu_base * phi_base
  beta_base  <- (1 - mu_base) * phi_base
  
  # Between-condition comparison (Treated vs Control: Delta = 1)
  eta_diff   <- as.numeric(1 * beta_hat)
  mu_diff    <- inv_logit(eta_diff)
  
  zeta_diff  <- as.numeric(1 * gamma_hat)
  phi_diff   <- exp(zeta_diff)
  
  alpha_diff <- max(mu_diff * phi_diff, 0.0001)
  beta_diff  <- max((1 - mu_diff) * phi_diff, 0.0001)
  
  p_grid <- seq(0.001, 0.999, length.out = 1000)
  
  # Etiquetas en dos líneas con centrado
  label_within  <- "Within-condition pairs\n(Control vs Control)"
  label_between <- "Between-condition pairs\n(Treated vs Control)"
  
  df_densities <- data.frame(
    p = rep(p_grid, 2),
    density = c(dbeta(p_grid, alpha_base, beta_base),
                dbeta(p_grid, alpha_diff, beta_diff)),
    pair_type = factor(
      rep(c(label_within, label_between), each = length(p_grid)),
      levels = c(label_within, label_between)
    )
  )
  
  p_beta <- ggplot(df_densities, aes(x = p, y = density, colour = pair_type, fill = pair_type)) +
    geom_line(linewidth = 0.8) +
    geom_area(alpha = 0.25, position = "identity") +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "gray30", linewidth = 0.5) +
    scale_colour_manual(values = setNames(c("#59A14F", "#E15759"), c(label_within, label_between))) +
    scale_fill_manual(values = setNames(c("#59A14F", "#E15759"), c(label_within, label_between))) +
    labs(
      title = "B. Estimated Beta distribution of winning probabilities",
      subtitle = paste0("Baseline: alpha = ", format(round(alpha_base, 3), nsmall = 3), 
                        ", beta = ", format(round(beta_base, 3), nsmall = 3), 
                        " | Differential: alpha = ", format(round(alpha_diff, 3), nsmall = 3), 
                        ", beta = ", format(round(beta_diff, 3), nsmall = 3)),
      x = "Cellular win probability (p)",
      y = "Probability density f(p)",
      colour = "Comparison type:",
      fill = "Comparison type:"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(face = "italic", size = 9.5, hjust = 0),
      legend.position = "bottom",
      legend.text = element_text(hjust = 0.5, size = 8.5)
    )
  
  # ---------------------------------------------------------------------------
  # Combination and export
  # ---------------------------------------------------------------------------
  combined_panel <- p_expr + p_beta + plot_layout(ncol = 2, widths = c(1.15, 1))
  
  out_path <- paste0("results/paradigm_panels/", output_filename, ".png")
  ggsave(out_path, plot = combined_panel, width = 12, height = 5.2, dpi = 300)
  cat("Panel saved successfully to:", out_path, "\n\n")
  
  return(combined_panel)
}

# ==============================================================================
# EXECUTION ACROSS THE 3 REPRESENTATIVE REGIMES
# ==============================================================================

# 1. Biological effect dominant regime
build_panel_from_rds(
  path_sim_data   = "results/loop_beta_shift/beta_shift_1_5_sim_data.rds",
  path_resultados = "results/loop_beta_shift/beta_shift_1_5_resultados.rds",
  output_filename = "Panel_1_Beta_Dominant",
  regime_title    = "Biological Effect Dominance",
  parameter_subtitle = "Beta = 1.5, Sigma Inter = 0.5, Sigma Intra = 0.5"
)

# 2. Inter-individual variability dominant regime
build_panel_from_rds(
  path_sim_data   = "results/loop_sigma_inter/sigma_inter_5_sim_data.rds",
  path_resultados = "results/loop_sigma_inter/sigma_inter_5_resultados.rds",
  output_filename = "Panel_2_Sigma_Inter_Dominant",
  regime_title    = "Inter-Individual Heterogeneity Dominance",
  parameter_subtitle = "Beta = 1.0, Sigma Inter = 5.0, Sigma Intra = 0.5"
)

# 3. Intra-individual variability dominant regime
build_panel_from_rds(
  path_sim_data   = "results/loop_sigma_intra/sigma_intra_5_sim_data.rds",
  path_resultados = "results/loop_sigma_intra/sigma_intra_5_resultados.rds",
  output_filename = "Panel_3_Sigma_Intra_Dominant",
  regime_title    = "Intra-Individual Cellular Noise Dominance",
  parameter_subtitle = "Beta = 1.0, Sigma Inter = 0.5, Sigma Intra = 5.0"
)