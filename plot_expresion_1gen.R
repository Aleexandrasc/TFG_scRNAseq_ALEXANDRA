###############################################################################
# GENERACIÓN DE PANELES (EXPRESIÓN CELULAR + DENSIDAD BETA ESTIMADA)
###############################################################################

library(ggplot2)
library(dplyr)
library(purrr)
library(patchwork)

source("07_functions.R")

dir.create("results/paradigm_panels", showWarnings = FALSE, recursive = TRUE)

# FUNCIÓN PARA EL PANEL A (Violines por UE para el gen analizado)

plot_panel_A <- function(sim_data, target_gene = "G1", subtitle_txt = "") {
  expr_vals <- as.numeric(sim_data$matrix_expresion_genica[target_gene, ])
  df_a <- sim_data$cell_metadata
  df_a$expression <- expr_vals
  
  ue_info <- sim_data$ue_metadata[, c("UE", "condition")]
  df_a <- merge(df_a, ue_info, by = "UE")
  
  n_ue <- length(unique(df_a$UE))
  df_a$UE <- factor(df_a$UE, levels = paste0("U", 1:n_ue))
  df_a$condition <- factor(df_a$condition, levels = c(0, 1), labels = c("Control", "Treated"))
  
  ggplot(df_a, aes(x = UE, y = expression, fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.65, colour = "black", linewidth = 0.35) +
    geom_boxplot(width = 0.12, fill = "white", colour = "black", outlier.shape = NA, alpha = 0.85) +
    scale_fill_manual(values = c("Control" = "#4E79A7", "Treated" = "#F28E2B")) +
    labs(
      title = paste0("A. Single-cell expression profiles (", target_gene, ")"),
      subtitle = subtitle_txt,
      x = "Experimental Unit (UE)",
      y = "Gene expression level",
      fill = "Condition:"
    ) +
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 16),
      plot.subtitle = element_text(face = "italic", size = 12, hjust = 0),
      axis.title = element_text(face = "bold", size = 13),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, colour = "black"),
      axis.text.y = element_text(size = 11, colour = "black"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "bottom"
    )
}


# FUNCIÓN PARA EL PANEL B (Densidades Beta)

plot_beta_density <- function(df) {
  alpha1 <- df$mu[1] * df$phi[1]
  beta1  <- (1 - df$mu[1]) * df$phi[1]
  
  alpha2 <- df$mu[2] * df$phi[2]
  beta2  <- (1 - df$mu[2]) * df$phi[2]
  
  sub_txt <- "Green: Within-condition pairs (Ties)\nSalmon: Between-condition pairs (Shifts)"
  
  ggplot(data.frame(x = c(0, 1)), aes(x = x)) +
    stat_function(
      fun = function(x) dbeta(x, alpha1, beta1),
      xlim = c(0, 1),
      linewidth = 1,
      colour = "green4"
    ) +
    stat_function(
      fun = function(x) dbeta(x, alpha2, beta2),
      xlim = c(0, 1),
      linewidth = 1,
      colour = "salmon"
    ) +
    coord_cartesian(xlim = c(0, 1), expand = FALSE) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = "B. Estimated Beta distribution",
      subtitle = sub_txt,
      x = "Cellular win probability (p)",
      y = "Probability density f(p)"
    ) +
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 16),
      plot.subtitle = element_text(face = "italic", size = 11, hjust = 0, lineheight = 1.1),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 12, colour = "black")
    )
}


# FUNCIÓN PARA EXTRAER PARÁMETROS DEL FIT

get_mu_phi <- function(fit) {
  beta_hat  <- fixef(fit)$cond
  gamma_hat <- fixef(fit)$disp
  delta     <- fit$obj$env$data$X[, 1]
  
  eta  <- beta_hat["Delta_condition"] * delta
  zeta <- gamma_hat["(Intercept)"] + gamma_hat["Delta_condition"] * delta
  
  data.frame(
    Delta_condition = delta,
    mu  = inv_logit(eta),
    phi = exp(zeta)
  )
}

# INPUTS

set.seed(123)

n_cells_ue <- rep(50, 20)
vars <- list(condition = c(rep(0, 10), rep(1, 10)))
pesos <- list(condition = rep(1, 100))

mu_beta0 <- 2
sigma_beta0 <- 0.01

mean_sigma_inter <- 0.5
sd_sigma_inter <- 0.03

mean_sigma_intra <- 0.5
sd_sigma_intra <- 0.03

mean_sigma_epsilon <- 0.5
sd_sigma_epsilon <- 0.03


# CASO 1: beta_cond = 1.5 (Biological Effect Dominance)

pesos <- list(condition = rep(1.5, 100))

sim_data <- simulate_dataset(
  vars, pesos,
  mu_beta0, sigma_beta0,
  mean_sigma_epsilon, sd_sigma_epsilon,
  mean_sigma_inter, sd_sigma_inter,
  mean_sigma_intra, sd_sigma_intra,
  n_cells_ue
)

sim_data <- build_differential_objects(sim_data)
resultados <- run_differential_analysis(sim_data)

tab_mu_phis <- get_mu_phi(resultados$fits[[1]])
tab_mu_phis <- tab_mu_phis[c("U1-U2", "U1-U11"), ]

# Gráficas
p_left_1  <- plot_panel_A(sim_data, "G1", "Beta = 1.5, Sigma Inter = 0.5, Sigma Intra = 0.5")
p_right_1 <- plot_beta_density(tab_mu_phis)
panel_1   <- p_left_1 + p_right_1 + plot_layout(ncol = 2, widths = c(1.2, 1))

ggsave("results/paradigm_panels/Beta_domina.pdf", plot = panel_1, width = 11, height = 4.8)
ggsave("results/paradigm_panels/Beta_domina.png", plot = panel_1, width = 11, height = 4.8, dpi = 300)
cat("Guardado Caso 1\n")

pesos <- list(condition = rep(1, 100))


# CASO 2: sigma_inter = 5 (Inter-Individual Dominance)

mean_sigma_inter <- 5

sim_data <- simulate_dataset(
  vars, pesos,
  mu_beta0, sigma_beta0,
  mean_sigma_epsilon, sd_sigma_epsilon,
  mean_sigma_inter, sd_sigma_inter,
  mean_sigma_intra, sd_sigma_intra,
  n_cells_ue
)

sim_data <- build_differential_objects(sim_data)
resultados <- run_differential_analysis(sim_data)

tab_mu_phis <- get_mu_phi(resultados$fits[[1]])
tab_mu_phis <- tab_mu_phis[c("U1-U2", "U1-U11"), ]

# Gráficas
p_left_2  <- plot_panel_A(sim_data, "G1", "Beta = 1.0, Sigma Inter = 5.0, Sigma Intra = 0.5")
p_right_2 <- plot_beta_density(tab_mu_phis)
panel_2   <- p_left_2 + p_right_2 + plot_layout(ncol = 2, widths = c(1.2, 1))

ggsave("results/paradigm_panels/Sigma_inter_domina.pdf", plot = panel_2, width = 11, height = 4.8)
ggsave("results/paradigm_panels/Sigma_inter_domina.png", plot = panel_2, width = 11, height = 4.8, dpi = 300)
cat("Guardado Caso 2\n")

mean_sigma_inter <- 0.5


# CASO 3: sigma_intra = 5 (Intra-Individual Dominance)

mean_sigma_intra <- 5

sim_data <- simulate_dataset(
  vars, pesos,
  mu_beta0, sigma_beta0,
  mean_sigma_epsilon, sd_sigma_epsilon,
  mean_sigma_inter, sd_sigma_inter,
  mean_sigma_intra, sd_sigma_intra,
  n_cells_ue
)

sim_data <- build_differential_objects(sim_data)
resultados <- run_differential_analysis(sim_data)

tab_mu_phis <- get_mu_phi(resultados$fits[[1]])
tab_mu_phis <- tab_mu_phis[c("U1-U2", "U1-U11"), ]

# Gráficas
p_left_3  <- plot_panel_A(sim_data, "G1", "Beta = 1.0, Sigma Inter = 0.5, Sigma Intra = 5.0")
p_right_3 <- plot_beta_density(tab_mu_phis)
panel_3   <- p_left_3 + p_right_3 + plot_layout(ncol = 2, widths = c(1.2, 1))

ggsave("results/paradigm_panels/Sigma_intra_domina.pdf", plot = panel_3, width = 11, height = 4.8)
ggsave("results/paradigm_panels/Sigma_intra_domina.png", plot = panel_3, width = 11, height = 4.8, dpi = 300)
cat("Guardado Caso 3\n")

mean_sigma_intra <- 0.5 
