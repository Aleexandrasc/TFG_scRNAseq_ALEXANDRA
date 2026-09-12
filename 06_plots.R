###############################################################################
# SCRIPT DE VISUALIZACIÓN GLOBAL - PLOTS
###############################################################################

library(ggplot2)
library(dplyr)
library(purrr)

# 1. Aseguramos que exista una carpeta unificada para los gráficos
dir_plots <- "results/plots"
if (!dir.exists(dir_plots)) {
  dir.create(dir_plots, recursive = TRUE)
  cat("¡Carpeta creada con éxito en:", dir_plots, "!\n")
}


list(loop = "1_sigma_inter",   label = "Sigma Inter (Biological noise)",  folder = "loop_sigma_inter",   suffix = "sigma_inter")

generar_y_guardar_plots <- function(loop_name, label_x, folder_name, file_suffix) {
  
  file_path <- paste0("results/", folder_name, "/all_metrics_", file_suffix, ".rds")
  
  if (!file.exists(file_path)) {
    cat("Aviso: No se encuentra el archivo en:", file_path, "\n")
    return(NULL)
  }
  
  cat("Procesando gráficos para:", loop_name, "\n")
  metrics <- readRDS(file_path)
  
  # Ordenamos el eje X
  metrics$loop_value <- factor(metrics$loop_value, levels = sort(unique(as.numeric(as.character(metrics$loop_value)))))
  
  ## Joaquin: ahora ya no hay que diferenciar panel de shift con panel de tie.
  p_bias <- ggplot(metrics,aes(x = loop_value, y = bias_win_fraction)) +
    geom_violin(trim = FALSE, alpha = 0.7, colour = "black", linewidth = 0.4) +
    geom_boxplot(
      width = 0.1,
      outlier.shape = NA,
      fill = "white",
      colour = "black",
      alpha = 0.8
    ) +
    labs(
      title = paste("Influence of", label_x, "on estimation bias"),
      x = label_x,
      y = "Model bias in Win Fraction predictions"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"),
      legend.position = "none"
    )#+ylim(-0.1,0.1)
  
  p_bias_fit <- ggplot(metrics,aes(x = loop_value, y = bias_fit)) +
    geom_violin(trim = FALSE, alpha = 0.7, colour = "black", linewidth = 0.4) +
    geom_boxplot(
      width = 0.1,
      outlier.shape = NA,
      fill = "white",
      colour = "black",
      alpha = 0.8
    ) +
    labs(
      title = paste("Influence of", label_x, "on estimation bias"),
      x = label_x,
      y = "Model bias (Max.likelihood point gradient)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"),
      legend.position = "none"
    )#+ylim(-0.1,0.1)
  
  p_stde <- ggplot(metrics,aes(x = loop_value, y = sd_error_win_fraction)) +
    geom_violin(trim = FALSE, alpha = 0.7, colour = "black", linewidth = 0.4) +
    geom_boxplot(
      width = 0.1,
      outlier.shape = NA,
      fill = "white",
      colour = "black",
      alpha = 0.8
    ) +
    labs(
      title = paste("Influence of", label_x, "on model standard error"),
      x = label_x,
      y = "Model Standard error in Win Fraction predictions"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"),
      legend.position = "none"
    )#+ylim(0,0.05)
    
    p_Fdr <- ggplot(metrics,aes(x = loop_value, y = -log10(Fdr))) +
      geom_violin(trim = FALSE, alpha = 0.7, colour = "black", linewidth = 0.4) +
      geom_boxplot(
        width = 0.1,
        outlier.shape = NA,
        fill = "white",
        colour = "black",
        alpha = 0.8
      ) +
      labs(
        title = paste("Influence of", label_x, "on statistical signal strength"),
        x = label_x,
        y = "-log10(Fdr) for condition effect"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic"),
        legend.position = "none"
      )#+ylim(0,0.05)

  
  ### Aqui falta un barplot que cuente el número de genes bajo 0.05 FDR.
  ### p_power_ok=(...)
  
  # Guardamos en results/plots/
  
  dir.create(paste0(dir_plots, "/",loop_name),recursive = TRUE)
  ggsave(paste0(dir_plots,"/",loop_name, "/Plot_Bias.png"), plot = p_bias, width = 9, height = 5, dpi = 300)
  ggsave(paste0(dir_plots,"/",loop_name, "/Plot_Bias_fit.png"), plot = p_bias_fit, width = 9, height = 5, dpi = 300)
  ggsave(paste0(dir_plots,"/",loop_name, "/Plot_stde.png"), plot = p_stde, width = 9, height = 5, dpi = 300)
  ggsave(paste0(dir_plots,"/",loop_name, "/Plot_Fdr.png"), plot = p_Fdr, width = 9, height = 5, dpi = 300)

}

# --- MAPEADO DE PARÁMETROS EXPERIMENTALES ---
config_loops <- list(
  list(loop = "1_sigma_inter",   label = "Sigma Inter (Biological noise)",  folder = "loop_sigma_inter",   suffix = "sigma_inter"),
  list(loop = "2_sigma_intra",   label = "Sigma Intra (Intra-UE cell noise)", folder = "loop_sigma_intra",   suffix = "sigma_intra"),
  list(loop = "3_sigma_epsilon", label = "Sigma Epsilon (Background noise)", folder = "loop_sigma_epsilon", suffix = "sigma_epsilon"),
  list(loop = "4_beta_shift",     label = "Effect Magnitude (Beta of Shift gene)", folder = "loop_beta_shift", suffix = "beta_shift"),
  list(loop = "5_n_ue",          label = "Number of Experimental Units (UEs)", folder = "loop_n_ue",          suffix = "n_ue"),
  list(loop = "6_n_cells_per_ue",label = "Number of Cells per UE",            folder = "loop_n_cells",       suffix = "n_cells")
)

config_loops <- list(
  list(loop = "1_sigma_inter",   label = "Sigma Inter (Biological noise)",  folder = "loop_sigma_inter",   suffix = "sigma_inter")
)

# Ejecutamos el pipeline gráfico completo
walk(config_loops, ~ generar_y_guardar_plots(.x$loop, .x$label, .x$folder, .x$suffix))

