###############################################################################
# El 'main' 
###############################################################################

{
    
    library(ggplot2)
    library(dplyr)
    library(purrr)
  
#IMPORTANTE: CADA VEZ QUE CAMBIE LA VERSIÓN CAMBIAR EL NÚMERO!!!
source("06_functions.R")

# RUTINA DE SIMULACIÓN
# Definimos los Inputs Primigenios

set.seed(123)

#1) Numero de células por unidad experimental
n_cells_ue <- rep(50, 20) 

#2) Diseño
vars <- list(
  condition = c(rep(0,10), rep(1,10)) #Nº de UE
)

# tiene que haber el mismo numero de elementos en los dos arrays anteriores

#3) Lo que quiero que pese el diseño
#de aqui podemos saber cuantos genes tiene cada célula

pesos <- list(
  condition = rep(1, 100)
)

#4) Definimos como quiero que sea mi ruido
#término basal
mu_beta0 = 2
sigma_beta0 = 0.01

# término inter
mean_sigma_inter = 0.5
sd_sigma_inter = 0.03

# término intra
mean_sigma_intra = 0.5
sd_sigma_intra = 0.03

# término ruido de fondo
mean_sigma_epsilon = 0.5
sd_sigma_epsilon = 0.03
}

###############################################################################
# LOOP 1: BARRIDO DE SIGMA_INTER
###############################################################################

{

dir.create("results", showWarnings = FALSE)
dir.create("results/loop_sigma_inter", showWarnings = FALSE)

valores_sigma_inter <- c(0.6, 0.7, 0.8, 0.9, 1, 5, 10)

all_metrics_sigma_inter <- list()

for (valor in valores_sigma_inter) {
  
  cat("Ejecutando sigma_inter =", valor, "\n")
  
  sim_data <- simulate_dataset(
    vars, pesos,
    mu_beta0, sigma_beta0,
    mean_sigma_epsilon, sd_sigma_epsilon,
    valor, sd_sigma_inter,
    mean_sigma_intra, sd_sigma_intra,
    n_cells_ue
  )
  
  sim_data <- build_differential_objects(sim_data)
  
  resultados <- run_differential_analysis(sim_data)
  
  metrics <- extract_performance_metrics(
    sim_data = sim_data,
    resultados = resultados,
    effect_name = "Delta_condition"
  )
  
  metrics$loop <- "sigma_inter"
  metrics$loop_value <- valor
      metrics$Fdr <- p.adjust(metrics$pvalue,method="BH")
  
  nombre_base <- paste0("sigma_inter_", valor)
  nombre_base <- gsub("\\.", "_", nombre_base)
  
  saveRDS(sim_data, paste0("results/loop_sigma_inter/", nombre_base, "_sim_data.rds"))
  saveRDS(resultados, paste0("results/loop_sigma_inter/", nombre_base, "_resultados.rds"))
  saveRDS(metrics, paste0("results/loop_sigma_inter/", nombre_base, "_metrics.rds"))
  write.csv(metrics, paste0("results/loop_sigma_inter/", nombre_base, "_metrics.csv"), row.names = FALSE)
  
  all_metrics_sigma_inter[[nombre_base]] <- metrics
}

all_metrics_sigma_inter <- do.call(rbind, all_metrics_sigma_inter)

saveRDS(all_metrics_sigma_inter, "results/loop_sigma_inter/all_metrics_sigma_inter.rds")
write.csv(all_metrics_sigma_inter, "results/loop_sigma_inter/all_metrics_sigma_inter.csv", row.names = FALSE)

head(all_metrics_sigma_inter)

}

###############################################################################
# LOOP 2: BARRIDO DE SIGMA_INTRA
###############################################################################

{
    dir.create("results", showWarnings = FALSE)
    dir.create("results/loop_sigma_intra", showWarnings = FALSE)

    valores_sigma_intra <- c(0.6, 0.7, 0.8, 0.9, 1, 5, 10)

    all_metrics_sigma_intra <- list()

    for (valor in valores_sigma_intra) {
      
      cat("Ejecutando sigma_intra =", valor, "\n")
      
      sim_data <- simulate_dataset(
        vars, pesos,
        mu_beta0, sigma_beta0,
        mean_sigma_epsilon, sd_sigma_epsilon,
        mean_sigma_inter, sd_sigma_inter,
        valor, sd_sigma_intra,
        n_cells_ue
      )
      
      sim_data <- build_differential_objects(sim_data)
      
      resultados <- run_differential_analysis(sim_data)
      
      metrics <- extract_performance_metrics(
        sim_data = sim_data,
        resultados = resultados,
        effect_name = "Delta_condition"
      )
      
      metrics$loop <- "sigma_intra"
      metrics$loop_value <- valor
      metrics$Fdr <- p.adjust(metrics$pvalue,method="BH")
      
      nombre_base <- paste0("sigma_intra_", valor)
      nombre_base <- gsub("\\.", "_", nombre_base)
      
      saveRDS(sim_data, paste0("results/loop_sigma_intra/", nombre_base, "_sim_data.rds"))
      saveRDS(resultados, paste0("results/loop_sigma_intra/", nombre_base, "_resultados.rds"))
      saveRDS(metrics, paste0("results/loop_sigma_intra/", nombre_base, "_metrics.rds"))
      write.csv(metrics, paste0("results/loop_sigma_intra/", nombre_base, "_metrics.csv"), row.names = FALSE)
      
      all_metrics_sigma_intra[[nombre_base]] <- metrics
    }

    all_metrics_sigma_intra <- do.call(rbind, all_metrics_sigma_intra)

    saveRDS(all_metrics_sigma_intra, "results/loop_sigma_intra/all_metrics_sigma_intra.rds")
    write.csv(all_metrics_sigma_intra, "results/loop_sigma_intra/all_metrics_sigma_intra.csv", row.names = FALSE)

    head(all_metrics_sigma_intra)
}

###############################################################################
# LOOP 3: BARRIDO DE SIGMA_EPSILON
###############################################################################

{
    dir.create("results", showWarnings = FALSE)
    dir.create("results/loop_sigma_epsilon", showWarnings = FALSE)

    valores_sigma_epsilon <- c(0.6, 0.7, 0.8, 0.9, 1, 5, 10)

    all_metrics_sigma_epsilon <- list()

    for (valor in valores_sigma_epsilon) {
      
      cat("Ejecutando sigma_epsilon =", valor, "\n")
      
      sim_data <- simulate_dataset(
        vars, pesos,
        mu_beta0, sigma_beta0,
        valor, sd_sigma_epsilon,
        mean_sigma_inter, sd_sigma_inter,
        mean_sigma_intra, sd_sigma_intra,
        n_cells_ue
      )
      
      sim_data <- build_differential_objects(sim_data)
      
      resultados <- run_differential_analysis(sim_data)
      
      metrics <- extract_performance_metrics(
        sim_data = sim_data,
        resultados = resultados,
        effect_name = "Delta_condition"
      )
      
      metrics$loop <- "sigma_epsilon"
      metrics$loop_value <- valor
      metrics$Fdr <- p.adjust(metrics$pvalue,method="BH") ##

      nombre_base <- paste0("sigma_epsilon_", valor)
      nombre_base <- gsub("\\.", "_", nombre_base)
      
      saveRDS(sim_data, paste0("results/loop_sigma_epsilon/", nombre_base, "_sim_data.rds"))
      saveRDS(resultados, paste0("results/loop_sigma_epsilon/", nombre_base, "_resultados.rds"))
      saveRDS(metrics, paste0("results/loop_sigma_epsilon/", nombre_base, "_metrics.rds"))
      write.csv(metrics, paste0("results/loop_sigma_epsilon/", nombre_base, "_metrics.csv"), row.names = FALSE)
      
      all_metrics_sigma_epsilon[[nombre_base]] <- metrics
    }

    all_metrics_sigma_epsilon <- do.call(rbind, all_metrics_sigma_epsilon)

    saveRDS(all_metrics_sigma_epsilon, "results/loop_sigma_epsilon/all_metrics_sigma_epsilon.rds")
    write.csv(all_metrics_sigma_epsilon, "results/loop_sigma_epsilon/all_metrics_sigma_epsilon.csv", row.names = FALSE)

    head(all_metrics_sigma_epsilon)
}

###############################################################################
# LOOP 4: BARRIDO DE LA MAGNITUD DEL EFECTO (BETA)
###############################################################################

{
    dir.create("results", showWarnings = FALSE)
    dir.create("results/loop_beta_shift", showWarnings = FALSE)

    # Los valores de Beta (fuerza del efecto)
    valores_beta <- c(-2,-1.5,-1,-0.5,-0.1,0,0.1,0.5, 1, 1.5,2)

    all_metrics_beta_shift <- list()

    for (valor in valores_beta) {
      
      cat("Ejecutando Beta del Shift =", valor, "\n")
      
      # Creamos los pesos dinámicos donde G1 toma el 'valor' de la Beta
      pesos_dinamicos <- list(
        condition = rep(valor, 100) 
      )
      
      # Simulación: pasamos "pesos_dinamicos" en lugar de los fijos de arriba,
      # y devolvemos sigma_beta0 a su valor de referencia (0.01)
      sim_data <- simulate_dataset(
        vars, pesos_dinamicos,
        mu_beta0, sigma_beta0, # sigma_beta0 se queda fijo en 0.01
        mean_sigma_epsilon, sd_sigma_epsilon,
        mean_sigma_inter, sd_sigma_inter,
        mean_sigma_intra, sd_sigma_intra,
        n_cells_ue
      )
      
      sim_data <- build_differential_objects(sim_data)
      
      resultados <- run_differential_analysis(sim_data)
      
      metrics <- extract_performance_metrics(
        sim_data = sim_data,
        resultados = resultados,
        effect_name = "Delta_condition"
      )
      
      metrics$loop <- "beta_shift"
      metrics$loop_value <- valor
      metrics$Fdr <- p.adjust(metrics$pvalue,method="BH")
      
      nombre_base <- paste0("beta_shift_", valor)
      nombre_base <- gsub("\\.", "_", nombre_base)
      
      saveRDS(sim_data, paste0("results/loop_beta_shift/", nombre_base, "_sim_data.rds"))
      saveRDS(resultados, paste0("results/loop_beta_shift/", nombre_base, "_resultados.rds"))
      saveRDS(metrics, paste0("results/loop_beta_shift/", nombre_base, "_metrics.rds"))
      write.csv(metrics, paste0("results/loop_beta_shift/", nombre_base, "_metrics.csv"), row.names = FALSE)
      
      all_metrics_beta_shift[[nombre_base]] <- metrics
    }

    all_metrics_beta_shift <- do.call(rbind, all_metrics_beta_shift)

    saveRDS(all_metrics_beta_shift, "results/loop_beta_shift/all_metrics_beta_shift.rds")
    write.csv(all_metrics_beta_shift, "results/loop_beta_shift/all_metrics_beta_shift.csv", row.names = FALSE)

    head(all_metrics_beta_shift)
}

###############################################################################
# LOOP 5: BARRIDO DEL NÚMERO DE UNIDADES EXPERIMENTALES (UEs)
###############################################################################

{
dir.create("results", showWarnings = FALSE)
dir.create("results/loop_n_ue", showWarnings = FALSE)

#valores_n_ue <- c(10, 20, 50, 100, 200, 500, 1000)
valores_n_ue <- c(10, 20, 40, 60, 80, 100)

all_metrics_n_ue <- list()

for (valor in valores_n_ue){
  
  cat("Ejecutando n_ue =", valor, "\n")
  
  # 1. Generamos dinámicamente las células por UE (fijadas a 50 como en la referencia)
  n_cells_ue_dinamico <- rep(50, valor)
  
  # 2. Generamos dinámicamente el diseño balanceado (mitad 0, mitad 1)
  vars_dinamico <- list(
    condition = c(rep(0, valor / 2), rep(1, valor / 2))
  )
  
  # Simulación con los parámetros dinámicos
  sim_data <- simulate_dataset(
    vars_dinamico, pesos,
    mu_beta0, sigma_beta0,
    mean_sigma_epsilon, sd_sigma_epsilon,
    mean_sigma_inter, sd_sigma_inter,
    mean_sigma_intra, sd_sigma_intra,
    n_cells_ue_dinamico
  )
  
  sim_data <- build_differential_objects(sim_data)
  
  resultados <- run_differential_analysis(sim_data)
  
  metrics <- extract_performance_metrics(
    sim_data = sim_data,
    resultados = resultados,
    effect_name = "Delta_condition"
  )
  
  metrics$loop <- "n_ue"
  metrics$loop_value <- valor
      metrics$Fdr <- p.adjust(metrics$pvalue,method="BH")
  
  nombre_base <- paste0("n_ue_", valor)
  nombre_base <- gsub("\\.", "_", nombre_base)
  
  saveRDS(sim_data, paste0("results/loop_n_ue/", nombre_base, "_sim_data.rds"))
  saveRDS(resultados, paste0("results/loop_n_ue/", nombre_base, "_resultados.rds"))
  saveRDS(metrics, paste0("results/loop_n_ue/", nombre_base, "_metrics.rds"))
  write.csv(metrics, paste0("results/loop_n_ue/", nombre_base, "_metrics.csv"), row.names = FALSE)
  
  all_metrics_n_ue[[nombre_base]] <- metrics
}

all_metrics_n_ue <- do.call(rbind, all_metrics_n_ue)

saveRDS(all_metrics_n_ue, "results/loop_n_ue/all_metrics_n_ue.rds")
write.csv(all_metrics_n_ue, "results/loop_n_ue/all_metrics_n_ue.csv", row.names = FALSE)

head(all_metrics_n_ue)

}

###############################################################################
# LOOP 6: BARRIDO DEL NÚMERO DE CÉLULAS POR UE
###############################################################################

{
dir.create("results", showWarnings = FALSE)
dir.create("results/loop_n_cells", showWarnings = FALSE)

#valores_n_cells <- c(100, 200, 500, 1000, 5000)
valores_n_cells <- c(10,20,50,100, 200, 500)

all_metrics_n_cells <- list()

for (valor in valores_n_cells) {
  
  cat("Ejecutando células por UE =", valor, "\n")
  
  # 1. Generamos dinámicamente la cantidad de células por UE para las 20 UEs de referencia
  n_cells_ue_dinamico <- rep(valor, 20)
  
  sim_data <- simulate_dataset(
    vars, pesos, # Aquí sí usamos la "vars" de referencia (10 vs 10)
    mu_beta0, sigma_beta0,
    mean_sigma_epsilon, sd_sigma_epsilon,
    mean_sigma_inter, sd_sigma_inter,
    mean_sigma_intra, sd_sigma_intra,
    n_cells_ue_dinamico
  )
  
  sim_data <- build_differential_objects(sim_data)
  
  resultados <- run_differential_analysis(sim_data)
  
  metrics <- extract_performance_metrics(
    sim_data = sim_data,
    resultados = resultados,
    effect_name = "Delta_condition"
  )
  
  metrics$loop <- "n_cells_per_ue"
  metrics$loop_value <- valor
      metrics$Fdr <- p.adjust(metrics$pvalue,method="BH")
  
  nombre_base <- paste0("n_cells_", valor)
  nombre_base <- gsub("\\.", "_", nombre_base)
  
  saveRDS(sim_data, paste0("results/loop_n_cells/", nombre_base, "_sim_data.rds"))
  saveRDS(resultados, paste0("results/loop_n_cells/", nombre_base, "_resultados.rds"))
  saveRDS(metrics, paste0("results/loop_n_cells/", nombre_base, "_metrics.rds"))
  write.csv(metrics, paste0("results/loop_n_cells/", nombre_base, "_metrics.csv"), row.names = FALSE)
  
  all_metrics_n_cells[[nombre_base]] <- metrics
}
}

{
all_metrics_n_cells <- do.call(rbind, all_metrics_n_cells)

saveRDS(all_metrics_n_cells, "results/loop_n_cells/all_metrics_n_cells.rds")
write.csv(all_metrics_n_cells, "results/loop_n_cells/all_metrics_n_cells.csv", row.names = FALSE)

head(all_metrics_n_cells)
}

# 1. Aseguramos que exista una carpeta unificada para los gráficos
dir_plots <- "results/plots"
if (!dir.exists(dir_plots)) {
  dir.create(dir_plots, recursive = TRUE)
  cat("¡Carpeta creada con éxito en:", dir_plots, "!\n")
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

# Ejecutamos el pipeline gráfico completo
walk(config_loops, ~ generar_y_guardar_plots(.x$loop, .x$label, .x$folder, .x$suffix))
