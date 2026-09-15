###############################################################################
# File de funciones
###############################################################################

################################################################################
# RUTINA DE SIMULACIÓN
################################################################################

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
  
  #BIAS PLOT
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
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),
      axis.title = element_text(face = "bold", size = 15),
      axis.text = element_text(size = 13, colour = "black"),
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
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),
      axis.title = element_text(face = "bold", size = 15),
      axis.text = element_text(size = 13, colour = "black"),
      legend.position = "none"
      )#+ylim(-0.1,0.1)
 
  #std error plot
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
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),
      axis.title = element_text(face = "bold", size = 15),
      axis.text = element_text(size = 13, colour = "black"),
      legend.position = "none"
    )#+ylim(0,0.05)
    
    #FDR PLOT
    # Calculamos el número de genes con FDR < 0.05 para el cartelito
    counts_fdr <- metrics %>%
      group_by(loop_value) %>%
      summarise(
        total = n(),
        sig_count = sum(Fdr < 0.05, na.rm = TRUE),
        y_pos = pmax(max(-log10(Fdr), na.rm = TRUE) *0.7 , 2.2),
        .groups = "drop"
      )
    
    p_Fdr <- ggplot(metrics,aes(x = loop_value, y = -log10(Fdr))) +
      geom_violin(trim = FALSE, alpha = 0.6, colour = "black", linewidth = 0.4) +
      geom_boxplot(
        width = 0.1,
        outlier.shape = NA,
        fill = "white",
        colour = "black",
        alpha = 0.8
      ) +
      # Línea horizontal a 0.05 FDR
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.6) +
      # Cartelito con genes significativos (sig_count/total)
      geom_text(
        data = counts_fdr,
        aes(x = loop_value, y = y_pos, label = paste0(sig_count, "/", total)),
        nudge_x = 0.15,   # Desplaza la etiqueta hacia la derecha del centro del violín
        hjust = 0,  # Alinea el texto a la izquierda de esa posición (para que crezca hacia la derecha)
        vjust = 1,        # Ancla el texto por su parte superior para que crezca hacia abajo
        size = 3.3,
        fontface = "bold",
        inherit.aes = FALSE
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) + # 18% de aire extra arriba
      labs(
        title = paste("Influence of", label_x, "on statistical signal strength"),
        x = label_x,
        y = "-log10(Fdr) for condition effect"
      ) +
      theme_bw(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),
        axis.title = element_text(face = "bold", size = 15),
        axis.text = element_text(size = 13, colour = "black"),
        legend.position = "none"
      )#+ylim(0,0.05)

  

  
  # Guardamos en results/plots/
  
  dir.create(paste0(dir_plots, "/",loop_name),recursive = TRUE)
  ggsave(paste0(dir_plots,"/",loop_name, "/Plot_Bias.png"), plot = p_bias, width = 9, height = 5, dpi = 300)
  ggsave(paste0(dir_plots,"/",loop_name, "/Plot_Bias_fit.png"), plot = p_bias_fit, width = 9, height = 5, dpi = 300)
  ggsave(paste0(dir_plots,"/",loop_name, "/Plot_stde.png"), plot = p_stde, width = 9, height = 5, dpi = 300)
  ggsave(paste0(dir_plots,"/",loop_name, "/Plot_Fdr.png"), plot = p_Fdr, width = 9, height = 5, dpi = 300)

}

# 1) Matriz diseño experimental (D)
generateMatrix_D <- function(n_ue, vars){
  UE <- paste0("U", seq_len(n_ue)) # me crea secuencia (U1, U2, ....) 
  #comprobaciones
  if (!is.list(vars)) {
    stop("vars debe ser una lista.")
  }
  if (length(vars) == 0) {
    stop("vars no puede estar vacía.")
  }
  lengths_vars <- sapply(vars, length)
  if (any(lengths_vars != n_ue)) {
    stop("Todas las variables deben tener la misma longitud que UE.")
  }
  
  D <- data.frame(UE = UE, vars, check.names = FALSE)
  return(D)
}

# 2) Matriz de efectos fijos (E)
generateMatrix_E <- function(n_g, pesos) {
  G <- paste0("G", seq_len(n_g))
  
  if (!is.list(pesos)) {
    stop("pesos debe ser una lista.")
  }
  
  if (length(pesos) == 0) {
    stop("pesos no puede estar vacía.")
  }
  
  lengths_pesos <- sapply(pesos, length)
  if (any(lengths_pesos != n_g)) {
    stop("Todos los vectores de pesos deben tener la misma longitud que n_g.")
  }
  
  E <- data.frame(G = G, pesos, check.names = FALSE)
  return(E)
}

# 3) Array del término basal 
generate_beta0 <- function(n_g, mu_beta0, sigma_beta0) {
  
  beta0 <- rnorm(n =
                   n_g, mean = mu_beta0, sd = sigma_beta0)
  
  return(beta0)
}

# 4.1) Función auxiliar: convierte media y sd reales de una lognormal
#hemos cambiado el input para que sea mas intuitivo y se metan la media y desviación estándar 'reales'
#las quiero transformar a meanlog y sdlog, yo me sé las formulas de la lognormal, despejo mu, sigma^2
#explicado mejor en apuntes
#Genero función que me los transforma para usarlos
convert_lognormal_params <- function(mean_real, sd_real) {
  
  if (mean_real <= 0) {
    stop("La media real debe ser positiva para una distribución lognormal.")
  }
  
  if (sd_real < 0) {
    stop("La desviación típica real debe ser positiva para una distribución lognormal.")
  }
  
  sdlog <- sqrt(log(1 + (sd_real^2 / mean_real^2)))
  meanlog <- log(mean_real) - (sdlog^2 / 2)
  
  return(list(meanlog = meanlog, sdlog = sdlog))
}

# 4.2) Array del ruido global, de fondo
generate_sigma_epsilon <- function(n_g, mean_sigma_epsilon, sd_sigma_epsilon) {
  
  params <- convert_lognormal_params(
    mean_sigma_epsilon,
    sd_sigma_epsilon
  )
  
  sigma_epsilon <- rlnorm( 
    n = n_g,
    meanlog = params$meanlog,
    sdlog = params$sdlog
  )
  
  return(sigma_epsilon)
}

# 5) Delta inter: array con un sigma por gen

generate_sigma_inter <- function(n_g, mean_sigma_inter, sd_sigma_inter) {
  
  params <- convert_lognormal_params(
    mean_sigma_inter,
    sd_sigma_inter
  )
  
  sigma_inter <- rlnorm(
    n = n_g,
    meanlog = params$meanlog,
    sdlog = params$sdlog
  )
  
  return(sigma_inter)
}

# 6) Delta intra: matriz con un sigma por gen y por ue

generate_sigma_intra <- function(n_g, n_ue, mean_sigma_intra, sd_sigma_intra) {
  
  params <- convert_lognormal_params(
    mean_sigma_intra,
    sd_sigma_intra
  )
  
  sigma_intra <- matrix(
    rlnorm(
      n = n_g * n_ue,
      meanlog = params$meanlog,
      sdlog = params$sdlog
    ),
    nrow = n_g,
    ncol = n_ue
  )
  
  return(sigma_intra)
}

# 7) construimos termino basal
build_basal_term <- function(beta0, n_cells_total) {
  
  
  basal_term <- matrix(beta0,
                       nrow = length(beta0),
                       ncol = n_cells_total)
  
  return(basal_term)
}

# 8) construimos efectos fijos

build_fixed_effects_term <- function(D,E, n_cells_ue) {
  #convierto las tablas en matrices de numeros
  D_numeric <- as.matrix(D[, -1])
  E_numeric <- as.matrix(E[, -1])
  
  fixed_effects_by_ue <- D_numeric %*% t(E_numeric)
  fixed_effects_by_ue <- t(fixed_effects_by_ue)
  
  #ahora la tengo que expandir, para cada ue, repito esa columna, 
  #tantas veces como células tenga esa ue
  indices <- rep(1:length(n_cells_ue), times = n_cells_ue)
  fixed_effects <- fixed_effects_by_ue[, indices]
  
  return(fixed_effects)
}

# 9) construimos delta inter
#tengo un array de sigmas, una por gen. Voy a construir una distr normal
#para cada gen. Y voy a sacar de cada distribución tantos valores como ue tenga

build_inter_term <- function(sigma_inter, n_cells_ue, n_g, n_ue) {
  
  # VALIDACIÓN: Impedir sigmas negativas
  if (any(sigma_inter < 0)) stop("Existen valores negativos en sigma_inter")
  
  # SUSTITUCIÓN: sapply genera una fila de ruidos para cada sigma
  # t() es para trasponer y que los genes sigan siendo filas
  delta_inter_by_ue <- t(sapply(sigma_inter, function(s) {
    rnorm(n = n_ue, mean = 0, sd = s)
  }))
  
  indices <- rep(1:n_ue, times = n_cells_ue)
  inter_term <- delta_inter_by_ue[, indices]
  return(inter_term)
}
# 10) construimos delta intra
#ahora tengo una matriz de sigmas, una por gen y por ue
#voy a sacar tantos valores (de una normal de media 0), como células haya en esa ue
#es exactamente lo mismo que los bucles, pero vectorizado.

build_intra_term <- function(sigma_intra, n_cells_ue, n_g, n_ue, n_cells_total) {
  
  if (any(sigma_intra < 0)) stop("Existen valores negativos en sigma_intra")

  # 1.Comprobamos que la matriz de sigmas tenga el tamaño correcto
  if (nrow(sigma_intra) != n_g || ncol(sigma_intra) != n_ue) {
    stop("Error: Las dimensiones de sigma_intra deben ser [n_g x n_ue].")
  }

  # 2.Creamos un vector de índices para "estirar" la matriz de sigmas.
  # Si la UE1 tiene 10 células, repetimos la primera columna de sigmas 10 veces.
  indices_celulas <- rep(1:n_ue, times = n_cells_ue)
  sigma_expandida <- sigma_intra[, indices_celulas]
  
  # 3. GENERACIÓN VECTORIZADA
  # Generamos una matriz de ruido estándar de golpe
  ruido_blanco <- matrix(rnorm(n_g * n_cells_total), 
                         nrow = n_g, 
                         ncol = n_cells_total)
  
  # 4. APLICACIÓN DE SIGMAS
  # Multiplicamos el ruido por las sigmas correspondientes elemento a elemento.
  intra_term <- ruido_blanco * sigma_expandida
  
  return(intra_term)
}

# 11) construimos epsilon 
#ahora tengo un array con un sigma por gen. No separo por ue, directamente
# saco tantos valores como células tenga para cada sigma del gen i (es ruido de fondo)

build_epsilon_term <- function(sigma_epsilon, n_cells_total, n_g) {
  
  if (any(sigma_epsilon < 0)) stop("Existen valores negativos en sigma_epsilon")
  
  # SUSTITUCIÓN: Aplicamos rnorm a cada valor del vector sigma_epsilon
  epsilon_term <- t(sapply(sigma_epsilon, function(s) {
    rnorm(n = n_cells_total, mean = 0, sd = s)
  }))
  
  return(epsilon_term)
}

#Con todo ya construido, vamos a generar la función que me da la matriz de expresión génica

simulate_dataset <- function(vars, pesos, 
                             mu_beta0, sigma_beta0,
                             mean_sigma_epsilon, sd_sigma_epsilon,
                             mean_sigma_inter, sd_sigma_inter,
                             mean_sigma_intra, sd_sigma_intra,
                             n_cells_ue) {
  
  n_ue <- length(n_cells_ue)
  n_g  <- length(pesos[[1]]) # Miramos el primer vector de la lista de pesos
  n_cells_total <- sum(n_cells_ue)
  
  # A) generamos inputs generales
  D <- generateMatrix_D(n_ue, vars)
  E <- generateMatrix_E(n_g, pesos)
  beta0 <- generate_beta0(n_g, mu_beta0, sigma_beta0)
  sigma_epsilon <- generate_sigma_epsilon(n_g, mean_sigma_epsilon, sd_sigma_epsilon)
  sigma_inter <- generate_sigma_inter(n_g, mean_sigma_inter, sd_sigma_inter)
  sigma_intra <- generate_sigma_intra(n_g, n_ue, mean_sigma_intra, sd_sigma_intra)
  
  # B) construimos las matrices que se suman
  basal_term <- build_basal_term(beta0, n_cells_total)
  fixed_term <- build_fixed_effects_term(D, E, n_cells_ue)
  inter_term <- build_inter_term(sigma_inter, n_cells_ue, n_g, n_ue)
  intra_term <- build_intra_term(sigma_intra, n_cells_ue, n_g, n_ue, n_cells_total)
  epsilon_term <- build_epsilon_term(sigma_epsilon, n_cells_total, n_g) 
  
  # C) construimos el data set
  # Eg,c = beta0 + Efectos + Delta_Inter + Delta_Intra + Epsilon
  matrix_expresion_genica <- basal_term + fixed_term + inter_term + intra_term + epsilon_term
  
  #ponemos nombres a las filas (genes) columnas (células) de la matriz:
  # nombres genes
  gene_names <- paste0("G", seq_len(n_g))
  
  # nombres células
  cell_names <- unlist(lapply(seq_len(n_ue), function(u) {
    paste0("U", u, "_C", seq_len(n_cells_ue[u]))
  }))
  
  rownames(matrix_expresion_genica) <- gene_names
  colnames(matrix_expresion_genica) <- cell_names
  
  # D) METADATOS
  cell_metadata <- data.frame(
    cell_id = paste0("C", seq_len(n_cells_total)),
    UE = rep(paste0("U", seq_len(n_ue)), times = n_cells_ue),
    row.names = cell_names
  )
  
  ue_metadata <- D
  rownames(ue_metadata) <- ue_metadata$UE
  
  return(list(matrix_expresion_genica = matrix_expresion_genica, 
              cell_metadata = cell_metadata,
              ue_metadata = ue_metadata,
              sigma_epsilon = sigma_epsilon,
              sigma_inter = sigma_inter,
              sigma_intra = sigma_intra))
}
#OUTPUTS: la matriz final Eg,c y la matriz 'diccionario' para saber a qué ue corresponde cada célula y la matrix D

################################################################################
# RUTINA DE INFERENCIA
################################################################################
build_differential_objects <- function(sim_data){
  
  # 1. Extraemos los nombres de las Unidades Experimentales (U1, U2...)
  # Los sacamos de ue_metadata
  ue_names <- sim_data$ue_metadata$UE #con el tercer $ elegimos la columna
  n_ue <- length(ue_names)
  
  # 2. Generamos todas las combinaciones posibles de 2 en 2
  # Esto nos da una matriz de 2 filas y donde cada columna es una pareja (U1-U2, U1-U3...) 
  parejas <- combn(ue_names, 2) 
  
  # 3. Creamos etiquetas para las columnas (ej: "U1-U2")
  #contamos el número de columnas que tenemos, es decir, el número de parejas
  #en realidad es lo mismo que hacer (n_ue*(n_ue-1))/2 
  n_parejas <- ncol(parejas) 
  nombres_parejas <- apply(parejas, 2, function(x) paste0(x[1], "-", x[2])) # me devuelve un vector con todas las parejas
  
  # 4. Preparamos las 3 matrices que necesitamos
    #matriz 1: Ncells,i: indica cuántas células tiene el primer ue de la pareja
    #matriz 2: Ncells,j: indica cuántas células tiene el segundo ue de la pareja
    #matriz 3: Ncells,i>j: para cada gen, tenemos que contar cuántas células del individuo i 
                          #tienen mas expresión que las ceúlas del individuo j
    #esta última la acabamos llamando Victorias (i>j)
    # y tmb he añadido la contraria: Derrotas (i<j)
    #cada una de las tres matrices tiene tamaño n_g*n_parejas
  
  #calculamos n_g
  n_g <- nrow(sim_data$matrix_expresion_genica)
  nombres_genes <- rownames(sim_data$matrix_expresion_genica) #nombres de G1, G2, ....
  
  # De momento las creamos vacías para luego rellenarlas
  #para crear un nuevo elemento en la lista de sim_data usas el $ y añades el nombre
  
  sim_data$Ncells_i <- matrix(0, nrow = n_g, ncol = n_parejas, dimnames = list(nombres_genes, nombres_parejas))
  sim_data$Ncells_j <- matrix(0, nrow = n_g, ncol = n_parejas, dimnames = list(nombres_genes, nombres_parejas))
  sim_data$Victorias_ij <- matrix(0, nrow = n_g, ncol = n_parejas, dimnames = list(nombres_genes, nombres_parejas))
  
  # Creamos matriz derrotas tmb 
  sim_data$Derrotas_ij <- matrix(0, nrow = n_g, ncol = n_parejas, dimnames = list(nombres_genes, nombres_parejas))
 
   # 5. Usamos 'apply' para recorrer las columnas de la matriz 'parejas' de golpe en lugar del for que había oensado
  # Cada 'p' será una columna (una pareja de nombres)
  
  resultados_lista <- apply(parejas, 2, function(p) {
    
    # A) Identificamos quiénes son i y j
    nombre_i <- p[1] #nombre de la primera pareja
    nombre_j <- p[2] #nombre de la segunda pareja
    
    # B) Filtramos las celdas (índices) usando los metadatos
    #en cell_metadata yo tenía en una columna, todas las células y en la otra la ue a la que se corresponden
    #lo que hace esto es comparar toda la matriz con el nombre_i que tengo, cuando estamos como en el primer 'buvle'
    # tendremos ue1, asi que devuelve true pues para las x primeras células qeu sea del ue1, y luego el which me pone 
    #las posiciones donde es true: ES DECIR, vemos en qué posiciones del cell_metadata estanos ahora
    cols_i <- which(sim_data$cell_metadata$UE == nombre_i)
    cols_j <- which(sim_data$cell_metadata$UE == nombre_j)
    
    # C) Extraemos la expresión de esos dos grupos
    #ya sabemos qué celulas queremos ahora comparar
    #vamos a la matrzi_expresion_génica que tiene todas las células en columnas, pero ahora con
    #cols_i y cols_j ya sabemos cuales tenemos que comparar
    #lo que hace esto es guardar esos valores de expresión génica de las celulas que me interesan en cada caso
    exp_i <- sim_data$matrix_expresion_genica[, cols_i, drop = FALSE]
    exp_j <- sim_data$matrix_expresion_genica[, cols_j, drop = FALSE]
    
    # D) Hacemos la comparación 'todos contra todos' por cada gen
    # Esto devuelve un vector con el número de veces que la expresión génica de (i>j) de cada gen para esta pareja
    #Me sirve para Victorias
    counts_v <- sapply(1:n_g, function(g) {
      sum(outer(exp_i[g, ], exp_j[g, ], ">"))
    })
    
    # Derrotas (Total - Victorias)
    # n_i * n_j es el número total de comparaciones por gen para esa pareja
    total_por_gen <- length(cols_i) * length(cols_j)
    counts_d <- total_por_gen - counts_v
    
    # Devolvemos un paquetito con los 4 datos de esta pareja
    return(list(
      ni = length(cols_i), # cuántas células tiene  el primer elemento de la pareja
      nj = length(cols_j), # cuántas células tiene  el segundo elemento de la pareja
      comp_v = counts_v, # cuántas veces i>j Victorias
      comp_d = counts_d
    ))
  })
  
  vec_ni <- sapply(resultados_lista, `[[`, "ni")
  vec_nj <- sapply(resultados_lista, `[[`, "nj")
  
  #Rellenamos las matrices
  sim_data$Ncells_i <- matrix(rep(vec_ni, each = n_g), nrow = n_g, dimnames = list(nombres_genes, nombres_parejas))
  sim_data$Ncells_j <- matrix(rep(vec_nj, each = n_g), nrow = n_g, dimnames = list(nombres_genes, nombres_parejas))
  
  #Ahora las de Victorias y derrotas
  res_v <- sapply(resultados_lista, `[[`, "comp_v")
  res_d <- sapply(resultados_lista, `[[`, "comp_d")
  
  # Si solo hay un gen, sapply devuelve vector, lo forzamos a matriz (me daba error sino)
  if (n_g == 1) {
    res_v <- matrix(res_v, nrow = 1)
    res_d <- matrix(res_d, nrow = 1)
  }
  
  sim_data$Victorias_ij <- res_v
  sim_data$Derrotas_ij  <- res_d

  
  #por si acaso vuelvo a poner los nombres
  rownames(sim_data$Victorias_ij) <- nombres_genes
  rownames(sim_data$Derrotas_ij)  <- nombres_genes
  colnames(sim_data$Victorias_ij) <- nombres_parejas
  colnames(sim_data$Derrotas_ij)  <- nombres_parejas
  
  
  # 7. Ahora necesitamos una cuarta matrix: la matriz de deltas de cada vars
  # n_parejas*num_vars
  # Construcción de la Matriz de Diseño Diferencial (Deltas)
  # Extraemos solo las columnas numéricas del diseño (quitamos la columna "UE")
  diseno_original <- sim_data$ue_metadata[, -1, drop = FALSE]
  
  # Usamos apply para recorrer las parejas y restar sus atributos
  delta_matrix <- apply(parejas, 2, function(p) {
    # Buscamos la fila del individuo i y del individuo j en ue_metadata
    fila_i <- which(sim_data$ue_metadata$UE == p[1])
    fila_j <- which(sim_data$ue_metadata$UE == p[2])
    
    # Restamos los valores: Trait(i) - Trait(j)
    return(as.numeric(diseno_original[fila_i, ] - diseno_original[fila_j, ]))
  })
  
  #Si solo hay una variable (ej. condition), 
  # apply devuelve vector, lo forzamos a matriz antes de trasponer
  if (is.null(dim(delta_matrix))) {
    delta_matrix <- matrix(delta_matrix, nrow = 1)
  }
  
  # Trasponemos para que las filas sean las parejas y las columnas los Traits
  delta_matrix <- t(delta_matrix)
  
  # Ponemos nombres claros a las columnas (ej: "Delta_condition")
  colnames(delta_matrix) <- paste0("Delta_", colnames(diseno_original))
  rownames(delta_matrix) <- nombres_parejas
  
  # Lo guardamos en el objeto final como un data.frame para que sea fácil de leer
  sim_data$Design_deltas <- as.data.frame(delta_matrix)
  
  return(sim_data)
  
}
#Para añadir dos columnas finales al output, necesitamos pasar a la escala de P
# Función auxiliar: inverse logit
# Transforma valores de escala logit a escala probabilidad
inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}
logit <- function(p) {
  log(p/(1-p))
}

run_differential_analysis <- function(sim_data, start = 1) {
  library(glmmTMB)
  # 1. Recordamos genes y sus nombres
  n_g <- nrow(sim_data$matrix_expresion_genica)
  gene_names <- rownames(sim_data$matrix_expresion_genica)
  # 2. Extraemos el diseño de las comparaciones (las Deltas)
  design <- sim_data$Design_deltas
  # 3. Creamos dos listas para guardar los resultados
  lista_resultados <- vector("list", n_g)
  lista_fits <- vector("list", n_g)
  # 4. Ajustamos el modelo para cada gen
  for (g in seq_len(n_g)) {
    # Preparamos el dataframe para EL GEN ACTUAL
    df_gen <- data.frame( wins = sim_data$Victorias_ij[g, ], losses = sim_data$Derrotas_ij[g, ], design )
    # Ajustamos el modelo Beta-Binomial
    nombres_deltas <- colnames(design)
    fit <- try(
      #glmmTMB(
      #formula = as.formula( paste( "cbind(wins, losses) ~ 0+", paste(nombres_deltas, collapse = " + ") ) ),
      #family = betabinomial(link = "logit"),
      #data = df_gen, start = list(beta = start)
      #),
      glmmTMB(
        formula = as.formula(
          paste(
            "cbind(wins, losses) ~ 0 +",
            paste(nombres_deltas, collapse = " + ")
          )
        ),
        dispformula = as.formula(
          paste(
            "~",
            paste(nombres_deltas, collapse = " + ")
          )
        ),
        family = betabinomial(link = "logit"),
        data = df_gen,
        start = list(beta = start),
        control = glmmTMBControl(
          optCtrl = list(
            step.max = 0.1,
            rel.tol = 1e-10,
            x.tol = 1e-10,
            iter.max = 100000,
            eval.max = 100000
          )
        )
      ),
      silent = TRUE )
    
    # Guardamos el objeto fit
    if (!inherits(fit, "try-error")) {
      lista_fits[[g]] <- fit
      # Extraemos la tabla de coeficientes
      tabla_coef <- summary(fit)$coefficients$cond
      tabla_coef <- as.data.frame(tabla_coef)
      tabla_coef$Estimate=fit$fit$par["beta"]
      tabla_coef$Probability <- inv_logit(tabla_coef$Estimate)
      tabla_coef$CI_low <- inv_logit( tabla_coef$Estimate - 1.96 * tabla_coef$`Std. Error` )
      tabla_coef$CI_high <- inv_logit( tabla_coef$Estimate + 1.96 * tabla_coef$`Std. Error` )
      lista_resultados[[g]] <- tabla_coef
    } else {
      lista_fits[[g]] <- NULL
      lista_resultados[[g]] <- NULL
    }
  }
  # Ponemos nombre a cada elemento para saber qué gen es
  names(lista_resultados) <- gene_names
  names(lista_fits) <- gene_names
  # Devolvemos ambos resultados
  return(list( resultados = lista_resultados, fits = lista_fits ))
}
################################################################################
# RUTINA DE VISUALIZACIÓN
################################################################################
library(ggplot2)

extract_performance_metrics <- function(sim_data, resultados,effect_name = "Delta_condition") {
  
  fits=resultados$fits
  resultados=resultados$resultados
  
  design <- sim_data$Design_deltas
  nombres_deltas <- colnames(design)
  
  formula_design <- as.formula(
    paste("~0+", paste(nombres_deltas, collapse = " + "))
  )
  
  X <- model.matrix(formula_design, data = design)
  
  gene_names <- rownames(sim_data$matrix_expresion_genica)
  
  metrics_list <- lapply(gene_names, function(gene) {
    
    tabla <- resultados[[gene]]
    fit=fits[[gene]]
    
    if (is.null(tabla)) {
      return(NULL)
    }
    
    if (!(effect_name %in% rownames(tabla))) {
      return(NULL)
    }
    
    {
      if(TRUE){
        w <- fit$obj$env$data$yobs
        n <- fit$obj$env$data$size
        l <- n - w
        
        ## Matriz de diseño para la media
        X <- fit$obj$env$data$X
        
        ## Parametros
        beta_hat  <- fixef(fit)$cond
        gamma_hat <- fixef(fit)$disp
        
        ## Media
        eta <- as.numeric(X %*% beta_hat)
        mu <- inv_logit(eta)
        
        ## Dispersion
        delta <- X[, "Delta_condition"]
        
        zeta <- gamma_hat["(Intercept)"] +
          gamma_hat["Delta_condition"] * delta
        
        phi <- exp(zeta)
        
        ## Parametros de la Beta
        alpha <- mu * phi
        beta_BB <- (1 - mu) * phi
        
        ## Score respecto a eta
        score_eta <- phi * mu * (1 - mu) *
          (
            digamma(w + alpha) - digamma(alpha) -
              digamma(l + beta_BB) + digamma(beta_BB)
          )
        
        ## Score respecto a beta
        grad_beta <- colSums(X * score_eta)
        
        bias_fit <- grad_beta
      }
      if(FALSE){
        w <- fit$obj$env$data$yobs
        n <- fit$obj$env$data$size
        l <- n - w
        
        X <- fit$obj$env$data$X
        
        ## Parámetros
        p <- ncol(X)
        
        beta_hat <- fit$fit$par[1:p]
        gamma_hat <- fit$fit$par[(p + 1):(2 * p)]
        
        ## Media
        eta <- as.numeric(X %*% beta_hat)
        mu <- inv_logit(eta)
        
        ## Dispersión
        zeta <- as.numeric(X %*% gamma_hat)
        phi <- exp(zeta)
        
        ## Parámetros de la Beta
        alpha <- mu * phi
        beta_BB <- (1 - mu) * phi
        
        ## Score respecto a eta
        score_eta <- phi * mu * (1 - mu) *
          (digamma(w + alpha) - digamma(alpha) -
             digamma(l + beta_BB) + digamma(beta_BB))
        
        ## Score respecto a beta
        score_beta <- X * score_eta
        grad_beta <- colSums(score_beta)
        
        ## Score respecto a phi
        score_phi <- mu *
          (digamma(w + alpha) - digamma(alpha)) +
          (1 - mu) *
          (digamma(l + beta_BB) - digamma(beta_BB)) -
          digamma(n + phi) +
          digamma(phi)
        
        ## Score respecto a log(phi)
        score_zeta <- phi * score_phi
        
        ## Score respecto a gamma
        score_gamma <- X * score_zeta
        grad_gamma <- colSums(score_gamma)
        
        ## Gradiente completo de la log-likelihood
        grad_score <- c(grad_beta, grad_gamma)
        
        bias_fit=grad_beta
        
      }
      if(FALSE){
        bias_fit=fit$obj$gr(fit$fit$par)[1]
      }
      
    }
    
    pred_frac <- inv_logit(eta)
    obs_frac <- w / n
    error_pairwise <- (pred_frac - obs_frac)
    bias_observaciones=mean(error_pairwise)
    
    sigma_vector <- (pred_frac - obs_frac)^2
    std_error=sqrt(sum(sigma_vector))/length(sigma_vector)
    
    
    pvalue <- tabla[effect_name, "Pr(>|z|)"]
    
    data.frame(
      gene = gene,
      estimate = tabla[effect_name, "Estimate"],
      std_error = tabla[effect_name, "Std. Error"],
      pvalue = pvalue,
      minus_log10_pvalue = -log10(pvalue),
      bias_fit = bias_fit,
      bias_win_fraction = bias_observaciones,
      sd_error_win_fraction = std_error
    )
  })
  
  do.call(rbind, metrics_list)
}
plot_raw_gene_1D <- function(sim_data, gene, variable) {
  
  if (!(gene %in% rownames(sim_data$matrix_expresion_genica))) {
    stop("El gen indicado no existe en matrix_expresion_genica.")
  }
  
  if (!(variable %in% colnames(sim_data$ue_metadata))) {
    stop(paste("La variable", variable, "no existe en ue_metadata. Variables disponibles:",
               paste(colnames(sim_data$ue_metadata), collapse = ", ")))
  }
  
  library(ggplot2)
  
  # Expresión del gen elegido
  #coge toda la fila de la matriz: matrix_expresion_genica y devuelve un array de valores de todas las expresiones
  expr <- as.numeric(sim_data$matrix_expresion_genica[gene, ])
  
  # Metadata celular
  #pegamos a la matriz metadata una columna mas que contenga la expresion
  df <- sim_data$cell_metadata
  df$expression <- expr
  
  # Añadimos la variable experimental desde ue_metadata
  ue_info <- sim_data$ue_metadata[, c("UE", variable)]
  colnames(ue_info) <- c("UE", "variable_value")
  
  #matriz con: cell_id   UE   expression   condition
  df <- merge(df, ue_info, by = "UE")
  
  # Detectamos si la variable es categórica o continua
  n_unique <- length(unique(df$variable_value))
  
  if (n_unique <= 5) {
    
    # Caso categórico/binario: violín + jitter
    df$variable_value <- as.factor(df$variable_value)
    
    p <- ggplot(df, aes(x = variable_value, y = expression)) +
      geom_violin(trim = FALSE) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
      labs(
        title = paste("Raw expression of", gene, "by", variable),
        x = variable,
        y = "Gene expression"
      ) +
      theme_minimal()
    
  } else {
    
    # Caso continuo: puntos + tendencia
    df$variable_value <- as.numeric(df$variable_value)
    
    p <- ggplot(df, aes(x = variable_value, y = expression)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(
        title = paste("Raw expression of", gene, "by", variable),
        x = variable,
        y = "Gene expression"
      ) +
      theme_minimal()
  }
  
  return(p)
  }
plot_raw_gene_1D_by_UE <- function(sim_data, gene, variable) {
  
  library(ggplot2)
  
  if (!(gene %in% rownames(sim_data$matrix_expresion_genica))) {
    stop("El gen indicado no existe en matrix_expresion_genica.")
  }
  
  if (!(variable %in% colnames(sim_data$ue_metadata))) {
    stop(paste("La variable", variable, "no existe en ue_metadata. Variables disponibles:",
               paste(colnames(sim_data$ue_metadata), collapse = ", ")))
  }
  
  # Expresión del gen elegido
  expr <- as.numeric(sim_data$matrix_expresion_genica[gene, ])
  
  # Metadata celular
  df <- sim_data$cell_metadata
  df$expression <- expr
  
  # Añadimos la variable experimental desde ue_metadata
  ue_info <- sim_data$ue_metadata[, c("UE", variable)]
  colnames(ue_info) <- c("UE", "variable_value")
  
  df <- merge(df, ue_info, by = "UE")
  
  # Ordenamos las UEs según el valor de la variable
  ue_order <- unique(ue_info$UE[order(ue_info$variable_value)])
  df$UE <- factor(df$UE, levels = ue_order)
  
  # Para colorear por condición/variable
  df$variable_value <- as.factor(df$variable_value)
  
  p <- ggplot(df, aes(x = UE, y = expression, fill = variable_value)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_jitter(
      aes(color = variable_value),
      width = 0.12,
      alpha = 0.4,
      size = 0.8
    ) +
    labs(
      title = paste("Raw expression of", gene, "by experimental unit"),
      subtitle = paste("Grouped by", variable),
      x = "Experimental unit",
      y = "Gene expression",
      fill = variable,
      color = variable
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

plot_gene_by_UE <- function(sim_data, gene) {
    library(ggplot2) # Comprobar que el gen existe
    if (!(gene %in% rownames(sim_data$matrix_expresion_genica))) { stop(paste("El gen", gene, "no existe.")) } # Extraemos la expresión del gen
    expr <- as.numeric(sim_data$matrix_expresion_genica[gene, ]) # Nombres de las células
    cell_names <- colnames(sim_data$matrix_expresion_genica) # Extraemos la UE de cada célula
    UE <- sub("_C.*", "", cell_names) # Dataframe para ggplot
    df <- data.frame( UE = UE, expression = expr ) # Ordenamos las UEs numéricamente
    df$UE <- factor( df$UE, levels = paste0("U", seq_len(length(unique(UE)))) ) # Violin plot
    p <- ggplot(df, aes(x = UE, y = expression)) +
    geom_violin(trim = FALSE) +
    labs( title = paste("Expression of", gene), x = "Experimental unit", y = "Gene expression" ) +
    theme_bw(base_size = 12)
    return(p) }
