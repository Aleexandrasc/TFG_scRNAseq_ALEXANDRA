# TFG_scRNAseq_ALEXANDRA
Pipeline en R para el análisis de expresión diferencial en scRNA-seq mediante comparaciones pareadas célula a célula y regresión Beta-Binomial (`glmmTMB`).


## Estructura del Repositorio

* **'06_functions.R'**: Contiene todas las funciones necesarias para el proyecto:
  * Simulación del dataset de expresión génica con estructura jerárquica de ruido ('simulate_dataset').
  * Construcción de las matrices de diseño diferencial ($\Delta X$) y recuento de victorias/derrotas célula a célula por parejas de UEs ('build_differential_objects').
  * Ajuste gen a gen mediante el modelo Beta-Binomial con dispersión estimada ('run_differential_analysis').
  * Extracción y cálculo de métricas de rendimiento (sesgo en predicciones, error estándar, FDR) junto con funciones de visualización.

* **'06_main.R'**: Script principal de ejecución. Define los parámetros basales del experimento y realiza los barridos sistemáticos de sensibilidad guardando los resultados y métricas:
  * Ruido biológico inter-individual ($\sigma_{\text{inter}}$).
  * Variabilidad intra-individual celular ($\sigma_{\text{intra}}$).
  * Ruido de fondo global ($\sigma_\epsilon$).
  * Magnitud del efecto biológico ($\beta$).
  * Número de unidades experimentales ($N_{UE}$) y número de células por UE ($N_{\text{cells}}$).

* **'06_plots.R'**: Script para automatizar la generación de gráficos globales a partir de los datos guardados en los barridos (gráficos de violín y cajas para sesgo, error estándar y $(-\log_{10}(\text{FDR}))$.

* **'plot_expresion_1gen.R'**: Script para generar los paneles detallados de inspección sobre un gen concreto (como $G_1$):
  * Distribución de la expresión celular por unidad experimental agrupada por condición (Control vs. Tratado).
  * Curvas de densidad teórica de la distribución Beta estimada para comparaciones dentro del mismo grupo vs. entre grupos.
