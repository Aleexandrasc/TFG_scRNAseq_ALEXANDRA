# TFG_scRNAseq_ALEXANDRA

Pipeline en R para el análisis de expresión diferencial en datos de secuenciación de ARN de célula única (scRNA-seq) mediante comparaciones pareadas célula a célula y modelado estadístico con regresión Beta-Binomial (`glmmTMB`).



## Estructura del Repositorio

* **`07_functions.R`**: Contiene el conjunto de funciones modulares que articulan el pipeline:
  * Simulación de perfiles de expresión génica mediante una estructura jerárquica de ruido biológico y técnico (`simulate_dataset`).
  * Construcción de matrices de diseño diferencial ($\Delta X$) y cómputo de victorias/derrotas célula a célula entre pares de unidades experimentales (`build_differential_objects`).
  * Ajuste estadístico gen a gen utilizando un modelo lineal generalizado Beta-Binomial (`run_differential_analysis`).
  * Evaluación de rendimiento inferencial: cálculo analítico de la función *score* / gradientes ($\nabla_{\boldsymbol{\beta}} \ell$), sesgo en probabilidades, error estándar y corrección por tasa de falso descubrimiento (FDR) (`extract_performance_metrics`).

* **`07_main.R`**: Script principal de ejecución que orquesta los análisis. Define los parámetros basales del experimento y ejecuta los barridos sistemáticos de sensibilidad paramétrica:
  * Ruido biológico inter-individual ($\sigma_{\text{inter}}$).
  * Variabilidad celular intra-individual ($\sigma_{\text{intra}}$).
  * Ruido estocástico de fondo ($\sigma_\epsilon$).
  * Magnitud del efecto biológico del tratamiento ($\beta$).
  * Cohorte de réplicas biológicas ($N_{\text{UE}}$) y profundidad de muestreo celular ($N_{\text{cells}}$).

* **`plot_expresion_1gen.R`**: Rutina complementaria orientada al diagnóstico visual e inspección detallada de genes individuales:
  * Distribución empírica de la expresión celular por unidad experimental desglosada por condición (Control frente a Tratado).
  * Curvas de densidad teórica de la distribución Beta estimada para pares intra-condición frente a pares inter-condición.
