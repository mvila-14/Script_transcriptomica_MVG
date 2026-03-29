# Script_transcriptomica_MVG

# Análisis Transcriptómico de Tumores Paratiroideos (ERα-negativos)

Este repositorio contiene el flujo de trabajo bioinformático para el estudio de la respuesta temprana (24h) a moduladores estrogénicos (DPN y OHT) en cultivos primarios de tumores paratiroideos.

# Contenido del Repositorio
* `script_final.R`: Análisis estadístico completo en R (DESeq2, PCA, Volcano Plots y GSEA).
* `DPN_response.gmt`: Firma molecular de respuesta a DPN (48h) utilizada para el análisis de enriquecimiento.
* `metadata.tsv`: Metadatos de las 24 muestras (Paciente, Agente y Tiempo).

# Pipeline Bioinformático
1. **Preprocesamiento (Bash):** Control de calidad con `FastQC`, alineamiento contra el cromosoma 21 con `HISAT2` y cuantificación con `HTSeq-count`.
2. **Análisis Diferencial (R):** Uso de `DESeq2` con un modelo multifactorial `~ patient + time + agent` para corregir la variabilidad inter-paciente.
3. **Enriquecimiento Funcional:** Análisis GSEA con el paquete `fgsea` para validar la respuesta coordinada a las 24h.

# Observaciones principales
* Identificación de RASD1 como gen clave en la respuesta a OHT a las 24h.
* Identificación de la caída de MMP3 tras el tratamiento con DPN.
* Validación mediante GSEA de que la firma de 48h ya está significativamente enriquecida a las 24h (NES = 1.61, p-adj = 0.001).
