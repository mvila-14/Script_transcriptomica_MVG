
#########################################################################
##          SCRIPT FINAL          ##
#########################################################################

# -----------------------------------------------------------------------
# PREGUNTA 4: DEGs (OHT vs CONTROL y DPN vs CONTROL)
# -----------------------------------------------------------------------

# 1. LIBRERÍAS (Todas las necesarias)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(clusterProfiler)
library(mygene)
library(pheatmap)
library(RColorBrewer)
library(enrichplot)

# 2. CARGA DE DATOS
counts <- read.table("rawcounts.tsv", header=TRUE, row.names=1, check.names=FALSE)
metadata <- read.table("metadata.tsv", header=TRUE, row.names=1, check.names=FALSE)

# Sincronización automática de columnas y filas 
counts <- counts[, rownames(metadata)]

# Convertir variables a factores 
metadata$patient <- factor(metadata$patient)
metadata$agent   <- factor(metadata$agent, levels = c("Control", "DPN", "OHT"))
metadata$time    <- factor(metadata$time)

# 3. OBJETO DESeq2 GLOBAL (24 muestras para máxima potencia estadística)
# El diseño bloquea el ruido del paciente y del tiempo para ver el fármaco
dds <- DESeqDataSetFromMatrix(countData = round(counts), 
                              colData = metadata, 
                              design = ~ patient + time + agent)

# Filtrado de genes de baja expresión
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 4. CONTROL DE CALIDAD: PCA 
# Usamos VST para normalizar los datos antes del gráfico
vsd <- vst(dds, blind=TRUE)

# Generamos el gráfico PCA agrupando por Agente y Tiempo 
plotPCA(vsd, intgroup=c("agent", "time")) + 
  theme_bw() + 
  geom_point(size=4) +
  ggtitle("PCA: Agrupamiento por Tratamiento y Tiempo (Global)")

# 5. ANÁLISIS DE EXPRESIÓN DIFERENCIAL
dds <- DESeq(dds)

# -----------------------------------------------------------------------
# PREGUNTA 4: DEGs (OHT vs Control y DPN vs Control)
# -----------------------------------------------------------------------

# Extraer resultados
res_OHT <- results(dds, contrast=c("agent","OHT","Control"), alpha=0.05)
res_DPN <- results(dds, contrast=c("agent","DPN","Control"), alpha=0.05)

# ANOTACIÓN: Traducir IDs de Ensembl a nombres de genes (Symbols)
genes_info <- queryMany(rownames(res_OHT), scopes="ensembl.gene", fields="symbol", species="human")
genes_info <- genes_info[!duplicated(genes_info$query), ]
symbol_map <- genes_info$symbol[match(rownames(res_OHT), genes_info$query)]

# Crear DataFrames anotados
res_OHT_df <- as.data.frame(res_OHT)
res_OHT_df$symbol <- symbol_map
res_DPN_df <- as.data.frame(res_DPN)
res_DPN_df$symbol <- symbol_map

# FILTRADO Y CONTEO (Usamos Log2FC > 0.5 para que salgan resultados significativos)
OHT_sig <- subset(res_OHT_df, padj < 0.05 & abs(log2FoldChange) > 0.5)
DPN_sig <- subset(res_DPN_df, padj < 0.05 & abs(log2FoldChange) > 0.5)

cat("\n==============================================\n")
cat("   RESULTADOS DE GENES SIGNIFICATIVOS (P4)    \n")
cat("==============================================\n")
cat("Nº genes significativos OHT vs Control:", nrow(OHT_sig), "\n")
cat("Nº genes significativos DPN vs Control:", nrow(DPN_sig), "\n")
cat("==============================================\n")

# Guardar tablas en CSV
write.csv(OHT_sig, "Resultados_OHT_Completo.csv")
write.csv(DPN_sig, "Resultados_DPN_Completo.csv")

# Volcanos Anotados
EnhancedVolcano(res_OHT_df, lab = res_OHT_df$symbol, x = "log2FoldChange", y = "padj", 
                pCutoff = 0.05, FCcutoff = 0.5, title = "OHT vs Control")

EnhancedVolcano(res_DPN_df, lab = res_DPN_df$symbol, x = "log2FoldChange", y = "padj", 
                pCutoff = 0.05, FCcutoff = 0.5, title = "DPN vs Control")

# -----------------------------------------------------------------------
# PREGUNTA 5: GSEA PARA DPN
# -----------------------------------------------------------------------

# 1. LIBRERÍAS
library(DESeq2)
library(fgsea)
library(ggplot2)

# 2. REPARACIÓN DEL OBJETO DDS (Para que reconozca el factor 'group')
# Creamos la columna combinada en el metadata
metadata$group <- factor(paste0(metadata$agent, "_", metadata$time))

# Volvemos a construir y ejecutar el dds con el nuevo diseño
# Usamos ~patient + group para bloquear el efecto individual de los pacientes
dds_fgsea <- DESeqDataSetFromMatrix(countData = round(counts), 
                                    colData = metadata, 
                                    design = ~ patient + group)

dds_fgsea <- DESeq(dds_fgsea[rowSums(counts(dds_fgsea)) >= 10, ])

# 3. EXTRACCIÓN DE RESULTADOS 
res_DPN_24h <- results(dds_fgsea, contrast=c("group", "DPN_24h", "Control_24h"))

# 4. PREPARACIÓN DEL RANKING
ranks_DPN <- res_DPN_24h$stat
names(ranks_DPN) <- rownames(res_DPN_24h)
ranks_DPN <- na.omit(ranks_DPN)
ranks_DPN <- sort(ranks_DPN, decreasing = TRUE)

# 5. CARGA DE RUTAS Y EJECUCIÓN DE fgsea
pathways_DPN <- gmtPathways("DPN_response.gmt")

set.seed(123) 
fgseaRes <- fgsea(pathways = pathways_DPN, 
                  stats = ranks_DPN,
                  minSize = 10,
                  maxSize = 500)

# 6. TABLA DE RESULTADOS
tabla_gsea <- as.data.frame(fgseaRes[, c("pathway", "pval", "padj", "NES", "size")])
print("--- TABLA DE RESULTADOS GSEA (P5) ---")
print(tabla_gsea)

# 7. GRÁFICO DE ENRIQUECIMIENTO (La montaña)
# Gráfico para DPN_perturbed
plotEnrichment(pathways_DPN[["DPN_perturbed"]], ranks_DPN) + 
  labs(title="GSEA: Firma DPN_perturbed (48h) a las 24h",
       subtitle="Tendencia de activación temprana",
       caption="Calculado mediante el estadístico de Wald (stat)") +
  theme_minimal()

