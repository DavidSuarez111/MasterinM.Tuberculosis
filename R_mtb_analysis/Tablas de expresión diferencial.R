#Expresiones diferenciales
#install.packages("openxlsx")
#install.packages("viridis")
library(openxlsx)
library(tidyverse)
library(conflicted)
library(viridis)
#by time 14

#BC391
BCvsBL_D14 <- dfExpresiónmasenBCqueBL #las sobrerepresentadas de BC

#BL323
BLvsBC_D14 <- dfExpresiónmasenBLqueBC #las sobrerepresentadas de BL

#by time 28
#BC391
BCvsBL_D28 <- dfExpresiónmasenBCqueBLD28 #sobrerepresentados en BCD28 vs BLD28
#BL323
BLvsBC_D28 <- dfExpresiónmasenBLqueBCD28 #sobrerepresentados en BLD28 vs BCD28

#14 vs 28 BC 
BC14vsBC28 <- dfExpresiónmasenBC14queBC28 #sobrerepresentadas en BC14 vs BC28
BC28vsBC14 <- dfExpresiónmasenBC28queBC14 #sobrerepresentadas en BC28 vs BC14
#14 vs 28 BL
BL14vsBL28 <- dfExpresiónmasenBL14queBL28 #sobrerepresentados en BL14 vs BL28
BL28vsBL14 <- dfExpresiónmasenBL28queBL14 #sobrerepresentados en BL28 vs BL14

#normalización puede jugar en contra aunque diga que una es mayor que la otra va a ser tomado como hay presencia en la muestra.
#Están sacadas de contexto.
workbookDEG <- createWorkbook()
addWorksheet(workbookDEG, "BCvsBL_D14")
writeData(workbookDEG, "BCvsBL_D14", BCvsBL_D14, rowNames = T)
addWorksheet(workbookDEG, "BLvsBC_D14")
writeData(workbookDEG, "BLvsBC_D14", BLvsBC_D14, rowNames = T)
addWorksheet(workbookDEG, "BCvsBL_D28")
writeData(workbookDEG, "BCvsBL_D28", BCvsBL_D28, rowNames = T)
addWorksheet(workbookDEG, "BLvsBC_D28")
writeData(workbookDEG, "BLvsBC_D28", BLvsBC_D28, rowNames = T)
addWorksheet(workbookDEG, "BC14vsBC28")
writeData(workbookDEG, "BC14vsBC28", BC14vsBC28, rowNames = T)
addWorksheet(workbookDEG, "BC28vsBC14")
writeData(workbookDEG, "BC28vsBC14", BC28vsBC14, rowNames = T)
addWorksheet(workbookDEG, "BL14vsBL28")
writeData(workbookDEG, "BL14vsBL28", BL14vsBL28, rowNames = T)
addWorksheet(workbookDEG, "BL28vsBL14")
writeData(workbookDEG, "BL28vsBL14", BL28vsBL14, rowNames = T)

saveWorkbook(workbookDEG, file = "BCvsBL_D14andD28.xlsx",  overwrite = TRUE)


#Datos ya con descripción DEG
#todo en referencia a comparación BL vs BC. 
#comparación al día 14
data_DEG_BCD14 <- read.xlsx("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/2. Transcriptomic Analysis Mastery in Immunology/BCvsBL_D14andD28_descrito.xlsx", sheet = 1)
data_DEG_BLD14 <- read.xlsx("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/2. Transcriptomic Analysis Mastery in Immunology/BCvsBL_D14andD28_descrito.xlsx", sheet = 2)
                    #colnames(data_DEG_BCD14)
#comparación al día 28
data_DEG_BCD28 <- read.xlsx("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/2. Transcriptomic Analysis Mastery in Immunology/BCvsBL_D14andD28_descrito.xlsx", sheet = 3)
data_DEG_BLD28 <- read.xlsx("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/2. Transcriptomic Analysis Mastery in Immunology/BCvsBL_D14andD28_descrito.xlsx", sheet = 4)
#comparación BCD14 vs BCD28
data_DEG_BCD14v28 <- read.xlsx("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/2. Transcriptomic Analysis Mastery in Immunology/BCvsBL_D14andD28_descrito.xlsx", sheet = 5)
data_DEG_BCD28v14 <- read.xlsx("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/2. Transcriptomic Analysis Mastery in Immunology/BCvsBL_D14andD28_descrito.xlsx", sheet = 6)
#comparación BLD14 vs BLD28
data_DEG_BLD14v28 <- read.xlsx("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/2. Transcriptomic Analysis Mastery in Immunology/BCvsBL_D14andD28_descrito.xlsx", sheet = 7)
data_DEG_BLD28v14 <- read.xlsx("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/2. Transcriptomic Analysis Mastery in Immunology/BCvsBL_D14andD28_descrito.xlsx", sheet = 8)
lista_data <- list(data_DEG_BCD14, data_DEG_BCD28, data_DEG_BLD14, data_DEG_BLD28, 
                   data_DEG_BCD14v28, data_DEG_BCD28v14, data_DEG_BLD14v28, data_DEG_BLD28v14)
str(lista_data)

funcion_columnas_deseadas <- function(df){
    df %>% 
    select(gen, baseMean, log2FoldChange, Nombre.de.gén) %>% 
    rename(base_mean = baseMean, log2fc = log2FoldChange, genename = Nombre.de.gén)
  }

bases_de_datos_filtradas <- sapply(lista_data, funcion_columnas_deseadas)
str(bases_de_datos_filtradas) #tiene 32 filas quiero dataframes de 4 columnas

separar_lista <- function(lista) {
  # Dividir la lista en 8 partes, cada una con 4 elementos
  lista_separada <- split(lista, ceiling(seq_along(lista) / 4))
  
  # Devolver la lista de listas separadas
  return(lista_separada)
}

bases_filtradas_en_lista <- separar_lista(bases_de_datos_filtradas)
str(bases_filtradas_en_lista)

view(bases_de_datos_filtradas)
BCD14 <- data.frame(bases_filtradas_en_lista[1])
BLD14 <- data.frame(bases_filtradas_en_lista[2])
BCD28 <- data.frame(bases_filtradas_en_lista[3])
BLD28 <- data.frame(bases_filtradas_en_lista[4])
BCD14vs28 <- data.frame(bases_filtradas_en_lista[5])
BCD28v14 <- data.frame(bases_filtradas_en_lista[6])
BLD14v28 <- data.frame(bases_filtradas_en_lista[7])
BLD28v14 <- data.frame(bases_filtradas_en_lista[8])
colnamesdf = c("gen", "base_mean", "log2fc", "geneproduct")
colnames(BCD14) <- colnamesdf
colnames(BLD14) <- colnamesdf
colnames(BCD28) <- colnamesdf
colnames(BLD28) <- colnamesdf
colnames(BCD14vs28) <- colnamesdf
colnames(BCD28v14) <- colnamesdf
colnames(BLD14v28) <- colnamesdf
colnames(BLD28v14) <- colnamesdf

lista_para_basesota <- list(BCD14, BLD14, BCD28, BLD28,  BCD14vs28, BCD28v14, BLD14v28, BLD28v14)

view(BCD14)

#hacer una basesota bien ordenada. #toca poner nuevas columnas por día, por cepa, por clasificación, por subclasificación
#objetivo gráficos de pie por clasificación o subclasificación



#Base de datos junta en proceso
#corregido_lista_basesota <- lapply(lista_para_basesota, function(df){
#  df %>% 
#    mutate(log2fc = as.numeric(log2fc))
#})
basesota <- bind_rows(lista_para_basesota)

#quizá anotar que viene de cual comparación para separarlo por separación ejemplo BCvsBL o BLvsBL y BCvsBC. 
#arreglar basesota. ponerle columna de cepa, columna de día, 
basesota$cepa = c(rep("BC", times = 20), rep("BL", times = 20), rep("BC", times = 20), rep("BL", times = 20),  
                  rep("BC", times = 40), rep("BL", times = 40))
basesota$dia = c(rep("D14", times = 40), rep("D28", times = 40), rep("D14", times = 20), rep("D28", times = 20), 
                  rep("D14", times = 20), rep("D28", times = 20))
basesota$ref_comp = c(rep("bcd14vsbld14", times = 40), rep("bcd28vsbld28", times = 40), rep("bcd14vsbcd28", times = 40), rep("bld14vsbld28", times = 40))

#ocupo conseguir información sobre clasificación y subclasificación supongo iniciaré con información no técnica
basesota$clasificacion = c() #160 fucking valores. si salieron. 
#basesota$subclasificacion = C() 

gene_product_unique <- unique(basesota$geneproduct)
print(gene_product_unique) #49 distintos productos en las 180 filas totales. 

view(basesota) 
str(basesota)

#BCD14, BLD14, BCD28, BLD28,  BCD14vs28, BCD28v14, BLD14v28, BLD28v14 #dfs ya filtrados

tipos_fv <- c("Proteínas Metabólicas", "Proteínas de Secreción", "Proteínas de Transcripción o Traducción", "Proteinas Adaptoras a Estrés", "Proteínas de Alteración Epigenética", "Proteínas Intercalantes de fagolisosoma", "Sideróforo", "Regulador Transcripcional", "RNA o procesamiento de RNA", "Lípidos", "Transposasa")
subtipo_fv <- c("Colesterol", "Grasas", "Glucosa","Nucleótidos", "Ag señuelo", "Proteína intercalante fagolisosomal", "RNA", "Alteración epigenética")


#formar columna vacía en df
basesota$clasificacion <- NA
View(basesota)
#basesota$geneproduct

#no me ha salido el for me iré a python así que haré un excel
basesotaDEG <- createWorkbook()
addWorksheet(basesotaDEG, "BaseDEG")
writeData(basesotaDEG, "BaseDEG", basesota, rowNames = T)
saveWorkbook(basesotaDEG, file = "BCvsBL_DEGppython.xlsx",  overwrite = TRUE)

##opción 2 unirlas directamente a una sola variable #pendiente a revisar; hice una lista con listas
#
#pendientes_a_clasificar <- basesota$geneproduct[is.na(basesota$clasificacion)]
#View(basesota)
#print(basesota)
#print(pendientes_a_clasificar)

#Importar excel 
basesota_clasificada <- read.xlsx("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/2. Transcriptomic Analysis Mastery in Immunology/genesDEGMtb.xlsx", sheet = 1)
view(basesota_clasificada)
colnames(basesota_clasificada) #Listo

#Resumen de la información para analizar

basesota_agrupada <- basesota_clasificada %>% 
  group_by(cepa, dia, ref_comp, clasificacion) %>%
  summarise(count = n(), .groups = "drop")
view(basesota_agrupada) #separar por ref comps para hacer los distintos pays. 

#separacion de basesota agrupada para poder comparar
#qué quiero comparar? 
#realmente las cepas entonces bc vs bl d14 y bcvsbld28, luego bl14 vs bl28 y bc14 vs bc28. 

basesota_agrupada_solo_dia14 = basesota_agrupada

#formar pie. #Revisar los colores, si ponerle lineas, conteos, considerar quitar datos. Revisar que si tengo los datos de entrada correctos porque parece replica el D14 y D28 de BL. 
ggplot(basesota_agrupada, aes(x = "", y = count, fill = clasificacion)) +
  geom_bar(stat = "identity", width = 1) +  # Create a stacked bar chart
  coord_polar("y") +  # Convert to pie chart
  facet_grid(cepa ~ dia ~ ref_comp) +  # Separate by strain and day
  labs(title = "Distribución de genes sobreexpresados por cepa y día",
       fill = "Grupos de proteínas") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))
#Se observa disminución de proteinas hipotéticas al dia28 en BC. Disminución de proteínas metabólicas sobreexpresadas. Menos proteínas de secreción, más proteínas de regulación transcripcional, proteínas de RNA o procesamiento y proteínas de transcripción y traducción. 
#Para contraste, retirar los hipotéticos. 
#En BL se observan las mismas tasas de grupos enzimáticos. Corroborando se encuentran los mismos tipos enzimáticos y los mismos genes. Lo siguiente es contrastar cantidades de cada una de estas enzimas. 
#
#Estos graficos de pay ya no los usaré porque los hice más organizados en python. 
#Opcion 2 estetizada se le pidió a chat gpt cambios estéticos
ggplot(basesota_agrupada, aes(x = "", y = count, fill = clasificacion)) +
  geom_bar(stat = "identity", width = 1)+#linetype = "solid") +  # Bordes blancos para separar segmentos
  coord_polar("y", start = 0) +  # Añade una rotación inicial para alinear el gráfico
  facet_grid(cepa ~ dia ~ ref_comp) +  # Separa por cepa y día
  labs(
    title = "Distribución de genes sobreexpresados por cepa y día",
    fill = "Grupos de proteínas"
  ) + # Cambiar paleta de colores a una más estética
  theme_minimal() +
  theme(
    axis.title = element_blank(),        # Quitar títulos de ejes
    axis.text = element_blank(),         # Quitar texto de ejes
    axis.ticks = element_blank(),        # Quitar marcas de ejes
    panel.grid = element_blank(),        # Quitar cuadrícula
    strip.text = element_text(size = 12, face = "bold"),  # Texto de facetas más claro y destacado
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centrar y resaltar título
    legend.position = "bottom",          # Mover leyenda a la parte inferior
    legend.title = element_text(size = 12, face = "bold"),  # Resaltar título de la leyenda
    legend.text = element_text(size = 10)  # Hacer texto de leyenda más legible
  )

#opcion 3 
ggplot(basesota_agrupada, aes(x = "", y = count, fill = clasificacion)) +
  geom_bar(stat = "identity", width = 1) +  # Create the pie chart
  coord_polar("y", start = 0) +  # Add rotation to align the chart
  facet_grid(cepa ~ dia) +  # Separate by cepa and day
  labs(
    title = "Distribución de genes sobreexpresados por cepa y día",
    fill = "Grupos de proteínas"
  ) +  # Add a more aesthetic color palette
  theme_minimal() +
  theme(
    axis.title = element_blank(),        # Remove axis titles
    axis.text = element_blank(),         # Remove axis text
    axis.ticks = element_blank(),        # Remove axis ticks
    panel.grid = element_blank(),        # Remove grid
    strip.text = element_text(size = 12, face = "bold"),  # Highlight facet text
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and highlight the title
    legend.position = "bottom",          # Move the legend to the bottom
    legend.title = element_text(size = 10, face = "bold"),  # Reduce the size of the legend title
    legend.text = element_text(size = 8),  # Reduce the size of legend text
    legend.key.width = unit(0.5, "inches"),  # Reduce the width of the legend keys
    legend.key.height = unit(0.5, "inches"),  # Make legend keys smaller in height
    plot.margin = margin(10, 10, 50, 10),  # Reduce the bottom margin for a larger pie chart area
    legend.margin = margin(t = 5, b = 5),  # Tighten the margin around the legend
    legend.box.spacing = unit(0.2, "cm"),  # Reduce spacing between legend items
    legend.box.margin = margin(t = 10)  # Add margin above the legend box
  ) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1))  # Use square legend keys

#Se está haciendo un análisis de los datos. Son los mismos tipos presentes en python. Otro seguimiento es asignar los valores de clasificación a clasificacion para pie llevando los grupos menores a 3 a otros ()

#Por el momento se harán gráficos de barras para comparar las cantidades de cada gen por tiempo en python y de ahí ya empezar a concluir y conectar las infos. 



