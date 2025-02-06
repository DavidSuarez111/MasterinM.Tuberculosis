##Este documento es para verificar calidad y controlar calidad de los archivos Fastq


#Espacio de trabajo 
getwd()
#setwd(pathBCd14)
#pathGenomasref <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Genomas referencia"

pathBCd14 <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/BC d14"
pathBLd14 <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/BL d14"

#Cargar paquetes y librerías #se intentó con fastqcr pero es solo para mac o linux
#install.packages("ShortRead")
#install.packages("magrittr")
#install.packages("scales")

#correr estas
library(ShortRead)
library(stringr)
library(magrittr)
library(tidyverse)
#library(dplyr)
library(scales)


# Revisión de calidad de archivos fastq de BC391 y BL323 ------------------


##revisión de calidad fastq de BC391 d14. Pasos? leer secuencias, ...?
#Acceder a los datos
filesBCd14 <- list.files(pathBCd14)
print(filesBCd14)
str(filesBCd14)
filesBLd14 <- list.files(pathBLd14)
print(filesBLd14)

##Análisis por Shortread por archivo
#Directorio de rutas a los archivos para BC391d14 y BL323d14
BCD14_1_12_dirs <- character()
BCD14_1_12_dirs <- paste0(pathBCd14, "/", filesBCd14)
BCD14_1_12_dirs[1]
BCD14_1_12_dirs

BLD14_1_12_dirs <- character()
BLD14_1_12_dirs <- paste0(pathBLd14, "/", filesBLd14)

#Conteos de Lecturas por archivo (mediante countFastq)
##countFastq
#de los countFastq se puede observar que los archivos están compuestos de 278K a 822K de lecturas y de 20.5 millones a 60.7 millones de nucleotidos por archivo. 
 

 countarch1BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[1]) #para la función countFastq se utiliza la dirección de los archivos con la terminación .fastq ergo las variables con solo un guión bajo ej. dirBCD14_1
 countarch2BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[2])
 countarch3BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[3])
 countarch4BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[4])
 countarch5BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[5])
 countarch6BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[6])
 countarch7BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[7])
 countarch8BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[8])
 countarch9BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[9])
countarch10BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[10])
countarch11BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[11])
countarch12BCD14 <- countFastq(dirPath = BCD14_1_12_dirs[12])


#Desarrollo de elegancia
#countarch_1a12_BCD14 <- sapply(BCD14_1_12_dirs, countFastq(dirPath = BCD14_1_12_dirs))

#Organizar los datos de Lecturas por archivo en una tabla
countarchBCD14list <- list(countarch1BCD14, countarch2BCD14, countarch3BCD14, countarch4BCD14, countarch5BCD14, countarch6BCD14, countarch7BCD14, countarch8BCD14, countarch9BCD14, countarch10BCD14, countarch11BCD14, countarch12BCD14)
print(countarchBCD14list)
view(countarchBCD14list) 

tibblecountBCD14 <- countarchBCD14list %>%
  map_dfr(~.x)
view(tibblecountBCD14)#organizado por hileras y columnas
rnBC <- row.names(tibblecountBCD14)

#Organizar la información y cambiar los nombres de las columnas
#TitleBCD14 <- "Datos de los archivos SRA de BC391 al día 14"
#NotadepieBCD14 <- "Fuente: Datos obtenidos de bioprojecto con número de acceso: "
dfcountBCD14_E1 <- tibblecountBCD14 %>% 
  rename(
  Lecturas = "records", 
  Nucleotidos = "nucleotides",
  Puntuaciones = "scores"
) %>% 
  select(- Puntuaciones) %>%
  mutate(Archivos = rnBC) %>%
  relocate(Archivos, .before = everything())
str(dfcountBCD14_E1)
typeof(dfcountBCD14_E1)
view(dfcountBCD14_E1) #Resultado: Tabla bien organizado de la composición de lecturas por archivo de la cepa BC391 al día 14 en la variable dfcountBCD14_E1

#Replica con Datos BL
countarch1BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[1])
countarch2BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[2])
countarch3BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[3])
countarch4BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[4])
countarch5BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[5])
countarch6BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[6])
countarch7BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[7])
countarch8BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[8])
countarch9BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[9])
countarch10BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[10])
countarch11BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[11])
countarch12BLD14 <- countFastq(dirPath = BLD14_1_12_dirs[12])

countarchBLD14list <- list(countarch1BLD14, countarch2BLD14, countarch3BLD14, countarch4BLD14, countarch5BLD14, countarch6BLD14, countarch7BLD14, countarch8BLD14, countarch9BLD14, countarch10BLD14, countarch11BLD14, countarch12BLD14)
print(countarchBLD14list)

tibblecountBLD14 <- countarchBLD14list %>%
  map_dfr(~.x) 
view(tibblecountBLD14) 
rnBL <-row.names(tibblecountBLD14)
str(tibblecountBLD14)

#TitleBLD14 <- "Datos de los archivos SRA de BL323 al día 14"
#NotadepieBLD14 <- "Fuente: Datos obtenidos de bioprojecto con número de acceso: "
dfcountBLD14_E1 <- tibblecountBLD14 %>% 
  rename(
    Lecturas = "records", 
    Nucleotidos = "nucleotides",
    Puntuaciones = "scores"
  ) %>% 
  select(- Puntuaciones) %>%  #Resultado de tabla organizada de lecturas por archivo de tipo beijing 323. Variable: dfcountBLD14_E1
  mutate(Archivos = rnBL) %>% 
  relocate(Archivos, .before = everything())
  
view(dfcountBLD14_E1) 
str(dfcountBLD14_E1)
#Resultados countfastq: Nombres de los archivos y número de lecturas por archivo. 
#BC391D14:
dfcountBCD14_E1
#BL323D14
dfcountBLD14_E1

#los he de unir en una tabla? #generar gráficos (poder de secuenciación, comparativa entre secuenciados)

groupBR <- rep(c("BR1", "BR2", "BR3"), each = 4)

dfBCg1 <- dfcountBCD14_E1 %>%
  select(-Nucleotidos) %>%
  mutate(Grupo = groupBR) %>% 
  transmute(Archivos = c("BR1T1", "BR1T2", "BR1T3", "BR1T4", "BR2T1", "BR2T2", "BR2T3", "BR2T4","BR3T1", "BR3T2", "BR3T3", "BR3T4"), 
            Lecturas, 
            Grupo)#Añadir separación por replicados biológicos para añadir color.
view(dfBCg1)
dfBCg2 <- dfcountBCD14_E1 %>% 
  select(-Lecturas) %>%
  mutate(Grupo = groupBR) %>% 
  transmute(Archivos = c("BR1T1", "BR1T2", "BR1T3", "BR1T4", "BR2T1", "BR2T2", "BR2T3", "BR2T4","BR3T1", "BR3T2", "BR3T3", "BR3T4"), 
            Nucleotidos, 
            Grupo) 
summary(dfBCg1) #20MB a 60.7 MB
summary(dfBCg2) #278K a 822K lecturas. 


ggplot(dfBCg1,aes(x= Archivos, y = Lecturas/1000, fill = Grupo)) +
  geom_col() +
   theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.x = element_blank()) +
  labs(y ="Lecturas (miles)", x ="Archivos (BC391)", title = "Número de Lecturas entre los Archivos de BC391 al Día 14") + #subtitle = "Número de Lecturas por Archivo") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ 
  scale_y_continuous(breaks = seq(0, 1000, by = 100), label= seq(0, 1000, by = 100)) +
  coord_cartesian(ylim = c(0, 1000), xlim = c(1,12))

ggplot(dfBCg2,aes(x= Archivos, y = Nucleotidos/1000000, fill = Grupo)) +
  geom_col() +
  labs(y ="Nucleótidos (Mb)", x ="Archivos (BC391)", title ="Megabases por Archivo de BC391 al Día 14") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
  #scale_y_continuous(breaks = seq(0, 1000, by = 100), label= seq(0, 1000, by = 100)) +
  #coord_cartesian(ylim = c(0, 1000), xlim = c(1,12))


#Replicar graficas para BL

groupBR <- rep(c("BR1", "BR2", "BR3"), each = 4)
dfBLg1 <- dfcountBLD14_E1 %>%
  select(-Nucleotidos) %>%
  mutate(Grupo = groupBR) %>% 
  transmute(Archivos = c("BR1T1", "BR1T2", "BR1T3", "BR1T4", "BR2T1", "BR2T2", "BR2T3", "BR2T4","BR3T1", "BR3T2", "BR3T3", "BR3T4"), 
            Lecturas, 
            Grupo)#Añadir separación por replicados biológicos para añadir color.
view(dfBCg1)
dfBLg2 <- dfcountBLD14_E1 %>% 
  select(-Lecturas) %>%
  mutate(Grupo = groupBR) %>% 
  transmute(Archivos = c("BR1T1", "BR1T2", "BR1T3", "BR1T4", "BR2T1", "BR2T2", "BR2T3", "BR2T4","BR3T1", "BR3T2", "BR3T3", "BR3T4"), 
            Nucleotidos, 
            Grupo) 
summary(dfBLg1)
summary(dfBLg2)


ggplot(dfBLg1,aes(x= Archivos, y = Lecturas/1000, fill = Grupo)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.x = element_blank()) +
  labs(y ="Lecturas (miles)", x ="Archivos (BL323)", title = "Número de Lecturas entre Archivos BL323 al Día 14") +#, subtitle = "Número de Lecturas por Archivo") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ 
  scale_y_continuous(breaks = seq(0, 1000, by = 100), label= seq(0, 1000, by = 100)) +
  coord_cartesian(ylim = c(0, 1000), xlim = c(1,12))

ggplot(dfBLg2,aes(x= Archivos, y = Nucleotidos/1000000, fill = Grupo)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y ="Nucleótidos (Mb)", x ="Archivos (BL323)", title = "Megabases por Archivo BL323 al Día 14") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
##Continuación. Abrir los archivos en formato ShortReadQ.  
#Para BCd14
arch1BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[1])
arch2BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[2])
arch3BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[3])
arch4BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[4])
arch5BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[5])
arch6BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[6])
arch7BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[7])
arch8BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[8])
arch9BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[9])
arch10BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[10])
arch11BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[11])
arch12BCD14 <- readFastq(dirPath = BCD14_1_12_dirs[12])

ListarchBCD14 <- list(arch1BCD14,  arch2BCD14, arch3BCD14, arch4BCD14,
                      arch5BCD14,  arch6BCD14, arch7BCD14, arch8BCD14, 
                      arch9BCD14, arch10BCD14,arch11BCD14,arch12BCD14)
ListarchBCD14[1]

#Para BLd14
 arch1BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[1])
 arch2BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[2])
 arch3BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[3])
 arch4BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[4])
 arch5BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[5])
 arch6BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[6])
 arch7BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[7])
 arch8BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[8])
 arch9BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[9])
arch10BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[10])
arch11BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[11])
arch12BLD14 <- readFastq(dirPath = BLD14_1_12_dirs[12])

ListarchBLD14 <- list(arch1BLD14, arch2BLD14, arch3BLD14, arch4BLD14,
                      arch5BLD14, arch6BLD14, arch7BLD14, arch8BLD14, 
                      arch9BLD14,arch10BLD14,arch11BLD14,arch12BLD14)


methods(class = "ShortRead") #funciones para utilizar en archivos de tipo ShortReadQ

#Revisión de calidad y edición para garantizar calidad. 
##A partir de aquí se hará una exploración estadística para verificar la calidad de los archivos y  
##definir los parámetros se harán los recortes de los archivos fastq/ShortReadq
#Tentativamente son 3 parámetros: longitud atípica de las lecturas (<74, tentativamente puede que esta no), eliminación de tails primeros 5-6 nt (baja calidad); Calidad >30 = 99.9% de veracidad de la secuencia = buena calidad (este último revisar si tiene pertinencia). 

#De Dr. Pato trimming de extremos basado en calidad. Calidad global del read por encima de un umbral. Quizá por longitud <30. 
arch1BCD14 #comenta la diversidad de longitudes de lecturas (50-76) y cantidad de lecturas

w1 <- width(arch1BCD14)
L1 <- length(arch1BCD14)
dfw1 <- as.data.frame(table(w1))
view(dfw1)
dfw1ysuma <- dfw1 %>% 
  mutate(Porcentaje = signif(Freq/sum(Freq)*100, digits = 2),
         SumaPorcientos = round(cumsum(Porcentaje), digits = 2)) %>%
  arrange(desc(w1))
view(dfw1ysuma) #si nos quedamos con las longitudes de lecturas de 76 a 73 nos quedamos con el 86% de los datos. 

df_hist_width <- dfw1ysuma %>% 
  select(w1, Freq) %>% 
  mutate(grupo = c("highlight", "highlight", "highlight", rep("normal", 24)), 
         Porcentaje = signif((Freq/sum(Freq))*100, digits = 2), 
         Porcentajeg = ifelse(Porcentaje >= 5, paste(Porcentaje, "%"),  "")) #%>%
  #transmute(w1, Freq, grupo, Porcentaje => 5)#filter(Freq >0)
view(df_hist_width)

ggplot(df_hist_width,aes(x= w1, y = Freq/1000, fill = grupo)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.x = element_blank(), legend.position = "none") +
  labs(y ="Conteo (miles)", x ="Longitud de Lecturas (pb)", title = "Longitud de Lecturas en Archivo 1 de BC391 al Día 14") +#, subtitle = "") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
scale_fill_manual(values = c("normal" = "gray", "highlight" = "red")) +
  geom_text(aes(label = paste(Porcentajeg)), vjust = -0.3, hjust= 0.6, size = 6)
#darle color, quizá para comunicar la contribución ejemplo si es menos del 1% en 1 color, si es más en otros o agregar números de porcentajes a los valores de y con mayor cobertura (74-76)

# Longitud de Lecturas BC391 ----------------------------------------------
fun_longitud_lecturasBC <- function(x){
  x1 <- width(x)
  x2 <- as.data.frame(table(x1))
  x3 <- x2 %>% 
    mutate(Porcentaje = signif(Freq/sum(Freq)*100, digits = 2),
           SumaPorcientos = round(cumsum(Porcentaje), digits = 2),
           Porcentajeg = ifelse(Porcentaje >= 5, paste(Porcentaje, "%"),  ""),
           Color = ifelse(Porcentaje > 5, "Highlight", "Normal")) %>%
    arrange(desc(x1))
}


#make an sapply to get all the widths; worked
#for(L in ListarchBCD14){
  RepsporLecsBC <- sapply(ListarchBCD14, FUN = fun_longitud_lecturasBC)

RepsporLecsBC 
  

#asignar a variables y unir todas a un solo data frame. Añadir valor de archivo proveniente y llevar a facet wrap. 
df_fw_BC <- bind_rows(RepsporLecsBC[,1], RepsporLecsBC[,2], RepsporLecsBC[,3], RepsporLecsBC[,4],
                      RepsporLecsBC[,5], RepsporLecsBC[,6], RepsporLecsBC[,7], RepsporLecsBC[,8],
                      RepsporLecsBC[,9],RepsporLecsBC[,10],RepsporLecsBC[,11],RepsporLecsBC[,12])
view(df_fw_BC)
  
df_BC_conid <- df_fw_BC %>% 
  mutate(Archivo = c(rep("BR1T1", 27), rep("BR1T2", 27), rep("BR1T3", 27), rep("BR1T4", 27),
                     rep("BR2T1", 27), rep("BR2T2", 27), rep("BR2T3", 27), rep("BR2T4", 27),
                     rep("BR3T1", 27), rep("BR3T2", 27), rep("BR3T3", 27), rep("BR3T4", 27)))
view(df_BC_conid)

#visualizar facet wrap

  df_hist_width_BC <- df_BC_conid %>% 
    select(x1, Freq, Archivo, Porcentajeg, Color) 
  view(df_hist_width_BC)
  
  ggplot(df_hist_width_BC,aes(x= x1, y = Freq/1000, fill = Color)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.x = element_blank(), legend.position = "none") +
    labs(y ="Conteo (miles)", x ="Longitud de Lecturas (pb)", title = "Longitud de Lecturas en Archivos BC391") +#, subtitle = "") +
    theme(plot.title = element_text(size = 30),
          plot.subtitle = element_text(size = 20),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25, angle = 90),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18))+
    facet_wrap(~Archivo) +
    scale_x_discrete(breaks = seq(50, 80, by = 5), label= seq(50, 80, by = 5)) +
  geom_text(aes(label = paste(Porcentajeg)), vjust = -0.3, hjust= 0.8, size = 4) +
    scale_fill_manual(values = c("Normal" = "gray", "Highlight" = "red"))
  #darle color, quizá para comunicar la contribución ejemplo si es menos del 1% en 1 color, si es más en otros o agregar números de porcentajes a los valores de y con mayor cobertura (74-76)

# Longitud de Lecturas BL323 ----------------------------------------------
  fun_longitud_lecturasBL <- function(x){
    x1 <- width(x)
    x2 <- as.data.frame(table(x1))
    x3 <- x2 %>% 
      mutate(Porcentaje = signif(Freq/sum(Freq)*100, digits = 2),
             SumaPorcientos = round(cumsum(Porcentaje), digits = 2),
             Porcentajeg = ifelse(Porcentaje >= 5, paste(Porcentaje, "%"),  ""),
             Color = ifelse(Porcentaje > 5, "Highlight", "Normal")) %>%
      arrange(desc(x1))
  }

RepsporLecsBL <- sapply(ListarchBLD14, FUN = fun_longitud_lecturasBL)
RepsporLecsBL 

#Replicar para BL323
  df_fw_BL <- bind_rows(RepsporLecsBL[,1], RepsporLecsBL[,2], RepsporLecsBL[,3], RepsporLecsBL[,4],
                        RepsporLecsBL[,5], RepsporLecsBL[,6], RepsporLecsBL[,7], RepsporLecsBL[,8],
                        RepsporLecsBL[,9], RepsporLecsBL[,10],RepsporLecsBL[,11],RepsporLecsBL[,12])
  df_fw_BL
  df_BL_conid <- df_fw_BL %>% 
    mutate(Archivo = c(rep("BR1T1", 27), rep("BR1T2", 27), rep("BR1T3", 27), rep("BR1T4", 27),
                       rep("BR2T1", 27), rep("BR2T2", 27), rep("BR2T3", 27), rep("BR2T4", 27),
                       rep("BR3T1", 27), rep("BR3T2", 27), rep("BR3T3", 27), rep("BR3T4", 27)))
  
  df_hist_width_BL <- df_BL_conid %>% 
    select(x1, Freq, Archivo, Porcentajeg, Color)
  view(df_hist_width_BL)
#Visualizar con Facet wrap
  #longitud de lecturas
  ggplot(df_hist_width_BL,aes(x= x1, y = Freq/1000, fill = Color)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.x = element_blank(), legend.position = "none") +
    labs(y ="Conteo (miles)", x ="Longitud de Lecturas (pb)", title = "Longitud de Lecturas en Archivos BL323")+ #, subtitle = "") +
    theme(plot.title = element_text(size = 30),
          plot.subtitle = element_text(size = 25),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 20, angle = 90),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18))+
    facet_wrap(~Archivo) +
    scale_x_discrete(breaks = seq(50, 80, by = 5), label= seq(50, 80, by = 5)) +
    geom_text(aes(label = paste(Porcentajeg)), vjust = -0.3, hjust= 0.8, size = 4) +
    scale_fill_manual(values = c("Normal" = "gray", "Highlight" = "red"))
  
#para BCD14

#para BLD14

#2do filtro calidades y tails
#del qa summary.
  
#Todos los archivos


qa_BC_BR1T1 <-  qa(arch1BCD14, lane = 1)
qa_BC_BR1T2 <-  qa(arch2BCD14, lane = 1)
qa_BC_BR1T3 <-  qa(arch3BCD14, lane = 1)
qa_BC_BR1T4 <-  qa(arch4BCD14, lane = 1)
qa_BC_BR2T1 <-  qa(arch5BCD14, lane = 1)
qa_BC_BR2T2 <-  qa(arch6BCD14, lane = 1)
qa_BC_BR2T3 <-  qa(arch7BCD14, lane = 1)
qa_BC_BR2T4 <-  qa(arch8BCD14, lane = 1)
qa_BC_BR3T1 <-  qa(arch9BCD14, lane = 1)
qa_BC_BR3T2 <- qa(arch10BCD14, lane = 1)
qa_BC_BR3T3 <- qa(arch11BCD14, lane = 1)
qa_BC_BR3T4 <- qa(arch12BCD14, lane = 1)

calidad_BC <- c(qa_BC_BR1T1, qa_BC_BR1T2, qa_BC_BR1T3, qa_BC_BR1T4,
                qa_BC_BR2T1, qa_BC_BR2T2, qa_BC_BR2T3, qa_BC_BR2T4,
                qa_BC_BR3T1, qa_BC_BR3T2, qa_BC_BR3T3, qa_BC_BR3T4)
calidad_BC

qa_BL_BR1T1 <-  qa(arch1BLD14, lane = 1)
qa_BL_BR1T1[["readQualityScore"]]
qa_BL_BR1T2 <-  qa(arch2BLD14, lane = 1)
qa_BL_BR1T3 <-  qa(arch3BLD14, lane = 1)
qa_BL_BR1T4 <-  qa(arch4BLD14, lane = 1)
qa_BL_BR2T1 <-  qa(arch5BLD14, lane = 1)
qa_BL_BR2T2 <-  qa(arch6BLD14, lane = 1)
qa_BL_BR2T3 <-  qa(arch7BLD14, lane = 1)
qa_BL_BR2T4 <-  qa(arch8BLD14, lane = 1)
qa_BL_BR3T1 <-  qa(arch9BLD14, lane = 1)
qa_BL_BR3T2 <- qa(arch10BLD14, lane = 1)
qa_BL_BR3T3 <- qa(arch11BLD14, lane = 1)
qa_BL_BR3T4 <- qa(arch12BLD14, lane = 1)

calidad_BL <- c(qa_BL_BR1T1, qa_BL_BR1T2, qa_BL_BR1T3, qa_BL_BR1T4,
                qa_BL_BR2T1, qa_BL_BR2T2, qa_BL_BR2T3, qa_BL_BR2T4,
                qa_BL_BR3T1, qa_BL_BR3T2, qa_BL_BR3T3, qa_BL_BR3T4)

str(calidad_BL) #Tipo lista
calidad_BL


# Calidad de Phred en BC391 y BL323 ---------------------------------------

función_BC_cal <- function(x){
  x1 <- x[["baseQuality"]] %>%
    select(score, count) %>% 
    mutate(Puntaje = -1:93)
  x2 <- x1[-1,] %>%
    mutate(Porcentaje = signif(((count/sum(count))*100), digits = 2)) %>% 
    select(-score) %>% 
    filter(count > 28)
}

función_BL_cal <- function(x){
  x1 <- x[["baseQuality"]] %>%
    select(score, count) %>% 
    mutate(Puntaje = -1:93)
  x2 <- x1[-1,] %>%
    mutate(Porcentaje = signif(((count/sum(count))*100), digits = 2)) %>% 
    select(-score) %>% 
    filter(count > 28)
}

#listofdfs_BC_qual
#listofdfs_BL_qual

#aplicar la función 
listofdfs_BC_qual <- sapply(calidad_BC, función_BC_cal)
listofdfs_BL_qual <- sapply(calidad_BL, función_BL_cal)

listofdfs_BC_qual 
listofdfs_BL_qual 

listofdfs_BC_qual_g1<- bind_rows(listofdfs_BC_qual[,1], listofdfs_BC_qual[,2], listofdfs_BC_qual[,3], listofdfs_BC_qual[,4],
                                 listofdfs_BC_qual[,5], listofdfs_BC_qual[,6], listofdfs_BC_qual[,7], listofdfs_BC_qual[,8],
                                 listofdfs_BC_qual[,9], listofdfs_BC_qual[,10],listofdfs_BC_qual[,11],listofdfs_BC_qual[,12])


listofdfs_BL_qual_g1<- bind_rows(listofdfs_BL_qual[,1], listofdfs_BL_qual[,2], listofdfs_BL_qual[,3], listofdfs_BL_qual[,4],
                                 listofdfs_BL_qual[,5], listofdfs_BL_qual[,6], listofdfs_BL_qual[,7], listofdfs_BL_qual[,8],
                                 listofdfs_BL_qual[,9],listofdfs_BL_qual[,10],listofdfs_BL_qual[,11],listofdfs_BL_qual[,12])

listofdfs_BC_qual_g1_fnot0 <- listofdfs_BC_qual_g1 %>% 
                              filter(Puntaje > 13) %>% 
                              mutate(Archivo = c(rep("BR1T1", 5), rep("BR1T2", 5), rep("BR1T3", 5), rep("BR1T4", 5),
                                                 rep("BR2T1", 5), rep("BR2T2", 5), rep("BR2T3", 5), rep("BR2T4", 5),
                                                 rep("BR3T1", 5), rep("BR3T2", 5), rep("BR3T3", 5), rep("BR3T4", 5)), 
                                     Conteo = count,
                                     LosPorcentajesChidos = ifelse(Porcentaje >= 12, paste(Porcentaje, "%"), ""),
                                     Color = ifelse(LosPorcentajesChidos >= 12, "Highlight", "Normal")) %>% #agregar id para facet wrap  
                              select(-count) %>% 
                              relocate(Archivo, Puntaje, Conteo, Porcentaje) 
listofdfs_BL_qual_g1_fnot0 <- listofdfs_BL_qual_g1 %>% 
                              filter(Puntaje > 13) %>%
                              mutate(Archivo = c(rep("BR1T1", 5), rep("BR1T2", 5), rep("BR1T3", 5), rep("BR1T4", 5),
                                                 rep("BR2T1", 5), rep("BR2T2", 5), rep("BR2T3", 5), rep("BR2T4", 5),
                                                 rep("BR3T1", 5), rep("BR3T2", 5), rep("BR3T3", 5), rep("BR3T4", 5)),
                                     Conteo = count,
                                     LosPorcentajesChidos = ifelse(Porcentaje >= 12, paste(Porcentaje, "%"),  ""),
                                     Color = ifelse(LosPorcentajesChidos >= 12, "Highlight", "Normal")) %>% 
                              select(-count) %>%
                              relocate(Archivo, Puntaje, Conteo, Porcentaje)  
view(listofdfs_BC_qual_g1_fnot0)
view(listofdfs_BL_qual_g1_fnot0)

#Visualizar calidades generales
ggplot(listofdfs_BC_qual_g1_fnot0, aes(x = Puntaje , y= Conteo/1000000,)) + #fill = group)) +
  geom_col() +
  labs(y ="Conteo (millones)", x ="Calidad (PQ)", title = "Conteo de Puntajes de Calidad de Bases Nucleotídicas BC391", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = c(14,21, 27, 32, 36), label= c(14, 21, 27, 32, 36)) +
  theme(axis.text.y = element_text(angle = 0), legend.position = "none") +
  scale_y_continuous(labels = comma, limits = c(0, 50, by = 5)) +
  geom_text(aes(label = paste(LosPorcentajesChidos)), vjust = -0.3, size = 3) +
  facet_wrap(~Archivo)#scale_fill_manual(values = c("normal" = "gray", "highlight" = "red")) 

ggplot(listofdfs_BL_qual_g1_fnot0, aes(x = Puntaje , y= Conteo/1000000,)) + #fill = group)) +
  geom_col() +
  labs(y ="Conteo (millones)", x ="Calidad (PQ)", title = "Conteo de Puntajes de Calidad de Bases Nucleotídicas BL323", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = c(14,21, 27, 32, 36), label= c(14, 21, 27, 32, 36)) +
  theme(axis.text.y = element_text(angle = 0), legend.position = "none") +
  scale_y_continuous(labels = comma, limits = c(0, 50, by = 5)) +
  geom_text(aes(label = paste(listofdfs_BC_qual_g1_fnot0$LosPorcentajesChidos)), vjust = -0.4, hjust = 0.3, size = 3) +
  facet_wrap(~Archivo)#scale_fill_manual(values = c("normal" = "gray", "highlight" = "red")) 


# Calidad por Ciclo BC391 -------------------------------------------------


#Por ciclo todos los datos
calidad_BC

fun_BC_percycle <- function(x){
  x1 <- x[["perCycle"]]
  x2 <- x1$quality %>% 
    select(-Quality, -lane) %>% 
    mutate(Puntajeacolor =  ifelse(Score == 14, "Mala Calidad",
                                   ifelse(Score == 21, "Baja Calidad",
                                          ifelse(Score == 27, "Calidad Neutra",
                                                 ifelse(Score == 32, "Calidad Aceptable", "Calidad Excelente")))))
}

listdfs_BC_percycle_g1 <- sapply(calidad_BC, fun_BC_percycle)
listdfs_BC_percycle_g1
#387, 387, 388, 384, 384, 381, 394, 392, 382, 380, 395, 392


BC_percycle <- bind_rows(listdfs_BC_percycle_g1[,1], listdfs_BC_percycle_g1[,2], listdfs_BC_percycle_g1[,3],
                         listdfs_BC_percycle_g1[,4], listdfs_BC_percycle_g1[,5], listdfs_BC_percycle_g1[,6],
                         listdfs_BC_percycle_g1[,7], listdfs_BC_percycle_g1[,8], listdfs_BC_percycle_g1[,9],
                         listdfs_BC_percycle_g1[,10],listdfs_BC_percycle_g1[,11],listdfs_BC_percycle_g1[,12])

BC_percycle_g1 <- BC_percycle %>%
                  mutate(Archivo = c(rep("BR1T1", 387), rep("BR1T2", 387), rep("BR1T3", 388), rep("BR1T4", 384),
                                                        rep("BR2T1", 384), rep("BR2T2", 381), rep("BR2T3", 394), rep("BR2T4", 392),
                                                        rep("BR3T1", 382), rep("BR3T2", 380), rep("BR3T3", 395), rep("BR3T4", 392)))
BC_percycle_g1$Puntajeacolor <- factor(BC_percycle_g1$Puntajeacolor, levels = c("Mala Calidad", "Baja Calidad", "Calidad Neutra", "Calidad Aceptable", "Calidad Excelente"))
view(BC_percycle_g1)

ggplot(BC_percycle_g1, aes(x = Cycle, y = Count/1000, fill = Puntajeacolor)) +
  geom_col(position = "stack") + 
  scale_fill_manual(values = c("Mala Calidad" = "darkred", "Baja Calidad" = "red", "Calidad Neutra" = "lightskyblue", "Calidad Aceptable" = "lightgreen", "Calidad Excelente"= "green"))+
  labs(y ="Conteo (miles)", x ="Ciclo (Posición Nt)", title = "Conteo de Calidad de nucleótidos por Ciclo en BC391", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = seq(0, 75, by = 15), label= seq(0, 75, by = 15)) +
  scale_y_continuous(breaks = seq(0, 700, by = 100), label= seq(0, 700, by = 100))+
  facet_wrap(~Archivo)
#cambiar eje x, liberar y, ordenar por factor de puntaje a color

#Calidad por lecturas BC
función_qualxLec_BC <- function(x){
  x1 <- x[["readQualityScore"]]
  x2 <- x1 %>% select(-lane, -type)
}

qual_Lec_BC <- sapply(calidad_BC, función_qualxLec_BC)
qual_Lec_BC[,1]
qual_LEC_BC_df <- bind_rows(qual_Lec_BC[,1], qual_Lec_BC[,2], qual_Lec_BC[,3], qual_Lec_BC[,4],
                            qual_Lec_BC[,5], qual_Lec_BC[,6], qual_Lec_BC[,7], qual_Lec_BC[,8],
                            qual_Lec_BC[,9], qual_Lec_BC[,10], qual_Lec_BC[,11], qual_Lec_BC[,12])
qual_LEC_BC_df

qual_LEC_BC_df_1 <- qual_LEC_BC_df %>% 
  mutate(Tag = c(rep("BR1T1", 512), rep("BR1T2", 512), rep("BR1T3", 512), rep("BR1T4", 512),  
                 rep("BR2T1", 512), rep("BR2T2", 512), rep("BR2T3", 512), rep("BR2T4", 512),
                 rep("BR3T1", 512), rep("BR3T2", 512), rep("BR3T3", 512), rep("BR3T4", 512)),
         Den_menor30 = ifelse(quality >= 30, 0, density))
Área_rechazoBC <- signif(sum(qual_LEC_BC_df_1$Den_menor30)/sum(qual_LEC_BC_df_1$density)*100, 3)
view(qual_LEC_BC_df_1)
sumofdensity <- sum(qual_LEC_BC_df$density) #falta poner un porcentaje de cuanto se retirará

ggplot(qual_LEC_BC_df_1, aes(quality)) +geom_density(fill = "lightblue") +
  geom_vline(xintercept = 30) + 
  labs(y ="Proporción (%)", x ="Calidad (PQ)", title = "Calidad de Lecturas en BC391", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  geom_text(x=29.5, y=0.04, label = paste(Área_rechazoBC, "%"), size = 3.5) 





# Calidad por ciclo BL ----------------------------------------------------

calidad_BL
fun_BL_percycle <- function(x){
  x1 <- x[["perCycle"]]
  x2 <- x1$quality %>% 
    select(-Quality, -lane) %>% 
    mutate(Puntajeacolor =  ifelse(Score == 14, "Mala Calidad",
                                   ifelse(Score == 21, "Baja Calidad",
                                          ifelse(Score == 27, "Calidad Neutra",
                                                 ifelse(Score == 32, "Calidad Aceptable", "Calidad Excelente")))))
}

listdfs_BL_percycle_g1 <- sapply(calidad_BL, fun_BL_percycle)

listdfs_BL_percycle_g1
#384, 388, 382, 385, 386, 382, 390, 392, 390, 382, 391, 391

BL_percycle <- bind_rows(listdfs_BL_percycle_g1[,1], listdfs_BL_percycle_g1[,2], listdfs_BL_percycle_g1[,3], 
                         listdfs_BL_percycle_g1[,4], listdfs_BL_percycle_g1[,5], listdfs_BL_percycle_g1[,6],
                         listdfs_BL_percycle_g1[,7], listdfs_BL_percycle_g1[,8], listdfs_BL_percycle_g1[,9],
                         listdfs_BL_percycle_g1[,10],listdfs_BL_percycle_g1[,11],listdfs_BL_percycle_g1[,12])


BL_percycle_g1 <- BL_percycle %>%
  mutate(Archivo = c(rep("BR1T1", 384), rep("BR1T2", 388), rep("BR1T3", 382), rep("BR1T4", 385),
                     rep("BR2T1", 386), rep("BR2T2", 382), rep("BR2T3", 390), rep("BR2T4", 392),
                     rep("BR3T1", 390), rep("BR3T2", 382), rep("BR3T3", 391), rep("BR3T4", 391))) 
BL_percycle_g1$Puntajeacolor <- factor(BL_percycle_g1$Puntajeacolor, levels = c("Mala Calidad", "Baja Calidad", "Calidad Neutra", "Calidad Aceptable", "Calidad Excelente"))


BL_percycle_g1

ggplot(BL_percycle_g1, aes(x = Cycle, y = Count/1000, fill = Puntajeacolor)) +
  geom_col(position = "stack") + 
  scale_fill_manual(values = c("Mala Calidad" = "darkred", "Baja Calidad" = "red", "Calidad Neutra" = "lightskyblue", "Calidad Aceptable" = "lightgreen", "Calidad Excelente"= "green"))+
  labs(y ="Conteo (miles)", x ="Ciclo (Posición Nt)", title = "Conteo de Calidad de nucleótidos por Ciclo en BL323", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = seq(0, 75, by = 15), label= seq(0, 75, by = 15)) +
  scale_y_continuous(breaks = seq(0, 500, by = 100), label= seq(0, 500, by = 100))+
  facet_wrap(~Archivo)

#Calidad por lecturas BL
función_qualxLec_BL <- function(x){
  x1 <- x[["readQualityScore"]]
  x2 <- x1 %>% select(-lane, -type)
}

qual_Lec_BL <- sapply(calidad_BL, función_qualxLec_BL)
qual_Lec_BL[,1]
qual_LEC_BL_df <- bind_rows(qual_Lec_BL[,1], qual_Lec_BL[,2], qual_Lec_BL[,3], qual_Lec_BL[,4],
                            qual_Lec_BL[,5], qual_Lec_BL[,6], qual_Lec_BL[,7], qual_Lec_BL[,8],
                            qual_Lec_BL[,9], qual_Lec_BL[,10],qual_Lec_BL[,11],qual_Lec_BL[,12])
qual_LEC_BL_df

qual_LEC_BL_df_1 <- qual_LEC_BL_df %>% 
  mutate(Tag = c(rep("BR1T1", 512), rep("BR1T2", 512), rep("BR1T3", 512), rep("BR1T4", 512),  
                 rep("BR2T1", 512), rep("BR2T2", 512), rep("BR2T3", 512), rep("BR2T4", 512),
                 rep("BR3T1", 512), rep("BR3T2", 512), rep("BR3T3", 512), rep("BR3T4", 512)),
         Den_menor30 = ifelse(quality >= 30, 0, density))
Área_rechazoBL <- signif(sum(qual_LEC_BL_df_1$Den_menor30)/sum(qual_LEC_BL_df_1$density)*100, 3)
view(qual_LEC_BL_df_1)

ggplot(qual_LEC_BL_df_1, aes(quality)) +geom_density(fill = "lightblue") +
  geom_vline(xintercept = 30) + 
  labs(y ="Proporción (%)", x ="Calidad (PQ)", title = "Calidad de Lecturas en BL323", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  geom_text(x=29.5, y=0.04, label = paste(Área_rechazoBL, "%"), size = 3.5)
#queda pendiente título del gráfico. 

# Filtrar BC391 -----------------------------------------------------------
#el filtrado elimina lecturas de promedio de calidad menor a 30 y remueve nucleótidos no específicos
#opciones a mejora 14 PQ -> N, codigo por 3 nucleotidos, trimTails.
ListarchBCD14
BCD14_1_filt <- funcion_lecturasmásde30ylimpieza(arch1BCD14)
BCD14_2_filt <- funcion_lecturasmásde30ylimpieza(arch2BCD14)
BCD14_3_filt <- funcion_lecturasmásde30ylimpieza(arch3BCD14)
BCD14_4_filt <- funcion_lecturasmásde30ylimpieza(arch4BCD14)
BCD14_5_filt <- funcion_lecturasmásde30ylimpieza(arch5BCD14)
BCD14_6_filt <- funcion_lecturasmásde30ylimpieza(arch6BCD14)
BCD14_7_filt <- funcion_lecturasmásde30ylimpieza(arch7BCD14)
BCD14_8_filt <- funcion_lecturasmásde30ylimpieza(arch8BCD14)
BCD14_9_filt <- funcion_lecturasmásde30ylimpieza(arch9BCD14)
BCD14_10_filt <- funcion_lecturasmásde30ylimpieza(arch10BCD14)
BCD14_11_filt <- funcion_lecturasmásde30ylimpieza(arch11BCD14)
BCD14_12_filt <- funcion_lecturasmásde30ylimpieza(arch12BCD14)



#Revisualizar quality percycle (postergado; debería impactar poco con el recorte/filtro tan esbelto)

#Escribir
writeFastq(BCD14_1_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_1filt.fastq", mode = "w")
writeFastq(BCD14_2_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_2filt.fastq", mode = "w")
writeFastq(BCD14_3_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_3filt.fastq", mode = "w")
writeFastq(BCD14_4_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_4filt.fastq", mode = "w")
writeFastq(BCD14_5_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_5filt.fastq", mode = "w")
writeFastq(BCD14_6_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_6filt.fastq", mode = "w")
writeFastq(BCD14_7_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_7filt.fastq", mode = "w")
writeFastq(BCD14_8_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_8filt.fastq", mode = "w")
writeFastq(BCD14_9_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_9filt.fastq", mode = "w")
writeFastq(BCD14_10_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_10filt.fastq", mode = "w")
writeFastq(BCD14_11_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_11filt.fastq", mode = "w")
writeFastq(BCD14_12_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14/BCD14_12filt.fastq", mode = "w")


# Filtrar BL323 -----------------------------------------------------------

#Filtrar BL323
#filterfastq, srFilter, trimTails, narrow, tables
ListarchBLD14
archivosfiltradosBLD14 <- sapply(ListarchBLD14, funcion_lecturasmásde30ylimpieza)

BLD14_1_filt <- funcion_lecturasmásde30ylimpieza(arch1BLD14)
BLD14_2_filt <- funcion_lecturasmásde30ylimpieza(arch2BLD14)
BLD14_3_filt <- funcion_lecturasmásde30ylimpieza(arch3BLD14)
BLD14_4_filt <- funcion_lecturasmásde30ylimpieza(arch4BLD14)
BLD14_5_filt <- funcion_lecturasmásde30ylimpieza(arch5BLD14)
BLD14_6_filt <- funcion_lecturasmásde30ylimpieza(arch6BLD14)
BLD14_7_filt <- funcion_lecturasmásde30ylimpieza(arch7BLD14)
BLD14_8_filt <- funcion_lecturasmásde30ylimpieza(arch8BLD14)
BLD14_9_filt <- funcion_lecturasmásde30ylimpieza(arch9BLD14)
BLD14_10_filt <- funcion_lecturasmásde30ylimpieza(arch10BLD14)
BLD14_11_filt <- funcion_lecturasmásde30ylimpieza(arch11BLD14)
BLD14_12_filt <- funcion_lecturasmásde30ylimpieza(arch12BLD14)


#Revisualizar por ciclo

#Escribir
 writeFastq(BCD14_1_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_1filt.fastq", mode = "w")
 writeFastq(BCD14_2_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_2filt.fastq", mode = "w")
 writeFastq(BCD14_3_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_3filt.fastq", mode = "w")
 writeFastq(BCD14_4_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_4filt.fastq", mode = "w")
 writeFastq(BCD14_5_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_5filt.fastq", mode = "w")
 writeFastq(BCD14_6_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_6filt.fastq", mode = "w")
 writeFastq(BCD14_7_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_7filt.fastq", mode = "w")
 writeFastq(BCD14_8_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_8filt.fastq", mode = "w")
 writeFastq(BCD14_9_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_9filt.fastq", mode = "w")
writeFastq(BCD14_10_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_10filt.fastq", mode = "w")
writeFastq(BCD14_11_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_11filt.fastq", mode = "w")
writeFastq(BCD14_12_filt, file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14/BLD14_12filt.fastq", mode = "w")








# Replicación por Pares ---------------------------------------------------
#Replicación con lecturas por pares. 2. 
#Rutas
pathBCd14_2 <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/BC d14 2"
pathBLd14_2 <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/BL d14 2"

filesBCd14_2 <- list.files(pathBCd14_2)
filesBLd14_2 <- list.files(pathBLd14_2)
filesBCd14_2
filesBLd14_2

BCD14_1_12_dirs_2 <- character()
BCD14_1_12_dirs_2 <- paste0(pathBCd14_2, "/", filesBCd14_2)
BCD14_1_12_dirs_2[1]
BCD14_1_12_dirs_2

BLD14_1_12_dirs_2 <- character()
BLD14_1_12_dirs_2 <- paste0(pathBLd14_2, "/", filesBLd14_2)

countarch1BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[1]) #para la función countFastq se utiliza la dirección de los archivos con la terminación .fastq ergo las variables con solo un guión bajo ej. dirBCD14_1
countarch2BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[2])
countarch3BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[3])
countarch4BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[4])
countarch5BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[5])
countarch6BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[6])
countarch7BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[7])
countarch8BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[8])
countarch9BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[9])
countarch10BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[10])
countarch11BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[11])
countarch12BCD14_2 <- countFastq(dirPath = BCD14_1_12_dirs_2[12])

countarchBCD14list_2 <- list(countarch1BCD14_2,  countarch2BCD14_2,  countarch3BCD14_2,  countarch4BCD14_2,
                             countarch5BCD14_2,  countarch6BCD14_2,  countarch7BCD14_2,  countarch8BCD14_2,
                             countarch9BCD14_2, countarch10BCD14_2, countarch11BCD14_2, countarch12BCD14_2)
print(countarchBCD14list_2)
view(countarchBCD14list_2) 

tibblecountBCD14_2 <- countarchBCD14list_2 %>%
  map_dfr(~.x)
view(tibblecountBCD14_2)#organizado por hileras y columnas
rnBC_2 <- row.names(tibblecountBCD14_2)

dfcountBCD14_E1_2 <- tibblecountBCD14_2 %>% 
  rename(
    Lecturas = "records", 
    Nucleotidos = "nucleotides",
    Puntuaciones = "scores"
  ) %>% 
  select(- Puntuaciones) %>%
  mutate(Archivos = rnBC_2) %>%
  relocate(Archivos, .before = everything())
str(dfcountBCD14_E1_2)
typeof(dfcountBCD14_E1_2)
view(dfcountBCD14_E1_2)

groupBR <- rep(c("BR1", "BR2", "BR3"), each = 4)

dfBCg1_2 <- dfcountBCD14_E1_2 %>%
  select(-Nucleotidos) %>%
  mutate(Grupo = groupBR) %>% 
  transmute(Archivos = c("BR1T1", "BR1T2", "BR1T3", "BR1T4", "BR2T1", "BR2T2", "BR2T3", "BR2T4","BR3T1", "BR3T2", "BR3T3", "BR3T4"), 
            Lecturas, 
            Grupo)#Añadir separación por replicados biológicos para añadir color.
view(dfBCg1_2)
dfBCg2_2 <- dfcountBCD14_E1_2 %>% 
  select(-Lecturas) %>%
  mutate(Grupo = groupBR) %>% 
  transmute(Archivos = c("BR1T1", "BR1T2", "BR1T3", "BR1T4", "BR2T1", "BR2T2", "BR2T3", "BR2T4","BR3T1", "BR3T2", "BR3T3", "BR3T4"), 
            Nucleotidos, 
            Grupo)

#plots BC_2
summary(dfBCg1_2) #20MB a 60.7 MB
summary(dfBCg2_2) #278K a 822K lecturas. 


ggplot(dfBCg1_2,aes(x= Archivos, y = Lecturas/1000, fill = Grupo)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.x = element_blank()) +
  labs(y ="Lecturas (miles)", x ="Archivos (BC391_2)", title = "Poder de Secuenciación entre Archivos BC391_2 Día 14", subtitle = "Número de Lecturas por Archivo") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ 
  scale_y_continuous(breaks = seq(0, 1000, by = 100), label= seq(0, 1000, by = 100)) +
  coord_cartesian(ylim = c(0, 1000), xlim = c(1,12))

ggplot(dfBCg2_2,aes(x= Archivos, y = Nucleotidos/1000000, fill = Grupo)) +
  geom_col() +
  labs(y ="Nucleótidos (Mb)", x ="Archivos (BC391_2)", title ="Megabases por Archivo de BC391_2 al Día 14") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

#replica para BL. Gráficos de lecturas y nucleotidos por archivo. 
countarch1BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[1]) #para la función countFastq se utiliza la dirección de los archivos con la terminación .fastq ergo las variables con solo un guión bajo ej. dirBCD14_1
countarch2BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[2])
countarch3BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[3])
countarch4BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[4])
countarch5BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[5])
countarch6BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[6])
countarch7BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[7])
countarch8BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[8])
countarch9BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[9])
countarch10BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[10])
countarch11BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[11])
countarch12BLD14_2 <- countFastq(dirPath = BLD14_1_12_dirs_2[12])

countarchBLD14list_2 <- list(countarch1BLD14_2,  countarch2BLD14_2,  countarch3BLD14_2,  countarch4BLD14_2,
                             countarch5BLD14_2,  countarch6BLD14_2,  countarch7BLD14_2,  countarch8BLD14_2,
                             countarch9BLD14_2, countarch10BLD14_2, countarch11BLD14_2, countarch12BLD14_2)
print(countarchBLD14list_2)
view(countarchBLD14list_2) 

tibblecountBLD14_2 <- countarchBLD14list_2 %>%
  map_dfr(~.x)
view(tibblecountBLD14_2)#organizado por hileras y columnas
rnBL_2 <- row.names(tibblecountBLD14_2)

dfcountBLD14_E1_2 <- tibblecountBLD14_2 %>% 
  rename(
    Lecturas = "records", 
    Nucleotidos = "nucleotides",
    Puntuaciones = "scores"
  ) %>% 
  select(- Puntuaciones) %>%
  mutate(Archivos = rnBL_2) %>%
  relocate(Archivos, .before = everything())
str(dfcountBLD14_E1_2)
typeof(dfcountBLD14_E1_2)
view(dfcountBLD14_E1_2)

groupBR <- rep(c("BR1", "BR2", "BR3"), each = 4)
dfBLg1_2 <- dfcountBLD14_E1_2 %>%
  select(-Nucleotidos) %>%
  mutate(Grupo = groupBR) %>% 
  transmute(Archivos = c("BR1T1", "BR1T2", "BR1T3", "BR1T4", "BR2T1", "BR2T2", "BR2T3", "BR2T4","BR3T1", "BR3T2", "BR3T3", "BR3T4"), 
            Lecturas, 
            Grupo)#Añadir separación por replicados biológicos para añadir color.
view(dfBCg1)
dfBLg2_2 <- dfcountBLD14_E1_2 %>% 
  select(-Lecturas) %>%
  mutate(Grupo = groupBR) %>% 
  transmute(Archivos = c("BR1T1", "BR1T2", "BR1T3", "BR1T4", "BR2T1", "BR2T2", "BR2T3", "BR2T4","BR3T1", "BR3T2", "BR3T3", "BR3T4"), 
            Nucleotidos, 
            Grupo) 
summary(dfBLg1_2)
summary(dfBLg2_2)


ggplot(dfBLg1_2,aes(x= Archivos, y = Lecturas/1000, fill = Grupo)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.x = element_blank()) +
  labs(y ="Lecturas (miles)", x ="Archivos (BL323_2)", title = "Número de Lecturas entre Archivos BL323_2")+#, subtitle = "Número de Lecturas por Archivo") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+ 
  scale_y_continuous(breaks = seq(0, 1000, by = 100), label= seq(0, 1000, by = 100)) +
  coord_cartesian(ylim = c(0, 1000), xlim = c(1,12))

ggplot(dfBLg2_2,aes(x= Archivos, y = Nucleotidos/1000000, fill = Grupo)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(y ="Nucleótidos (Mb)", x ="Archivos (BL323_2)", title = "Megabases por Archivo BL323_2") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))

#Readfastq BC_2
 arch1BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[1])
 arch2BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[2])
 arch3BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[3])
 arch4BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[4])
 arch5BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[5])
 arch6BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[6])
 arch7BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[7])
 arch8BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[8])
 arch9BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[9])
arch10BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[10])
arch11BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[11])
arch12BCD14_2 <- readFastq(dirPath = BCD14_1_12_dirs_2[12])

ListarchBCD14_2 <- list(arch1BCD14_2,  arch2BCD14_2, arch3BCD14_2, arch4BCD14_2,
                        arch5BCD14_2,  arch6BCD14_2, arch7BCD14_2, arch8BCD14_2, 
                        arch9BCD14_2, arch10BCD14_2,arch11BCD14_2,arch12BCD14_2)
ListarchBCD14_2[1]

fun_longitud_lecturasBC_2 <- function(x){
  x1 <- width(x)
  x2 <- as.data.frame(table(x1))
  x3 <- x2 %>% 
    mutate(Porcentaje = signif(Freq/sum(Freq)*100, digits = 2),
           SumaPorcientos = round(cumsum(Porcentaje), digits = 2),
           Porcentajeg = ifelse(Porcentaje >= 5, paste(Porcentaje, "%"),  ""),
           Color = ifelse(Porcentaje > 5, "Highlight", "Normal")) %>%
    arrange(desc(x1))
}
RepsporLecsBC_2 <- sapply(ListarchBCD14_2, FUN = fun_longitud_lecturasBC_2)

RepsporLecsBC_2 

df_fw_BC_2 <- bind_rows(RepsporLecsBC_2[,1], RepsporLecsBC_2[,2], RepsporLecsBC_2[,3], RepsporLecsBC_2[,4],
                        RepsporLecsBC_2[,5], RepsporLecsBC_2[,6], RepsporLecsBC_2[,7], RepsporLecsBC_2[,8],
                        RepsporLecsBC_2[,9], RepsporLecsBC_2[,10],RepsporLecsBC_2[,11],RepsporLecsBC_2[,12])
view(df_fw_BC_2)

df_BC_conid_2 <- df_fw_BC_2 %>% 
  mutate(Archivo = c(rep("BR1T1_2", 27), rep("BR1T2_2", 27), rep("BR1T3_2", 27), rep("BR1T4_2", 27),
                     rep("BR2T1_2", 27), rep("BR2T2_2", 27), rep("BR2T3_2", 27), rep("BR2T4_2", 27),
                     rep("BR3T1_2", 27), rep("BR3T2_2", 27), rep("BR3T3_2", 27), rep("BR3T4_2", 27)))
view(df_BC_conid_2)

#visualizar facet wrap

df_hist_width_BC_2 <- df_BC_conid_2 %>% 
  select(x1, Freq, Archivo, Porcentajeg, Color) 
view(df_hist_width_BC)

ggplot(df_hist_width_BC_2,aes(x= x1, y = Freq/1000, fill = Color)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.x = element_blank(), legend.position = "none") +
  labs(y ="Conteo (miles)", x ="Longitud de Lecturas (pb)", title = "Longitud de Lecturas Archivos BC391_2", subtitle = "") +
  theme(plot.title = element_text(size = 30),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25, angle = 90),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))+
  facet_wrap(~Archivo) +
  scale_x_discrete(breaks = seq(50, 80, by = 5), label= seq(50, 80, by = 5)) +
  geom_text(aes(label = paste(Porcentajeg)), vjust = -0.3, hjust= 0.8, size = 4) +
  scale_fill_manual(values = c("Normal" = "gray", "Highlight" = "red"))


#Replica readfastq BL_2
 arch1BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[1])
 arch2BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[2])
 arch3BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[3])
 arch4BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[4])
 arch5BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[5])
 arch6BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[6])
 arch7BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[7])
 arch8BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[8])
 arch9BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[9])
arch10BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[10])
arch11BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[11])
arch12BLD14_2 <- readFastq(dirPath = BLD14_1_12_dirs_2[12])

ListarchBLD14_2 <- list(arch1BLD14_2,  arch2BLD14_2, arch3BLD14_2, arch4BLD14_2,
                        arch5BLD14_2,  arch6BLD14_2, arch7BLD14_2, arch8BLD14_2, 
                        arch9BLD14_2, arch10BLD14_2,arch11BLD14_2,arch12BLD14_2)
ListarchBLD14_2[1]

fun_longitud_lecturasBL_2 <- function(x){
  x1 <- width(x)
  x2 <- as.data.frame(table(x1))
  x3 <- x2 %>% 
    mutate(Porcentaje = signif(Freq/sum(Freq)*100, digits = 2),
           SumaPorcientos = round(cumsum(Porcentaje), digits = 2),
           Porcentajeg = ifelse(Porcentaje >= 5, paste(Porcentaje, "%"),  ""),
           Color = ifelse(Porcentaje > 5, "Highlight", "Normal")) %>%
    arrange(desc(x1))
}

RepsporLecsBL_2 <- sapply(ListarchBLD14_2, FUN = fun_longitud_lecturasBL_2)
RepsporLecsBL_2

df_fw_BL_2 <- bind_rows(RepsporLecsBL_2[,1], RepsporLecsBL_2[,2], RepsporLecsBL_2[,3], RepsporLecsBL_2[,4],
                        RepsporLecsBL_2[,5], RepsporLecsBL_2[,6], RepsporLecsBL_2[,7], RepsporLecsBL_2[,8],
                        RepsporLecsBL_2[,9], RepsporLecsBL_2[,10],RepsporLecsBL_2[,11],RepsporLecsBL_2[,12])
df_fw_BL_2
df_BL_conid_2 <- df_fw_BL_2 %>% 
  mutate(Archivo = c(rep("BR1T1_2", 27), rep("BR1T2_2", 27), rep("BR1T3_2", 27), rep("BR1T4_2", 27),
                     rep("BR2T1_2", 27), rep("BR2T2_2", 27), rep("BR2T3_2", 27), rep("BR2T4_2", 27),
                     rep("BR3T1_2", 27), rep("BR3T2_2", 27), rep("BR3T3_2", 27), rep("BR3T4_2", 27)))

df_hist_width_BL_2 <- df_BL_conid_2 %>% 
  select(x1, Freq, Archivo, Porcentajeg, Color)
view(df_hist_width_BL)
#Visualizar con Facet wrap
#longitud de lecturas
ggplot(df_hist_width_BL_2,aes(x= x1, y = Freq/1000, fill = Color)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.x = element_blank(), legend.position = "none") +
  labs(y ="Conteo (miles)", x ="Longitud de Lecturas (pb)", title = "Longitud de Lecturas Archivos BL323_2", subtitle = "") +
  theme(plot.title = element_text(size = 30),
        plot.subtitle = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))+
  facet_wrap(~Archivo) +
  scale_x_discrete(breaks = seq(50, 80, by = 5), label= seq(50, 80, by = 5)) +
  geom_text(aes(label = paste(Porcentajeg)), vjust = -0.3, hjust= 0.8, size = 4) +
  scale_fill_manual(values = c("Normal" = "gray", "Highlight" = "red"))


# QA_2 --------------------------------------------------------------------
#Por pares estos son los .fastq2. 
#Calidad de Bases nucleotídicas y calidad de BN por ciclo.

#BC_2 qa
qa_BC_BR1T1_2 <-  qa(arch1BCD14_2, lane = 1)
qa_BC_BR1T2_2 <-  qa(arch2BCD14_2, lane = 1)
qa_BC_BR1T3_2 <-  qa(arch3BCD14_2, lane = 1)
qa_BC_BR1T4_2 <-  qa(arch4BCD14_2, lane = 1)
qa_BC_BR2T1_2 <-  qa(arch5BCD14_2, lane = 1)
qa_BC_BR2T2_2 <-  qa(arch6BCD14_2, lane = 1)
qa_BC_BR2T3_2 <-  qa(arch7BCD14_2, lane = 1)
qa_BC_BR2T4_2 <-  qa(arch8BCD14_2, lane = 1)
qa_BC_BR3T1_2 <-  qa(arch9BCD14_2, lane = 1)
qa_BC_BR3T2_2 <- qa(arch10BCD14_2, lane = 1)
qa_BC_BR3T3_2 <- qa(arch11BCD14_2, lane = 1)
qa_BC_BR3T4_2 <- qa(arch12BCD14_2, lane = 1)

calidad_BC_2 <- c(qa_BC_BR1T1_2, qa_BC_BR1T2_2, qa_BC_BR1T3_2, qa_BC_BR1T4_2,
                  qa_BC_BR2T1_2, qa_BC_BR2T2_2, qa_BC_BR2T3_2, qa_BC_BR2T4_2,
                  qa_BC_BR3T1_2, qa_BC_BR3T2_2, qa_BC_BR3T3_2, qa_BC_BR3T4_2)
calidad_BC_2

función_BC_cal_2 <- function(x){
  x1 <- x[["baseQuality"]] %>%
    select(score, count) %>% 
    mutate(Puntaje = -1:93)
  x2 <- x1[-1,] %>%
    mutate(Porcentaje = signif(((count/sum(count))*100), digits = 2)) %>% 
    select(-score) %>% 
    filter(count > 28)
}
listofdfs_BC_qual_2 <- sapply(calidad_BC_2, función_BC_cal_2)
listofdfs_BC_qual_2 


listofdfs_BC_qual_g1_2<- bind_rows(listofdfs_BC_qual_2[,1], listofdfs_BC_qual_2[,2], listofdfs_BC_qual_2[,3], listofdfs_BC_qual_2[,4],
                                   listofdfs_BC_qual_2[,5], listofdfs_BC_qual_2[,6], listofdfs_BC_qual_2[,7], listofdfs_BC_qual_2[,8],
                                   listofdfs_BC_qual_2[,9], listofdfs_BC_qual_2[,10],listofdfs_BC_qual_2[,11],listofdfs_BC_qual_2[,12])
listofdfs_BC_qual_g1_2
listofdfs_BC_qual_g1_fnot0_2 <- listofdfs_BC_qual_g1_2 %>% 
  filter(Puntaje > 13) %>% 
  mutate(Archivo = c(rep("BR1T1_2", 5), rep("BR1T2_2", 5), rep("BR1T3_2", 5), rep("BR1T4_2", 5),
                     rep("BR2T1_2", 5), rep("BR2T2_2", 5), rep("BR2T3_2", 5), rep("BR2T4_2", 5),
                     rep("BR3T1_2", 5), rep("BR3T2_2", 5), rep("BR3T3_2", 5), rep("BR3T4_2", 5)), 
         Conteo = count,
         LosPorcentajesChidos = ifelse(Porcentaje >= 12, paste(Porcentaje, "%"), ""),
         Color = ifelse(LosPorcentajesChidos >= 12, "Highlight", "Normal")) %>% #agregar id para facet wrap  
  select(-count) %>% 
  relocate(Archivo, Puntaje, Conteo, Porcentaje) 
view(listofdfs_BC_qual_g1_fnot0)

ggplot(listofdfs_BC_qual_g1_fnot0_2, aes(x = Puntaje , y= Conteo/1000000,)) + #fill = group)) +
  geom_col() +
  labs(y ="Conteo (millones)", x ="Calidad (PQ)", title = "Conteo de Calidades de Bases Nucleotídicas BC391_2", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = c(14,21, 27, 32, 36), label= c(14, 21, 27, 32, 36)) +
  theme(axis.text.y = element_text(angle = 0), legend.position = "none") +
  scale_y_continuous(labels = comma, limits = c(0, 50, by = 5)) +
  geom_text(aes(label = paste(LosPorcentajesChidos)), vjust = -0.3, size = 3) +
  facet_wrap(~Archivo)#scale_fill_manual(values = c("normal" = "gray", "highlight" = "red")) 

fun_BC_percycle_2 <- function(x){
  x1 <- x[["perCycle"]]
  x2 <- x1$quality %>% 
    select(-Quality, -lane) %>% 
    mutate(Puntajeacolor =  ifelse(Score == 14, "Mala Calidad",
                                   ifelse(Score == 21, "Baja Calidad",
                                          ifelse(Score == 27, "Calidad Neutra",
                                                 ifelse(Score == 32, "Calidad Aceptable", "Calidad Excelente")))))
}

listdfs_BC_percycle_g1_2 <- sapply(calidad_BC_2, fun_BC_percycle_2)
listdfs_BC_percycle_g1_2
#411, 415, 414, 407, 376, 376, 388, 386, 377, 376, 384, 383
sum(411, 415, 414, 407, 376, 376, 388, 386, 377, 376, 384, 383
)
#visualizar por ciclos BC391_2
BC_percycle_2 <- bind_rows(listdfs_BC_percycle_g1_2[,1], listdfs_BC_percycle_g1_2[,2], listdfs_BC_percycle_g1_2[,3],
                           listdfs_BC_percycle_g1_2[,4], listdfs_BC_percycle_g1_2[,5], listdfs_BC_percycle_g1_2[,6],
                           listdfs_BC_percycle_g1_2[,7], listdfs_BC_percycle_g1_2[,8], listdfs_BC_percycle_g1_2[,9],
                           listdfs_BC_percycle_g1_2[,10],listdfs_BC_percycle_g1_2[,11],listdfs_BC_percycle_g1_2[,12])
BC_percycle_2

BC_percycle_g1_2 <- BC_percycle_2 %>%
  mutate(Archivo = c(rep("BR1T1_2", 411), rep("BR1T2_2", 415), rep("BR1T3_2", 414), rep("BR1T4_2", 407),
                     rep("BR2T1_2", 376), rep("BR2T2_2", 376), rep("BR2T3_2", 388), rep("BR2T4_2", 386),
                     rep("BR3T1_2", 377), rep("BR3T2_2", 376), rep("BR3T3_2", 384), rep("BR3T4_2", 383)))
BC_percycle_g1_2$Puntajeacolor <- factor(BC_percycle_g1_2$Puntajeacolor, levels = c("Mala Calidad", "Baja Calidad", "Calidad Neutra", "Calidad Aceptable", "Calidad Excelente"))
view(BC_percycle_g1_2)

ggplot(BC_percycle_g1_2, aes(x = Cycle, y = Count/1000, fill = Puntajeacolor)) +
  geom_col(position = "stack") + 
  scale_fill_manual(values = c("Mala Calidad" = "darkred", "Baja Calidad" = "red", "Calidad Neutra" = "lightskyblue", "Calidad Aceptable" = "lightgreen", "Calidad Excelente"= "green"))+
  labs(y ="Conteo (miles)", x ="Ciclo (Posición Nt)", title = "Conteo de Calidad de Nucleótidos por Ciclo en BC391_2", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = seq(0, 75, by = 15), label= seq(0, 75, by = 15)) +
  scale_y_continuous(breaks = seq(0, 700, by = 100), label= seq(0, 700, by = 100))+
  facet_wrap(~Archivo)
#cambiar eje x, liberar y, ordenar por factor de puntaje a color

#Filtrar 
#filterfastq, srFilter, trimTails, narrow, tables

#Revisualizar por ciclo

#Escribir
#writeFastq


#BL_2 qa
qa_BL_BR1T1_2 <-  qa(arch1BLD14_2, lane = 1)
qa_BL_BR1T2_2 <-  qa(arch2BLD14_2, lane = 1)
qa_BL_BR1T3_2 <-  qa(arch3BLD14_2, lane = 1)
qa_BL_BR1T4_2 <-  qa(arch4BLD14_2, lane = 1)
qa_BL_BR2T1_2 <-  qa(arch5BLD14_2, lane = 1)
qa_BL_BR2T2_2 <-  qa(arch6BLD14_2, lane = 1)
qa_BL_BR2T3_2 <-  qa(arch7BLD14_2, lane = 1)
qa_BL_BR2T4_2 <-  qa(arch8BLD14_2, lane = 1)
qa_BL_BR3T1_2 <-  qa(arch9BLD14_2, lane = 1)
qa_BL_BR3T2_2 <- qa(arch10BLD14_2, lane = 1)
qa_BL_BR3T3_2 <- qa(arch11BLD14_2, lane = 1)
qa_BL_BR3T4_2 <- qa(arch12BLD14_2, lane = 1)

calidad_BL_2 <- c(qa_BL_BR1T1_2, qa_BL_BR1T2_2, qa_BL_BR1T3_2, qa_BL_BR1T4_2,
                  qa_BL_BR2T1_2, qa_BL_BR2T2_2, qa_BL_BR2T3_2, qa_BL_BR2T4_2,
                  qa_BL_BR3T1_2, qa_BL_BR3T2_2, qa_BL_BR3T3_2, qa_BL_BR3T4_2)
calidad_BL_2

función_BL_cal_2 <- function(x){
  x1 <- x[["baseQuality"]] %>%
    select(score, count) %>% 
    mutate(Puntaje = -1:93)
  x2 <- x1[-1,] %>%
    mutate(Porcentaje = signif(((count/sum(count))*100), digits = 2)) %>% 
    select(-score) %>% 
    filter(count > 28)
}

listofdfs_BL_qual_2 <- sapply(calidad_BL_2, función_BL_cal_2)
listofdfs_BL_qual_2 
listofdfs_BL_qual_g1_2<- bind_rows(listofdfs_BL_qual_2[,1], listofdfs_BL_qual_2[,2], listofdfs_BL_qual_2[,3], listofdfs_BL_qual_2[,4],
                                   listofdfs_BL_qual_2[,5], listofdfs_BL_qual_2[,6], listofdfs_BL_qual_2[,7], listofdfs_BL_qual_2[,8],
                                   listofdfs_BL_qual_2[,9], listofdfs_BL_qual_2[,10],listofdfs_BL_qual_2[,11],listofdfs_BL_qual_2[,12])
listofdfs_BL_qual_g1_fnot0_2 <- listofdfs_BL_qual_g1_2 %>% 
  filter(Puntaje > 13) %>%
  mutate(Archivo = c(rep("BR1T1_2", 5), rep("BR1T2_2", 5), rep("BR1T3_2", 5), rep("BR1T4_2", 5),
                     rep("BR2T1_2", 5), rep("BR2T2_2", 5), rep("BR2T3_2", 5), rep("BR2T4_2", 5),
                     rep("BR3T1_2", 5), rep("BR3T2_2", 5), rep("BR3T3_2", 5), rep("BR3T4_2", 5)),
         Conteo = count,
         LosPorcentajesChidos = ifelse(Porcentaje >= 12, paste(Porcentaje, "%"),  ""),
         Color = ifelse(LosPorcentajesChidos >= 12, "Highlight", "Normal")) %>% 
  select(-count) %>%
  relocate(Archivo, Puntaje, Conteo, Porcentaje)  

view(listofdfs_BL_qual_g1_fnot0_2)

ggplot(listofdfs_BL_qual_g1_fnot0_2, aes(x = Puntaje , y= Conteo/1000000,)) + #fill = group)) +
  geom_col() +
  labs(y ="Conteo (millones)", x ="Calidad (PQ)", title = "Conteo de Calidades de Bases Nucleotídicas BL323_2", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = c(14,21, 27, 32, 36), label= c(14, 21, 27, 32, 36)) +
  theme(axis.text.y = element_text(angle = 0), legend.position = "none") +
  scale_y_continuous(labels = comma, limits = c(0, 50, by = 5)) +
  geom_text(aes(label = paste(listofdfs_BC_qual_g1_fnot0$LosPorcentajesChidos)), vjust = -0.4, hjust = 0.3, size = 3) +
  facet_wrap(~Archivo)#scale_fill_manual(values = c("normal" = "gray", "highlight" = "red")) 

fun_BL_percycle_2 <- function(x){
  x1 <- x[["perCycle"]]
  x2 <- x1$quality %>% 
    select(-Quality, -lane) %>% 
    mutate(Puntajeacolor =  ifelse(Score == 14, "Mala Calidad",
                                   ifelse(Score == 21, "Baja Calidad",
                                          ifelse(Score == 27, "Calidad Neutra",
                                                 ifelse(Score == 32, "Calidad Aceptable", "Calidad Excelente")))))
}

listdfs_BL_percycle_g1_2 <- sapply(calidad_BL_2, fun_BL_percycle_2)
listdfs_BL_percycle_g1_2
#404, 396, 407, 399, 376, 376, 381, 385, 377, 376, 384, 381

#visualizar por ciclos BL323_2
BL_percycle_2 <- bind_rows(listdfs_BL_percycle_g1_2[,1], listdfs_BL_percycle_g1_2[,2], listdfs_BL_percycle_g1_2[,3], 
                           listdfs_BL_percycle_g1_2[,4], listdfs_BL_percycle_g1_2[,5], listdfs_BL_percycle_g1_2[,6],
                           listdfs_BL_percycle_g1_2[,7], listdfs_BL_percycle_g1_2[,8], listdfs_BL_percycle_g1_2[,9],
                           listdfs_BL_percycle_g1_2[,10],listdfs_BL_percycle_g1_2[,11],listdfs_BL_percycle_g1_2[,12])


BL_percycle_g1_2 <- BL_percycle_2 %>%
  mutate(Archivo = c(rep("BR1T1_2", 404), rep("BR1T2_2", 396), rep("BR1T3_2", 407), rep("BR1T4_2", 399),
                     rep("BR2T1_2", 376), rep("BR2T2_2", 376), rep("BR2T3_2", 381), rep("BR2T4_2", 385),
                     rep("BR3T1_2", 377), rep("BR3T2_2", 376), rep("BR3T3_2", 384), rep("BR3T4_2", 381))) 
BL_percycle_g1_2$Puntajeacolor <- factor(BL_percycle_g1_2$Puntajeacolor, levels = c("Mala Calidad", "Baja Calidad", "Calidad Neutra", "Calidad Aceptable", "Calidad Excelente"))
BL_percycle_g1_2

ggplot(BL_percycle_g1_2, aes(x = Cycle, y = Count/1000, fill = Puntajeacolor)) +
  geom_col(position = "stack") + 
  scale_fill_manual(values = c("Mala Calidad" = "darkred", "Baja Calidad" = "red", "Calidad Neutra" = "lightskyblue", "Calidad Aceptable" = "lightgreen", "Calidad Excelente"= "green"))+
  labs(y ="Conteo (miles)", x ="Ciclo (Posición Nt)", title = "Conteo de Calidad de nucleótidos por Ciclo en BL323_2", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = seq(0, 75, by = 15), label= seq(0, 75, by = 15)) +
  scale_y_continuous(breaks = seq(0, 500, by = 100), label= seq(0, 500, by = 100))+
  facet_wrap(~Archivo)

#Filtrar
#filterfastq, srFilter, trimTails, narrow, tables

# Filtrar BC rev --------------------------------------------------------------
#sapply y rev

ListarchBCD14_2

read_more_30 <- function(x){
  x0<- rowMeans(as(quality(x), "matrix"), na.rm=T) >= 30
  return(x[x0])#x4<- x[x2(x)]
  #writeFastq()
  #filtered <- x[Read_more30(x)]
}


funcion_lecturasmásde30ylimpieza <- function(x){
  return(clean(read_more_30(x)))
}



#datosfiltradosBCD14_2 <- sapply(ListarchBCD14_2, funcion_lecturasmásde30ylimpieza)
  BCD14_1_filt_2 <- arch1BCD14_2 #funcion_lecturasmásde30ylimpieza(arch1BCD14_2)
  BCD14_2_filt_2 <- arch2BCD14_2 #funcion_lecturasmásde30ylimpieza(arch2BCD14_2)
  BCD14_3_filt_2 <- arch3BCD14_2 #funcion_lecturasmásde30ylimpieza(arch3BCD14_2)
  BCD14_4_filt_2 <- arch4BCD14_2 #funcion_lecturasmásde30ylimpieza(arch4BCD14_2)
  BCD14_5_filt_2 <- arch5BCD14_2 #funcion_lecturasmásde30ylimpieza(arch5BCD14_2)
  BCD14_6_filt_2 <- arch6BCD14_2 #funcion_lecturasmásde30ylimpieza(arch6BCD14_2)
  BCD14_7_filt_2 <- arch7BCD14_2 #funcion_lecturasmásde30ylimpieza(arch7BCD14_2)
  BCD14_8_filt_2 <- arch8BCD14_2 #funcion_lecturasmásde30ylimpieza(arch8BCD14_2)
  BCD14_9_filt_2 <- arch9BCD14_2 #funcion_lecturasmásde30ylimpieza(arch9BCD14_2)
BCD14_10_filt_2 <-  arch10BCD14_2 #funcion_lecturasmásde30ylimpieza(arch10BCD14_2)
BCD14_11_filt_2 <-  arch11BCD14_2 #funcion_lecturasmásde30ylimpieza(arch11BCD14_2)
BCD14_12_filt_2 <-  arch12BCD14_2 #funcion_lecturasmásde30ylimpieza(arch12BCD14_2)

#Reverse complement
complementoreverso1 <- ShortReadQ(
               id = ShortRead::id(BCD14_1_filt_2),
  sread = reverseComplement(sread(BCD14_1_filt_2)), 
                quality = quality(BCD14_1_filt_2)
)
complementoreverso2 <- ShortReadQ(
               id = ShortRead::id(BCD14_2_filt_2),
  sread = reverseComplement(sread(BCD14_2_filt_2)), 
                quality = quality(BCD14_2_filt_2)
)
complementoreverso3 <- ShortReadQ(
               id = ShortRead::id(BCD14_3_filt_2),
  sread = reverseComplement(sread(BCD14_3_filt_2)), 
                quality = quality(BCD14_3_filt_2)
)
complementoreverso4 <- ShortReadQ(
               id = ShortRead::id(BCD14_4_filt_2),
  sread = reverseComplement(sread(BCD14_4_filt_2)), 
                quality = quality(BCD14_4_filt_2)
)
complementoreverso5 <- ShortReadQ(
               id = ShortRead::id(BCD14_5_filt_2),
  sread = reverseComplement(sread(BCD14_5_filt_2)), 
                quality = quality(BCD14_5_filt_2)
)
complementoreverso6 <- ShortReadQ(
  id = ShortRead::id(BCD14_6_filt_2),
  sread = reverseComplement(sread(BCD14_6_filt_2)), 
  quality = quality(BCD14_6_filt_2)
)
complementoreverso7 <- ShortReadQ(
  id = ShortRead::id(BCD14_7_filt_2),
  sread = reverseComplement(sread(BCD14_7_filt_2)), 
  quality = quality(BCD14_7_filt_2)
)
complementoreverso8 <- ShortReadQ(
  id = ShortRead::id(BCD14_8_filt_2),
  sread = reverseComplement(sread(BCD14_8_filt_2)), 
  quality = quality(BCD14_8_filt_2)
)
complementoreverso9 <- ShortReadQ(
  id = ShortRead::id(BCD14_9_filt_2),
  sread = reverseComplement(sread(BCD14_9_filt_2)), 
  quality = quality(BCD14_9_filt_2)
)
complementoreverso10 <- ShortReadQ(
  id = ShortRead::id(BCD14_10_filt_2),
  sread = reverseComplement(sread(BCD14_10_filt_2)), 
  quality = quality(BCD14_10_filt_2)
)
complementoreverso11 <- ShortReadQ(
  id = ShortRead::id(BCD14_11_filt_2),
  sread = reverseComplement(sread(BCD14_11_filt_2)), 
  quality = quality(BCD14_11_filt_2)
)
complementoreverso12 <- ShortReadQ(
  id = ShortRead::id(BCD14_12_filt_2),
  sread = reverseComplement(sread(BCD14_12_filt_2)), 
  quality = quality(BCD14_12_filt_2)
)
#Combinados #no ha funcionado; quizá checando hacer vectores y tejer ambos archivos pero mucho código. 
 combined_fastq1BC <- ShortReadQ(
   id = c(ShortRead::id(BCD14_1_filt), ShortRead::id(complementoreverso1)),
   sread = c(sread(BCD14_1_filt), sread(complementoreverso1)), 
   quality = c(quality(BCD14_1_filt), quality(complementoreverso1)))
 combined_fastq2BC <- ShortReadQ(
   id = ShortRead::c(id(BCD14_2_filt), complementoreverso2),
   sread = c(sread(BCD14_2_filt), sread(complementoreverso2)), 
   quality = c(quality(BCD14_2_filt), quality(complementoreverso2))) 
 combined_fastq3BC <- c(BCD14_3_filt, complementoreverso3) 
 combined_fastq3BC
 
 combined_fastq4BC <- c(BCD14_4_filt, complementoreverso4) 
 combined_fastq5BC <- c(BCD14_5_filt, complementoreverso5) 
 combined_fastq6BC <- c(BCD14_6_filt, complementoreverso6) 
 combined_fastq7BC <- c(BCD14_7_filt, complementoreverso7) 
 combined_fastq8BC <- c(BCD14_8_filt, complementoreverso8) 
 combined_fastq9BC <- c(BCD14_9_filt, complementoreverso9) 
combined_fastq10BC <- c(BCD14_10_filt, complementoreverso10) 
combined_fastq11BC <- c(BCD14_11_filt, complementoreverso11) 
combined_fastq12BC <- c(BCD14_12_filt, complementoreverso12) 

#writeFastq por unión
pathcombinados <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/"
writeFastq(combined_fastq1BC, file = paste(pathcombinados, "BCD14 y 2/BCD14_1_1y2"))

#writeFastq complemento reverso de archivos 
pathBCd14_2_filt <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14_2"
  writeFastq(complementoreverso1, file = paste(pathBCd14_2_filt, "/BCD14_2_1sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso2, file = paste(pathBCd14_2_filt, "/BCD14_2_2sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso3, file = paste(pathBCd14_2_filt, "/BCD14_2_3sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso4, file = paste(pathBCd14_2_filt, "/BCD14_2_4sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso5, file = paste(pathBCd14_2_filt, "/BCD14_2_5sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso6, file = paste(pathBCd14_2_filt, "/BCD14_2_6sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso7, file = paste(pathBCd14_2_filt, "/BCD14_2_7sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso8, file = paste(pathBCd14_2_filt, "/BCD14_2_8sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso9, file = paste(pathBCd14_2_filt, "/BCD14_2_9sinfilt.fastq", sep = ""), mode = "w")
writeFastq(complementoreverso10, file = paste(pathBCd14_2_filt, "/BCD14_2_10sinfilt.fastq", sep = ""), mode = "w")
writeFastq(complementoreverso11, file = paste(pathBCd14_2_filt, "/BCD14_2_11sinfilt.fastq", sep = ""), mode = "w")
writeFastq(complementoreverso12, file = paste(pathBCd14_2_filt, "/BCD14_2_12sinfilt.fastq", sep = ""), mode = "w")

# Filtrar BL rev --------------------------------------------------------------
#sapply y rev unir. 
#ListarchBLD14_2
#datosfiltradosBLD14_2 <- sapply(ListarchBLD14_2, funcion_lecturasmásde30ylimpieza)
 BLD14_1_filt_2 <- arch1BLD14_2#funcion_lecturasmásde30ylimpieza(arch1BLD14_2)
 BLD14_2_filt_2 <- arch2BLD14_2#funcion_lecturasmásde30ylimpieza(arch2BLD14_2)
 BLD14_3_filt_2 <- arch3BLD14_2#funcion_lecturasmásde30ylimpieza(arch3BLD14_2)
 BLD14_4_filt_2 <- arch4BLD14_2#funcion_lecturasmásde30ylimpieza(arch4BLD14_2)
 BLD14_5_filt_2 <- arch5BLD14_2#funcion_lecturasmásde30ylimpieza(arch5BLD14_2)
 BLD14_6_filt_2 <- arch6BLD14_2#funcion_lecturasmásde30ylimpieza(arch6BLD14_2)
 BLD14_7_filt_2 <- arch7BLD14_2#funcion_lecturasmásde30ylimpieza(arch7BLD14_2)
 BLD14_8_filt_2 <- arch8BLD14_2#funcion_lecturasmásde30ylimpieza(arch8BLD14_2)
 BLD14_9_filt_2 <- arch9BLD14_2#funcion_lecturasmásde30ylimpieza(arch9BLD14_2)
BLD14_10_filt_2 <- arch10BLD14_2#funcion_lecturasmásde30ylimpieza(arch10BLD14_2)
BLD14_11_filt_2 <- arch11BLD14_2#funcion_lecturasmásde30ylimpieza(arch11BLD14_2)
BLD14_12_filt_2 <- arch12BLD14_2#funcion_lecturasmásde30ylimpieza(arch12BLD14_2)

#Revisualizar por ciclo

#complemento Reverso
complementoreverso1BL <- ShortReadQ(
               id = ShortRead::id(BLD14_1_filt_2),
  sread = reverseComplement(sread(BLD14_1_filt_2)), 
                quality = quality(BLD14_1_filt_2)
)

complementoreverso2BL <- ShortReadQ(
               id = ShortRead::id(BLD14_2_filt_2),
  sread = reverseComplement(sread(BLD14_2_filt_2)), 
                quality = quality(BLD14_2_filt_2)
)
complementoreverso3BL <- ShortReadQ(
               id = ShortRead::id(BLD14_3_filt_2),
  sread = reverseComplement(sread(BLD14_3_filt_2)), 
                quality = quality(BLD14_3_filt_2)
)
complementoreverso4BL <- ShortReadQ(
               id = ShortRead::id(BLD14_4_filt_2),
  sread = reverseComplement(sread(BLD14_4_filt_2)), 
                quality = quality(BLD14_4_filt_2)
)
complementoreverso5BL <- ShortReadQ(
               id = ShortRead::id(BLD14_5_filt_2),
  sread = reverseComplement(sread(BLD14_5_filt_2)), 
                quality = quality(BLD14_5_filt_2)
)
complementoreverso6BL <- ShortReadQ(
               id = ShortRead::id(BLD14_6_filt_2),
  sread = reverseComplement(sread(BLD14_6_filt_2)), 
                quality = quality(BLD14_6_filt_2)
)
complementoreverso7BL <- ShortReadQ(
               id = ShortRead::id(BLD14_7_filt_2),
  sread = reverseComplement(sread(BLD14_7_filt_2)), 
                quality = quality(BLD14_7_filt_2)
)
complementoreverso8BL <- ShortReadQ(
               id = ShortRead::id(BLD14_8_filt_2),
  sread = reverseComplement(sread(BLD14_8_filt_2)), 
                quality = quality(BLD14_8_filt_2)
)
complementoreverso9BL <- ShortReadQ(
               id = ShortRead::id(BLD14_9_filt_2),
  sread = reverseComplement(sread(BLD14_9_filt_2)), 
                quality = quality(BLD14_9_filt_2)
)
complementoreverso10BL <- ShortReadQ(
               id = ShortRead::id(BLD14_10_filt_2),
  sread = reverseComplement(sread(BLD14_10_filt_2)), 
                quality = quality(BLD14_10_filt_2)
)
complementoreverso11BL <- ShortReadQ(
               id = ShortRead::id(BLD14_11_filt_2),
  sread = reverseComplement(sread(BLD14_11_filt_2)), 
                quality = quality(BLD14_11_filt_2)
)
complementoreverso12BL <- ShortReadQ(
               id = ShortRead::id(BLD14_12_filt_2),
  sread = reverseComplement(sread(BLD14_12_filt_2)), 
                quality = quality(BLD14_12_filt_2)
)

#escribir por unión *no ha jalado
c()
 combined_fastq1BL <- c(BLD14_1_filt, complementoreverso1BL) 
 class(combined_fastq1BL)
 combined_fastq2BL <- c(BLD14_2_filt, complementoreverso2BL) 
 combined_fastq3BL <- c(BLD14_3_filt, complementoreverso3BL) 
 combined_fastq4BL <- c(BLD14_4_filt, complementoreverso4BL) 
 combined_fastq5BL <- c(BLD14_5_filt, complementoreverso5BL) 
 combined_fastq6BL <- c(BLD14_6_filt, complementoreverso6BL) 
 combined_fastq7BL <- c(BLD14_7_filt, complementoreverso7BL) 
 combined_fastq8BL <- c(BLD14_8_filt, complementoreverso8BL) 
 combined_fastq9BL <- c(BLD14_9_filt, complementoreverso9BL) 
combined_fastq10BL <- c(BLD14_10_filt, complementoreverso10BL) 
combined_fastq11BL <- c(BLD14_11_filt, complementoreverso11BL) 
combined_fastq12BL <- c(BLD14_12_filt, complementoreverso12BL) 
#writeFastq combinado

#writeFastq complemento reverso BL
pathBLd14_2_filt <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14_2"
  writeFastq(complementoreverso1BL, file = paste(pathBLd14_2_filt, "/BLD14_2_1sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso2BL, file = paste(pathBLd14_2_filt, "/BLD14_2_2sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso3BL, file = paste(pathBLd14_2_filt, "/BLD14_2_3sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso4BL, file = paste(pathBLd14_2_filt, "/BLD14_2_4sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso5BL, file = paste(pathBLd14_2_filt, "/BLD14_2_5sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso6BL, file = paste(pathBLd14_2_filt, "/BLD14_2_6sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso7BL, file = paste(pathBLd14_2_filt, "/BLD14_2_7sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso8BL, file = paste(pathBLd14_2_filt, "/BLD14_2_8sinfilt.fastq", sep = ""), mode = "w")
  writeFastq(complementoreverso9BL, file = paste(pathBLd14_2_filt, "/BLD14_2_9sinfilt.fastq", sep = ""), mode = "w")
writeFastq(complementoreverso10BL, file = paste(pathBLd14_2_filt, "/BLD14_2_10sinfilt.fastq", sep = ""), mode = "w")
writeFastq(complementoreverso11BL, file = paste(pathBLd14_2_filt, "/BLD14_2_11sinfilt.fastq", sep = ""), mode = "w")
writeFastq(complementoreverso12BL, file = paste(pathBLd14_2_filt, "/BLD14_2_12sinfilt.fastq", sep = ""), mode = "w")



# Explorativo -------------------------------------------------------------

#Por lecturas calificadas todos los datos
#Ejemplo por otros. 
#df <- qa[["readQualityScore"]]
#ShortRead:::.plotReadQuality(df[df$type=="read",])


view(qa_BC_BR1T1[["readQualityScore"]]) 
ensayo_rQS_BC1 <- qa_BC_BR1T1[["readQualityScore"]]

ensayo_BC1_f <- ensayo_rQS_BC1 %>% 
  select(-lane, -type)
ensayo_BC1_f #en notación científica. Abajo de 30 de calidad hay a la menos 5-menos 7. Se pueden retirar. 

ensayo_BC1_f1 <- ensayo_BC1_f1 %>% 
  mutate(
  )

sum(ensayo_BC1_f$density)

ggplot(ensayo_BC1_f, aes(quality)) +geom_density(fill = "lightblue")

#consultar con pato 

#what about which qualities to stay with to make the filters. Este es ejemplo de solo 1 archivo / molde
QABCD141_E1 <- qa(arch1BCD14, lane = 1) #debería extenderlo a todos los archivos
QABCD141_E1
quality(arch1BCD14)
encoding(quality(arch1BCD14))

BQ_BCD141 <- QABCD141_E1[["baseQuality"]] %>%   
  select(score, count) %>% 
  mutate(Puntaje = -1:93)
str(BQ_BCD141)
view(BQ_BCD141)
BQ_BCD141_2 <- BQ_BCD141[-1,] %>% 
  mutate(Porcentaje = signif(((count/sum(count))*100), digits = 2), 
         group = c(rep("normal", 32), "highlight", rep("normal", 3),
                   "highlight", rep("normal", 57))) %>% 
  select(-score) %>% 
  filter(count > 28)  
view(BQ_BCD141_2)

ggplot(BQ_BCD141_2, aes(x = Puntaje , y= count/1000000, fill = group)) +
  geom_col() +
  labs(y ="Conteo (millones)", x ="Calidad (PQ)", title = "Conteo de Puntajes de Calidad de Bases Nucleotídicas Archivo 1 BC391", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = c(14,21, 27, 32, 36), label= c(14, 21, 27, 32, 36)) +
  theme(axis.text.y = element_text(angle = 0), legend.position = "none") +
  scale_y_continuous(labels = comma, limits = c(0, 50, by = 5)) +
  geom_text(aes(label = paste(BQ_BCD141_2$Porcentaje, "%")), vjust = -0.1, size = 10) +
  scale_fill_manual(values = c("normal" = "gray", "highlight" = "red")) 
#El 95% de los nucleotidos fueron secuenciados con una calidad >30 es decir de 99.9% de precisión. 
#Al usar todos los datos ocuparé agruparlos por estructura similar en un marco de datos, ordenados por nombre de archivo y obtener los promedios. Interesante ver como juntar graficamente la información. 

#conteo de calidades por ciclo
PC_BCD141 <- QABCD141_E1[["perCycle"]]
head(PC_BCD141)
str(PC_BCD141)

PC_BCD141_bc <- PC_BCD141$baseCall %>% 
  select(-lane) %>% 
  rename(
    BaseNitrogenada = "Base"
  )
head(PC_BCD141_bc)
head(PC_BCD141_bc)
PC_BCD141_q <- PC_BCD141$quality %>% 
  select(-Quality, -lane) %>% 
  mutate(Puntajeacolor =  ifelse(Score == 14, "Mala Calidad",
                                 ifelse(Score == 21, "Baja Calidad",
                                        ifelse(Score == 27, "Calidad Neutra",
                                               ifelse(Score == 32, "Calidad Aceptable", "Calidad Excelente")))))
head(PC_BCD141_q)
str(PC_BCD141_q)
view(PC_BCD141_q)
PC_BCD141_q$Puntajeacolor <- 
  factor(PC_BCD141_q$Puntajeacolor, levels = c("Mala Calidad", "Baja Calidad", "Calidad Neutra", "Calidad Aceptable", "Calidad Excelente"))

factor(PC_BCD141_q$Puntajeacolor, levels = c("Mala Calidad", "Baja Calidad", "Calidad Neutra", "Calidad Aceptable", "Calidad Excelente"))

ggplot(PC_BCD141_q, aes(x = Cycle, y = Count/1000, fill = Puntajeacolor)) +
  geom_col(position = "stack") + 
  scale_fill_manual(values = c("Mala Calidad" = "darkred", "Baja Calidad" = "red", "Calidad Neutra" = "lightskyblue", "Calidad Aceptable" = "lightgreen", "Calidad Excelente"= "green"))+
  labs(y ="Conteo (miles)", x ="Ciclo (Posición Nt)", title = "Conteo de Calidad de nucleótidos por Ciclo en BC391", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = seq(0, 80, by = 5), label= seq(0, 80, by = 5)) +
  scale_y_continuous(breaks = seq(0, 700, by = 100), label= seq(0, 700, by = 100))#+
#coord_cartesian(ylim = c(0, 1000), xlim = c(1,12))
#gráfica de calidad (puntaje y conteo) por ciclo. 
#No hay mayor cantidad de secuencias de mala calidad al inicio, ni de baja calidad, solo casi ausencia de calidad excelente en los primeros ciclos. Considerar no usar el filtro o eliminar nadamás los nucleótidos de calidad baja, que le pasaría a los archivos? 
ggplot(PC_BCD141_bc, aes(x = Cycle, y = Count/1000, fill = BaseNitrogenada)) +
  geom_col(position = "stack") + 
  scale_fill_manual(values = c("A" = "black", "C" = "red", "G" = "lightskyblue", "T" = "lightgreen", "Calidad Excelente"= "green"))+
  labs(y ="Conteo (miles)", x ="Ciclo (Posición Nt)",title = "Conteo de Bases Nitrogenadas por Ciclo", subtitle = "Basado en Calidad de Phred") +
  theme(plot.title = element_text(size = 25),
        plot.subtitle = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, angle = 90),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))+
  scale_x_continuous(breaks = seq(0, 80, by = 5), label= seq(0, 80, by = 5)) +
  scale_y_continuous(breaks = seq(0, 700, by = 100), label= seq(0, 700, by = 100))#+
#bueno se puede ver que no hay un poli A claro en los primeros ciclos de las lecturas. 

q1 <- quality(arch1BCD14)
head(q1)
print(q1)

#Filtrar (explorativo)
#Resultados limpiar no nucleotidos solo quitó 20 lecturas, retirar por bases en general menor a 20 quitó 600K lecturas por lo tanto no se pueden retirar las lecturas por los pares de bases de mala calidad (solo retiré las de 14 PQ); 
#A hacer solo por avg quality. 
#srFilter(x, minQuality = 20) & avg quality (reads)

#código ejemplo: filtered_fq <- srFilter(fq, minQuality = 20, maxLength = 100, trim = TRUE, trimStart = 5, trimEnd = 5)
#filter_function <- function(sr) {
# Custom filtering logic
#avg_qual <- mean(as.numeric(quality(sr)))
#return(avg_qual >= 30)
#}
#filtered_fq <- srFilter(fq, filter = filter_function)

arch1BCD14 #no deja revisar desde aquí las calidades por readQualityScore solo en objetos qa.
tables(arch1BCD14)
#arch1assy0 <- quality()
arch1assy1 <- clean(arch1BCD14)
arch1assy2 <- as(quality(arch1assy1), "matrix")
quality_threshold <- 20
arch1assy3 <- arch1assy2 > quality_threshold
arch1assyNA <- is.na(arch1assy3)

## Count NAs per row
#na_per_row <- rowSums(is.na(arch1assy3))
#view(na_per_row)
#table(na_per_row)
## Count NAs per column
#na_per_column <- colSums(is.na(arch1assy3))
#
## Plotting NAs per row
#barplot(na_per_row, main = "NA Count Per Row", xlab = "Row", ylab = "Number of NAs", col = "lightblue")
#
## Plotting NAs per column
#barplot(na_per_column, main = "NA Count Per Column", xlab = "Column", ylab = "Number of NAs", col = "lightcoral")

table(arch1assy)
archy1assy3<- arch1assy[apply(arch1assy2,1, function(row) any(row, na.rm = T))]
arch1assy4 <- arch1BCD14[quality(arch1BCD14) ]


tables(arch1assy)
arch1assy0 <- 
filtros <- srFilterGenerator(filter = filter_byread_more30)

arch1assyfiltassay <- arch1BCD14[goodq(arch1BCD14)]  #regresa 152K lecturas de 760K. Es decir que casi todas las lecturas tienen alguna base menor a 20. y eso que son solo el 3% de las bases.


#goodq <- srFilter(
arch1BCD14[goodq(arch1BCD14)] 
  goodq <- function(x) {
  apply(as(quality(x), "matrix") > 20)
}
  
  
#filter_byread_more30 <- function(x){
#  apply(rowMeans(as(quality(x), "matrix"), 1, min, na.rm = TRUE))  > 30
#  #return(x[avg_qual >= 30])
#}

arch1BCD14 #primero limpiar y luego aplicar read_more_30_

table(rowMeans(as(quality(arch1BCD14), "matrix"), na.rm=T) >= 30)

#filter by average quality
read_more_30 <- function(x){
  x0<- rowMeans(as(quality(x), "matrix"), na.rm=T) >= 30
        return(x[x0])
}


#fun_filt1_goodbases <- srFilterGenerator(function(x){
#  apply(as(quality(x), "matrix"), 1, min, na.rm=TRUE) > 20}, .name= "Goodbases")
  

  #as(quality(arch1assy), "matrix")
ejemplo1 <- filter_reads_basemore20(arch1BCD14, 20)
str(ejemplo1) #se pierde mucho con bases arriba de 20. 
ejemplo4 <- funcion_lecturasmásde30ylimpieza(arch1BCD14)


#Sliding filter thresholds (20, 25, 30)
#ejemplo por IA para revisar. 
slidingWindowFilter <- function(shortreadq, window_size = 4, quality_threshold = 20) {
  # Extract quality scores from the ShortReadQ object
  quality <- as(quality(quality(shortreadq)), "matrix")
  
  # Calculate mean quality score per sliding window
  mean_quality <- apply(quality, 1, function(scores) {
    sapply(seq_along(scores), function(i) {
      # Take the average of the sliding window
      if (i + window_size - 1 <= length(scores)) {
        mean(scores[i:(i + window_size - 1)])
      } else {
        NA
      }
    })
  })
  
  # Find the position where the mean quality drops below the threshold
  cut_positions <- apply(mean_quality, 2, function(means) {
    cutoff <- which(means < quality_threshold)
    if (length(cutoff) == 0) {
      NA  # No cut needed
    } else {
      min(cutoff)  # First position below threshold
    }
  })
  
  # Trim reads at the calculated positions
  trimmed_reads <- trimTailw(shortreadq, cut_positions)
  return(trimmed_reads)
}

# Example usage
# Load a ShortReadQ object (replace 'path/to/file.fastq' with your file path)
fq <- readFastq("path/to/your/fastqfile.fastq")
filtered_fq <- slidingWindowFilter(fq, window_size = 4, quality_threshold = 20)

#tailw filter
#ejemplo por IA
# Function to apply tail filtering
tailFilter <- function(shortreadq, quality_threshold = 20, min_length = 36) {
  # Extract quality scores from the ShortReadQ object
  quality <- as(quality(quality(shortreadq)), "matrix")
  
  # Identify the position from the end where quality drops below the threshold
  cut_positions <- apply(quality, 1, function(scores) {
    below_threshold <- which(rev(scores) < quality_threshold)
    if (length(below_threshold) == 0) {
      NA  # No trimming needed
    } else {
      max(below_threshold)  # Farthest position below threshold from the end
    }
  })
  
  # Calculate trimming positions
  trim_positions <- ncol(quality) - cut_positions
  
  # Trim reads at the calculated positions ensuring minimum length
  trimmed_reads <- trimTailw(shortreadq, trim_positions, minlen = min_length)
  return(trimmed_reads)
}

# Example usage
tail_filtered_fq <- tailFilter(fq, quality_threshold = 20, min_length = 36)

#filter_reads_basemore20 <- function(fastq_obj, quality_threshold) {
#  # Extract the quality scores
#  quality_scores <- quality(fastq_obj)
#  
#  # Convert quality scores to a matrix
#  quality_matrix <- as(quality_scores, "matrix")
#  
#  # Define a logical matrix where TRUE means the score is above the threshold
#  filtered_matrix <- quality_matrix > quality_threshold
#  
#  # Keep reads with all bases above the quality threshold
#  # You can adjust the criteria as needed (e.g., any() for partial filtering)
#  filtered_reads <- fastq_obj[apply(filtered_matrix, 1, all, na.rm = TRUE)]
#  
#  return(filtered_reads)
#}
