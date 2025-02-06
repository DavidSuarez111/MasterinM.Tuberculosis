#El alineamiento falló y se optó por actualizar R
#install.packages("installr")
#library(installr)
#updateR()
#se corrió desde Rgui. 

# Librerías ---------------------------------------------------------------

R.Version()$version.string

#Para paquetes con conflictos
devtools::install_github("r-lib/conflicted")

#Cargar Paquetes
install.packages("BiocManager")
install.packages("ggplot2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("Rsamtools")
BiocManager::install("Rsubread")
BiocManager::install("ShortRead")
BiocManager::install("DESeq2")

#install.packages("Rsamtools")
#install.packages("Rsubread", update = T)
#install.packages("ShortRead")
install.packages("EnhancedVolcano")

library(Rsamtools)
library(Rsubread)
library(ShortRead)
library(DESeq2)
library(stringr)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(conflicted)
library(pheatmap)
#library(EnhancedVolcano)


# Directorios -------------------------------------------------------------


#Indexar un genoma Rsubread
getwd()
genome_path <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia"
list.dirs(genome_path)
list.files(genome_path)

BC391path <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BC 391/BC 391.fasta"
BC391index <- paste("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BC 391", "/BC391_index", sep = "")
BLpath <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BL 323/BL 323.fasta"
BLindex <- paste("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BL 323", "/BL_index", sep = "")
Musmusculuspath <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/Mus musculus/GCA_000001635.9_GRCm39_genomic.fna"
Musmusculusindex <- paste("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/Mus musculus", "/Mus_index", sep= "")

#AnnotationFiles
#musmusculus <- file?
list.files("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BC 391")
list.files("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BL 323")
BCGFF <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BC 391/BC 391.gff3"
BLGFF <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BL 323/BL 323.gff3"
#MusGFF <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia" #aún no lo tengo ok.
  
buildindex( #BC391
  basename = BC391index, 
  reference = BC391path)
buildindex( #BL
  basename = BLindex, 
  reference = BLpath)
buildindex(#Musmusculus
  basename = Musmusculusindex,
  reference = Musmusculuspath)

list.files(path = dirname(BC391index), pattern = basename(BC391index))
list.files(path = dirname(BLindex), pattern = basename(BLindex))
list.files(path = dirname(Musmusculusindex), pattern = basename(Musmusculusindex))

#Alinear los datos Rsamtools 
#rutaaarchivosBC391D14 <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14"
#rutaaarchivosBL323D14 <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14"

rutaaarchivosBC391D14 <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/BC d14"
rutaaarchivosBL323D14 <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/BL d14"

#BCD14_1:12filt
#BLD14_1:12filt

#rutaaarchivosBC391D14_2 <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BCD14_2"
#rutaaarchivosBL323D14_2 <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/Filtrados/BLD14_2"

#rutaaarchivosBC391D14_2 <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/BC d14 2" #cr"
#rutaaarchivosBL323D14_2 <- "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/Datos limpios fastq/BL d14 2 cr"
#no se ocupan los secundarios el bioproyecto dice que son de lectura simple no por pares. 

#BCD14_2_1:12
#BLD14_2_1:12

rutaoutputBAMBCD14 <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM"
rutaoutputBAMBLD14 <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM"

list.files(rutaoutputBAMBCD14)
list.files(rutaoutputBAMBLD14)
#Variables genomas indexados
BC391index
BLindex

# Alineamientos BAM -------------------------------------------------------

#Mapeo con genoma de ratón para corroborar que el alineamiento principal es a ratón. {naturaleza de librería RNAseq} 
paste(rutaaarchivosBC391D14, "/SRR BR1T1 BC D14.fastq", sep = "")
align( #1 -Mus musculus
  index = Musmusculusindex,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR1T1 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/MmD14_1.bam",
  nthreads = 4,
  input_format = "FASTQ"
) 
#96.7% de transcritos alineados a ratón. La naturaleza indica que no se separó RNA específico bacteriano. Contrario a lo expresado en 2021 por lisis diferencial revela transcriptoma de M. tb. o será que a pesar de eso? La verdad sus resultados eran inminentes entonces todo indica que es esto.

#Alineamientos Archivo 1 al 12. BC391. 
align( #1 BC
  index = BLindex,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR1T1 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_1_alineadoaBL.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #2 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR1T2 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_2.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #3 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR1T3 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_3.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #4 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR1T4 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_4.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #5 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR2T1 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_5.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #6 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR2T2 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_6.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #7 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR2T3 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_7.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #8 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR2T4 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_8.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #9 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR3T1 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_9.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #10 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR3T2 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_10.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #11 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR3T3 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_11.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #12 BC
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D14, "/SRR BR3T4 BC D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/OneDrive/Documentos/1. Maestría/Tesis Maestría/AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_12.bam",
  nthreads = 4,
  input_format = "FASTQ"
)

#Alineamientos Archivo 1 al 12. BL323
align( #1 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR1T1 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_1_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #2 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR1T2 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_2_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #3 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR1T3 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_3_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #4 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR1T4 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_4_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #5 -BL
  index = BLindex,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR2T1 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_5_alineadoaBL2.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #6 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR2T2 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_6_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #7 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR2T3 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_7_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #8 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR2T4 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_8_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #9 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR3T1 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_9_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #10 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR3T2 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_10_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #11 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR3T3 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_11_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #12 -BL
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D14, "/SRR BR3T4 BL D14.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_12_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)


#Resumen estadístico de alineamientos. 
ArchivosBC <- c("BC1", "BC2", "BC3", "BC4","BC5", "BC6", "BC7", "BC8", "BC9", "BC10", "BC11", "BC12")
MappedBC <- c(805, 728, 932, 799, 243, 241, 261, 248, 173, 155, 192, 188)
TotalreadsBC <- c(760323, 727193, 822406, 751707, 360094, 315519, 350527, 356507, 323216, 278112, 312284, 315266)
alineamientoDFBC <- data.frame(NombreArchivo = ArchivosBC, Mapeados = MappedBC, LecturasTotales = TotalreadsBC)
sumaTRBC <- sum(TotalreadsBC)
sumaMRBC <- sum(MappedBC)


ArchivosBL <- c("BL1", "BL2", "BL3", "BL4","BL5", "BL6", "BL7", "BL8", "BL9", "BL10", "BL11", "BL12")
MappedBL <- c(3995, 3711, 4633, 3992, 1218, 1041, 1258, 1167, 3230, 2857, 3297, 3271)
TotalreadsBL <- c(486048, 465565, 526771, 480462, 372425, 327610, 364503, 368597, 487787, 423531, 472250, 479957)
alineamientoDFBL <- data.frame(NombreArchivo = ArchivosBL, Mapeados = MappedBL, LecturasTotales = TotalreadsBL)
sumaTRBL <- sum(TotalreadsBL)
sumaMRBL <- sum(MappedBL)


AlineamientoBC1 <- alineamientoDFBC %>% mutate(Porcentaje = signif(((MappedBC/TotalreadsBC)*100), digits = 3))
AlineamientoBL1 <- alineamientoDFBL %>% mutate(Porcentaje = signif(((MappedBL/TotalreadsBL)*100), digits = 3))
alineamientoDFBL
promedioBC <- mean(AlineamientoBC1$Porcentaje) #0.078
promedioBL <- mean(AlineamientoBL1$Porcentaje) #0.61; concuerda con lo reportado en el artículo por pando y colaboradores 2022. "Of the total reads, 0.6% mapped to MTB genomes." Close Related Drug-Resistance Beijing...
promediolecturasBC <- mean(AlineamientoBC1$Mapeados[1:12])
promediolecturasBL <- mean(AlineamientoBL1$Mapeados[1:12])
promediolecturasBL/promediolecturasBC #(7 veces más y los primeros 4 archivos BC son más grandes)
resumenBCyBLmapeo <- rbind(alineamientoDFBC, alineamientoDFBL) #normalizando serían de 155 a 932 lecturas BC vs 1041 a 4633. 
promediolecturastotalesBC <- mean(AlineamientoBC1$LecturasTotales)
promediolecturastotalesBL <- mean(AlineamientoBL1$LecturasTotales)
resumenBCyBL_PorcMapeo <- resumenBCyBLmapeo %>% mutate(PorcentajedeMapeo = c(AlineamientoBC1$Porcentaje, AlineamientoBL1$Porcentaje), Resumen = c("Promedios de Alineamiento en Porcentaje", paste0(" BC: ",signif(round(promedioBC, digits = 2))), 
                                                                                                                                                    paste0(" BL: ",signif(round(promedioBL, digits = 2))),"",
                                                                                                                                                    "Rango de lecturas Mapeadas", "BC: 155 a 932", "BL: 1041 a 4633", "", "Promedio Lecturas", "BC: 414", "Promedio Lecturas", "BL: 2806", "7x librería", "",
                                                                                                                                                    "Promedio Lecturas Totales", "BC: 472763", "BL: 437959", "7% de diferencia en tamaño", rep("", times = 6)))
view(resumenBCyBL_PorcMapeo)

# Matriz de conteo --------------------------------------------------------
#Archivos en: 
rutaoutputBAMBCD14
list.files(rutaoutputBAMBCD14)

rutaoutputBAMBLD14
list.files(rutaoutputBAMBLD14)

#Ok tengo Bam files. Primero se ocupan indexar generando archivos .bai mediante samtools. Para indexar se ocupan ordenar los archivos BAM para tener las lecturas no mapeadas al final usando sortBam(.bam, destination = "")
#Ordenamiento de archivos .bam con sortBam()
#BCD14
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_1.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD1", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_2.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD2", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_3.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD3", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_4.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD4", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_5.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD5", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_6.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD6", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_7.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD7", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_8.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD8", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_9.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD9", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_10.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD10", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_11.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD11", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/BCD14_12.bam", destination = paste(rutaoutputBAMBCD14, "/sortBCD12", sep = ""))

#sort BLD14
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_1_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD1aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_2_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD2aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_3_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD3aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_4_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD4aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_5_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD5aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_6_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD6aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_7_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD7aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_8_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD8aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_9_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD9aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_10_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD10aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_11_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD11aBC", sep = ""))
sortBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/BLD14_12_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD14, "/sortBLD12aBC", sep = ""))

#Indexado de archivos. .baL a .bai. 
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD1.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD2.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD3.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD4.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD5.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD6.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD7.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD8.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD9.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD10.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD11.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD12.bam")

indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD1aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD2aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD3aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD4aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD5aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD6aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD7aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD8aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD9aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD10aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD11aBC.bam")
indexBam("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD12aBC.bam")
         
#ubisBCbai <- c("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD1.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD2.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD3.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD4.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD5.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD6.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD7.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD8.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD9.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD10.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD11.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD12.bam.bai") 
#ubisBLbai <- c("C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD1.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD2.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD3.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD4.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD5.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD6.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD7.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD8.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD9.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD10.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD11.bam.bai",
#               "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD12.bam.bai")
#no se ha cambiado los bais el nombre al alinear a BC pero creo que no se ocupa. 
str(ubisBCbai)
print(ubisBCbai[1])

BCGTF <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BC 391/BC 391.gff3"
BLGTF <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BL 323/BL 323.gff3"
#BLGTF no se usará. Para mantener la referencia en solo 1 genoma. 

#Tabla de conteos por featurecounts (siguiente)
#ensayo 1 tabla de conteo (funcionó) pal 14 si acaso (este no)
BC391_1_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD1.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
#llevar a todos los archivos BC391 y BL323 y ver seguimientos. 
#BC391
BC391_1_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD1.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 

BC391_2_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD2.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391_3_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD3.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391_4_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD4.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391_5_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD5.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391_6_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD6.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391_7_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD7.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391_8_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD8.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391_9_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD9.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391_10_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD10.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391_11_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD11.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391_12_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM/sortBCD12.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
#BL323
BL323_1_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD1aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323_2_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD2aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323_3_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD3aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323_4_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD4aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323_5_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD5aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323_6_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD6aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323_7_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD7aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323_8_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD8aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323_9_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD9aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGTF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323_10_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD10aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                             annot.ext = BCGTF,
                             isGTFAnnotationFile = TRUE,
                             GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                             GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                             useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                             nthreads = 4) 
BL323_11_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD11aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                             annot.ext = BCGTF,
                             isGTFAnnotationFile = TRUE,
                             GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                             GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                             useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                             nthreads = 4) 
BL323_12_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD12aBC.bam", 
                             annot.ext = BCGTF,
                             isGTFAnnotationFile = TRUE,
                             GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                             GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                             useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                             nthreads = 4) 
#listBC1fc <- c(BC391_1_tc, BC391_2_tc, BC391_3_tc, BC391_4_tc,
#               BC391_5_tc, BC391_6_tc, BC391_7_tc, BC391_8_tc,
#               BC391_9_tc, BC391_10_tc, BC391_11_tc, BC391_12_tc)
#
#listBL1fc <- c(BL323_1_tc, BL323_2_tc,  BL323_3_tc,  BL323_4_tc,
#               BL323_5_tc, BL323_6_tc,  BL323_7_tc,  BL323_8_tc,
#               BL323_9_tc, BL323_10_tc, BL323_11_tc, BL323_12_tc)

#realmente ocupo la subsección de counts
BC391_1_tc_c <- BC391_1_tc$counts
BC391_2_tc_c <-BC391_2_tc$counts 
BC391_3_tc_c <-BC391_3_tc$counts 
BC391_4_tc_c <-BC391_4_tc$counts
BC391_5_tc_c <-BC391_5_tc$counts
BC391_6_tc_c <-BC391_6_tc$counts 
BC391_7_tc_c <-BC391_7_tc$counts 
BC391_8_tc_c <-BC391_8_tc$counts
BC391_9_tc_c <-BC391_9_tc$counts 
BC391_10_tc_c<-BC391_10_tc$counts
BC391_11_tc_c<-BC391_11_tc$counts
BC391_12_tc_c<-BC391_12_tc$counts
  
BL323_1_tc_c <-BL323_1_tc$counts 
BL323_2_tc_c <-BL323_2_tc$counts
BL323_3_tc_c <-BL323_3_tc$counts
BL323_4_tc_c <-BL323_4_tc$counts
BL323_5_tc_c <-BL323_5_tc$counts
BL323_6_tc_c <-BL323_6_tc$counts
BL323_7_tc_c <-BL323_7_tc$counts
BL323_8_tc_c <-BL323_8_tc$counts
BL323_9_tc_c <-BL323_9_tc$counts
BL323_10_tc_c<-BL323_10_tc$counts
BL323_11_tc_c<-BL323_11_tc$counts
BL323_12_tc_c<-BL323_12_tc$counts


listBC1fc_c <- list(BC391_1_tc_c,  BC391_2_tc_c, BC391_3_tc_c,  BC391_4_tc_c,
                 BC391_5_tc_c,  BC391_6_tc_c, BC391_7_tc_c,  BC391_8_tc_c,
                 BC391_9_tc_c, BC391_10_tc_c,BC391_11_tc_c, BC391_12_tc_c)
str(listBC1fc_c)

listBL1fc_c <- list(BL323_1_tc_c,  BL323_2_tc_c,  BL323_3_tc_c,  BL323_4_tc_c,
                 BL323_5_tc_c,  BL323_6_tc_c,  BL323_7_tc_c,  BL323_8_tc_c,
                 BL323_9_tc_c, BL323_10_tc_c, BL323_11_tc_c, BL323_12_tc_c)


BC391_1_tc_c[,]
view(BC391_1_tc_c)
str(BC391_1_tc_c)
view(listBC1fc_c[2])
str(listBC1fc_c[1])

genesBC <- row.names(BC391_1_tc_c) #4002 (DEA)
genesBL <- row.names(BL323_1_tc_c) #3883 #no se usará para mantener misma referencia y poder normalizar y comparar.






procFC <- function(x){
  Rn <- row.names(x)
  y <- tibble(x)
  resdf <- as.data.frame(y) %>% mutate(genename = Rn) %>% select(genename, everything())
  return(resdf)
} #En esta función quiero integrar los genes como una variable, cambiar el nombre de los conteos a conteo y especie, filtrar genes no relevantes (será?), 
#para aplicar con todos los datos

procesBC <- sapply(listBC1fc_c, procFC)

dfBCd14 <- data.frame(
  genesBC = genesBC, 
  BC1T1 = procesBC[,1],
  BC1T2 = procesBC[,2],
  BC1T3 = procesBC[,3],
  BC1T4 = procesBC[,4],
  BC2T1 = procesBC[,5],
  BC2T2 = procesBC[,6],
  BC2T3 = procesBC[,7],
  BC2T4 = procesBC[,8],
  BC3T1 = procesBC[,9],
  BC3T2 = procesBC[,10],
  BC3T3 = procesBC[,11],
  BC3T4 = procesBC[,12]
)
#clipr::write_clip(colnames(dfBCd14))

dfBCD14_2 <- dfBCd14 %>% select(
  genesBC,
  BC1T1.sortBCD1.bam,
  BC1T2.sortBCD2.bam,
  BC1T3.sortBCD3.bam,
  BC1T4.sortBCD4.bam,
  BC2T1.sortBCD5.bam,
  BC2T2.sortBCD6.bam,
  BC2T3.sortBCD7.bam,
  BC2T4.sortBCD8.bam,
  BC3T1.sortBCD9.bam,
  BC3T2.sortBCD10.bam,
  BC3T3.sortBCD11.bam,
  BC3T4.sortBCD12.bam
)

dfBC_c <- dfBCD14_2 %>% rename(
  Conteo_BC1T1 = BC1T1.sortBCD1.bam,
  Conteo_BC1T2 = BC1T2.sortBCD2.bam,
  Conteo_BC1T3 = BC1T3.sortBCD3.bam,
  Conteo_BC1T4 = BC1T4.sortBCD4.bam,
  Conteo_BC2T1 = BC2T1.sortBCD5.bam,
  Conteo_BC2T2 = BC2T2.sortBCD6.bam,
  Conteo_BC2T3 = BC2T3.sortBCD7.bam,
  Conteo_BC2T4 = BC2T4.sortBCD8.bam,
  Conteo_BC3T1 = BC3T1.sortBCD9.bam,
  Conteo_BC3T2 = BC3T2.sortBCD10.bam,
  Conteo_BC3T3 = BC3T3.sortBCD11.bam,
  Conteo_BC3T4 = BC3T4.sortBCD12.bam
)
view(dfBC_c)

procesBL <- sapply(listBL1fc_c, procFC)

dfBLd14 <- data.frame(
  genesBC2 = genesBC, #no se usaron los de BL para normalizar y comparar posteriormente. 
  BL1T1 = procesBL[,1],
  BL1T2 = procesBL[,2],
  BL1T3 = procesBL[,3],
  BL1T4 = procesBL[,4],
  BL2T1 = procesBL[,5],
  BL2T2 = procesBL[,6],
  BL2T3 = procesBL[,7],
  BL2T4 = procesBL[,8],
  BL3T1 = procesBL[,9],
  BL3T2 = procesBL[,10],
  BL3T3 = procesBL[,11],
  BL3T4 = procesBL[,12]
)
view(dfBLd14)

dfBLD14_2 <- dfBLd14 %>% select(
  genesBC2,
  BL1T1.sortBLD1aBC.bam,
  BL1T2.sortBLD2aBC.bam,
  BL1T3.sortBLD3aBC.bam,
  BL1T4.sortBLD4aBC.bam,
  BL2T1.sortBLD5aBC.bam,
  BL2T2.sortBLD6aBC.bam,
  BL2T3.sortBLD7aBC.bam,
  BL2T4.sortBLD8aBC.bam,
  BL3T1.sortBLD9aBC.bam,
  BL3T2.sortBLD10aBC.bam,
  BL3T3.sortBLD11aBC.bam,
  BL3T4.sortBLD12aBC.bam
)
view(dfBLD14_2)
#clipr::write_clip(colnames(dfBLd14))

dfBL_c <- dfBLD14_2 %>% rename(
  Conteo_BL1T1 = BL1T1.sortBLD1aBC.bam,
  Conteo_BL1T2 = BL1T2.sortBLD2aBC.bam,
  Conteo_BL1T3 = BL1T3.sortBLD3aBC.bam,
  Conteo_BL1T4 = BL1T4.sortBLD4aBC.bam,
  Conteo_BL2T1 = BL2T1.sortBLD5aBC.bam,
  Conteo_BL2T2 = BL2T2.sortBLD6aBC.bam,
  Conteo_BL2T3 = BL2T3.sortBLD7aBC.bam,
  Conteo_BL2T4 = BL2T4.sortBLD8aBC.bam,
  Conteo_BL3T1 = BL3T1.sortBLD9aBC.bam,
  Conteo_BL3T2 = BL3T2.sortBLD10aBC.bam,
  Conteo_BL3T3 = BL3T3.sortBLD11aBC.bam,
  Conteo_BL3T4 = BL3T4.sortBLD12aBC.bam
)

#Resultados tablas de gen y conteo por cepa: 
view(dfBC_c)
view(dfBL_c)
dfBL

#filtrado de los data frames: 
# (not used) filtered_counts <- fc_result$counts[rowSums(fc_result$counts >= 50) > 0, ],
#better approach: df$TotalCounts <- rowSums(df[, -1]); filtered_df <- df[df$TotalCounts >= 50, ]

#!!no filtrar no correr las siguientes lineas porque afectan la normalización, añaden 1 columna innecesaria. 
#dfBC_c$conteofilas <- rowSums(dfBC_c[,-1])
#dfBC_filtrado <- dfBC_c[dfBC_c$conteofilas >= 1,]
#dfBC_filtrado <- dfBC_c[dfBC_c$conteofilas >= ?]
view(dfBC_filtrado) #de (4002 genes a solo 466 expresados) posiblemente suficiente filtro por naturaleza de datos.

#dfBL_c$conteofilas <- rowSums(dfBL_c[,-1])
#dfBL_filtrado <- dfBL_c[dfBL_c$conteofilas >= 1,] 
view(dfBL_filtrado) #de 3883 a 1368 filas. No usaré diferente filtro me parece injusto, mal diseño.
#usando BC fue de 4002 a 1412 filas. 

#Resultados de conteos filtrados
#dfBC_filtrado
#dfBL_filtrado

#Unir los marcos de datos en 1 tabla. 
dfBCyBL <- cbind(dfBC_c, dfBL_c)
colnames(dfBCyBL)
dfBCyBLconteoslimpios <- dfBCyBL %>% select(-genesBC2)
view(dfBCyBLconteoslimpios) #4002 filas
dfBCyBLconteoslimpios$conteofilas <- rowSums(dfBCyBLconteoslimpios[,-1])
dfBCyBLfiltrado <- dfBCyBLconteoslimpios[dfBCyBLconteoslimpios$conteofilas >= 1,]
view(dfBCyBLfiltrado) #1596 filas. 
dfBCyBLfilt <- dfBCyBLfiltrado %>% select(-conteofilas)
muestras = colnames(dfBCyBLfilt[,-1])
view(dfBCyBLfilt)
#mas estricto el filtrado? 

#Resultado: 
dfBCyBLfilt

#. Generar un marco de datos de metadatos. 
metadatosBCyBL <- data.frame(
  IDdemuestras = muestras,
  cepa = rep(c("BC391", "BL323"), each = 12),
  replicadoB = rep(1:3, each = 4, times = 2),
  replicadoT = rep(1:4, times = 6)
)

view(metadatosBCyBL)
#REsultado: 
metadatosBCyBL #es el diseño del experimento, que condiciones hay. Que varía entre cada muestra. 
#RAW PCA 
dfBCyBLfilt
metadatosBCyBL


# Normalización y organización   ------------------------------------------
#1. Entre replicados técnicos se pueden emplear conteos crudos en el PCA. Solo distingue cercanía entre RT.
#2. Se emplearán 2 métodos de normalización idóneos/compatibles con expresión diferencial por las librerías: DESeq sizefactors y EdgeR TMM. 
#3. El diseño será de 2x2. 2 métodos de normalización y 2 cantidades de muestras (12 y 24: para análisis intra especie e interespecie respectivamente)


#edgeR (TMM). útil para diferencias dentro de lo normalmente expresado no para outliers. 
#código ejemplo primario: #edgeR (calcNormFactors corrige entre tamaños de librerías)
# Create DGEList object: y <- DGEList(counts = count_matrix, group = sample_info$condition)
# Calculate normalization factors: y <- calcNormFactors(y)
# Access normalized counts: normalized_counts <- cpm(y, normalized.lib.sizes = TRUE)

    #TMM; 12 y 24. (v1)

#PCA (EdgeR)

#DESeq2 (sizeFactors)
#Código ejemplo primario: #DESeq (sizeFactors automáticamente corrige entre diferencias en profundidad de secuenciado)
#dds <- DESeqDataSetFromMatrix(countData = count_matrix,
#                              colData = sample_info,  # Sample information
#                              design = ~ condition)   # Experimental design
# Run DESeq to normalize: dds <- DESeq(dds)
# Access normalized counts: normalized_counts <- counts(dds, normalized = TRUE)
   

    #sizeFactors; 12 y 24. (v2)
#12 y 12 (intracepa)
BCdf <- dfBCyBLfilt[,2:13]
row.names(BCdf) <- dfBCyBLfilt[,1]
BCdf$conteofilas <- rowSums(BCdf)
BCdffilt <- BCdf[BCdf$conteofilas >= 1,]
BCdffilt1 <- BCdffilt %>% select(-conteofilas)
muestrasBC <- colnames(BCdf[,2:13])
metadatosBC <- data.frame(
  IDdemuestras = muestrasBC,
  replicadoB = rep(1:3, each = 4),
  replicadoT = rep(1:4, times = 3)
)
metadatosBC$replicadoB <- as.factor(metadatosBC$replicadoB)
metadatosBC$replicadoT <- as.factor(metadatosBC$replicadoT)
ddsBC <- DESeqDataSetFromMatrix(countData = BCdffilt1,
                                   colData = metadatosBC,
                                   design = ~ replicadoB + replicadoT)
ddsBCnormalizado <- DESeq(ddsBC, fitType = "mean")
str(ddsBCnormalizado)
view(ddsBCnormalizado)

conteos_normalizadosBC <- counts(ddsBCnormalizado, normalized = T)  
str(conteos_normalizadosBC)
view(conteos_normalizadosBC)
BCvst <- varianceStabilizingTransformation(ddsBCnormalizado, blind = F) #vst mantiene linear los conteos bajos y aplica un logaritmo a los conteos altos. Para mantener los números en la misma escala. 
accederanormalizadosBC <- assay(BCvst)
view(accederanormalizadosBC)
BC_PCA_RbRT <- plotPCA(BCvst, intgroup = "replicadoB", returnData = T)
str(BC_PCA_RbRT)
plotPCA(BCvst, intgroup = "replicadoB") 
pca_ggplotBC <- ggplot(BC_PCA_RbRT, aes(x = PC1, y = PC2, color = replicadoB)) +
  geom_point(size = 3) +  # Adjust size as needed
  labs(title = "PCA de expresión genética normalizada de archivos BC", subtitle ="Agrupamiento por replicados Biológicos") +
  theme_minimal()
print(pca_ggplotBC)

BLdf <- dfBCyBLfilt[,14:25]
row.names(BLdf) <- dfBCyBLfilt[,1]
muestrasBL <- colnames(BLdf[,1:12])
metadatosBL <- data.frame(
  IDdemuestras = muestrasBL,
  replicadoB = rep(1:3, each = 4),
  replicadoT = rep(1:4, times = 3)
)
metadatosBL$replicadoB <- as.factor(metadatosBL$replicadoB)
metadatosBL$replicadoT <- as.factor(metadatosBL$replicadoT)
ddsBL <- DESeqDataSetFromMatrix(countData = BLdf,
                                colData = metadatosBL,
                               design = ~ replicadoB + replicadoT)
ddsBLnormalizado <- DESeq(ddsBL)
conteos_normalizadosBL <- counts(ddsBLnormalizado, normalized = T)
view(conteos_normalizadosBL)
BLvst <- varianceStabilizingTransformation(ddsBLnormalizado, blind = F)
view(BLvst)
BLver <- assay(BLvst)
view(BLver)
BL_PCA_RbRT <- plotPCA(BLvst, intgroup = "replicadoB", returnData = T)
str(BL_PCA_RbRT)
plotPCA(BLvst, intgroup = "replicadoB") #grupo 2 falla, considerar remover? 
pca_ggplotBL <- ggplot(BL_PCA_RbRT, aes(x = PC1, y = PC2, color = replicadoB)) +
  geom_point(size = 3) +  # Adjust size as needed
  labs(title = "PCA de expresión genética normalizada de archivos BL", subtitle ="Agrupamiento por replicados Biológicos") +
  theme_minimal()
print(pca_ggplotBL)

#24
view(dfBCyBLfilt)
metadatosBCyBL #gen más expresado gene- "BJM02_06990" (ribosomal 23 y 16s)

#matriz de conteo sin nombre de genes como variable o columna. 
conteogenes <- dfBCyBLfilt[,-1]
rownames(conteogenes) <- dfBCyBLfilt[,1]
view(conteogenes)

ddsBCyBL <- DESeqDataSetFromMatrix(countData = conteogenes,
                                   colData = metadatosBCyBL,
                                   design = ~ cepa + replicadoB)
str(ddsBCyBL)
ddsBCyBL2 <- DESeqDataSetFromMatrix(countData = conteogenes,
                                   colData = metadatosBCyBL,
                                   design = ~ replicadoB + cepa)
ddsBCyBL3 <- DESeqDataSetFromMatrix(countData = conteogenes,
                                   colData = metadatosBCyBL,
                                   design = ~ replicadoT + cepa)
NormalizadoBCyBL <- DESeq(ddsBCyBL)
head(counts(NormalizadoBCyBL, normalized = T))

NormalizadoBCyBL2 <- DESeq(ddsBCyBL2)
NormalizadoBCyBL3 <- DESeq(ddsBCyBL3)

conteos_normalizadosBCyBL <- counts(NormalizadoBCyBL, normalized = T)
view(conteos_normalizadosBCyBL)
conteos_normalizadosBCyBL2 <- counts(NormalizadoBCyBL2, normalized = T)
conteos_normalizadosBCyBL3 <- counts(NormalizadoBCyBL3, normalized = T)
head(conteos_normalizadosBCyBL)

BCyBLvst <- varianceStabilizingTransformation(NormalizadoBCyBL, blind = F)
head(BCyBLvst)
BCyBLvst2 <- varianceStabilizingTransformation(NormalizadoBCyBL2, blind = F)
BCyBLvst3 <- varianceStabilizingTransformation(NormalizadoBCyBL3, blind = F)
#PCA (DESeq)

PCA1cRbRT <- plotPCA(BCyBLvst, intgroup = c("cepa", "replicadoB", "replicadoT"), returnData = T) #ambas cepas

PCAggcepa <- plotPCA(BCyBLvst, intgroup = "cepa", returnData = T) #se puede ver una amplia dispersión intracepa en BC391, los datos por BL323 están agrupados. #No se está considerando el tamaño de la librería original sino de la librería secundaria las lecturas totales unidas a Mtv para normalizar. La naturaleza de los datos afecta la normalización. La cantidad de lecturas mapeadas a tb por lecturas totales suguieren una mayor expresión de RNA por BL323 al día 14.  
plotPCA(BCyBLvst, intgroup = "cepa") #se puede ver una amplia dispersión intracepa en BC391, los datos por BL323 están agrupados. #No se está considerando el tamaño de la librería original sino de la librería secundaria las lecturas totales unidas a Mtv para normalizar. La naturaleza de los datos afecta la normalización. La cantidad de lecturas mapeadas a tb por lecturas totales suguieren una mayor expresión de RNA por BL323 al día 14.  
plotPCA(BCyBLvst, intgroup = "replicadoB")
plotPCA(BCyBLvst, intgroup = c("cepa", "replicadoB"))
plotPCA(BCyBLvst, intgroup = c("cepa", "replicadoB", "replicadoT"))
        
PCA1cRbRT2rBc <- plotPCA(BCyBLvst2, intgroup = c("replicadoB", "cepa"), returnData = T) #ambas cepas
plotPCA(BCyBLvst2, intgroup = c("replicadoB", "cepa"))#, returnData = T) #ambas cepas
PCA1cRbRT3 <- plotPCA(BCyBLvst3, intgroup = c("replicadoT"), returnData = T) #ambas cepas
PCA1cRbRT3 <- plotPCA(BCyBLvst3, intgroup = c("replicadoT"))#, returnData = T) #ambas cepas


#str(PCA1cRbRT)
#str(PCA1forgg)
#head(PCA1cRbRT)
#plotPCA(BCyBLvst, intgroup = c("ReplicadoB", "ReplicadoT")) #ambas cepas
#plotPCA(BCyBLvst, intgroup = "ReplicadoT") #ambas cepas

#pasar a ggplot
PCA1forgg <- PCA1cRbRT$data
PCA1forgg$CombinedGroup <- interaction(PCA1forgg$cepa, PCA1forgg$replicadoB, PCA1forgg$replicadoT)
PCA1forgg$replicadoB <- as.factor(PCA1forgg$replicadoB)
PCA1forgg$cepa <- as.factor(PCA1forgg$cepa)
PCA1forgg$replicadoT <- as.factor(PCA1forgg$replicadoT)
pca_ggplot <- ggplot(PCA1forgg, aes(x = PC1, y = PC2, color = cepa, shape = replicadoB)) +
  geom_point(size = 3) +  # Adjust size as needed
  labs(title = "PCA de archivos BC y BL", subtitle = "Color Por Cepa, Forma por Replicado Biológico") +
  theme_minimal()
pca_ggplot0 <- ggplot(PCAggcepa, aes(x = PC1, y = PC2, color = cepa)) +
  geom_point(size = 3) +  # Adjust size as needed
  labs(title = "PCA de archivos BC y BL", subtitle = "Color Por Cepa") +
  theme_minimal()
pca_ggplot2 <- ggplot(PCA1forgg, aes(x = PC1, y = PC2, color = cepa, shape = replicadoT)) +
  geom_point(size = 3) + 
  labs(title = "PCA: Color by Strain, Shape by Technical Replicate") +
  theme_minimal()

#pca_ggplot3_combinado <- ggplot(PCA1forgg, aes(x = PC1, y = PC2, color = CombinedGroup)) +
#  geom_point(size = 3) +
#  labs(title = "PCA with Combined Strain, Biological and Technical Replicates") +
#  theme_minimal() #muy poco agradable de ver. 

print(pca_ggplot0)
print(pca_ggplot)
print(pca_ggplot2)
print(pca_ggplot3_combinado)

#Análisis Transcriptómico. 
#Input objetos DESeq

view(NormalizadoBCyBL)
NormalizadoBCyBL
ddsBCnormalizado
ddsBLnormalizado

#colapsar los replicados técnicos. 
#función collapseReplicates(dds, groupby = dds$bio_rep) #utiliza DESeq() cool.
#design(dds_collapsed) <- ~ strain
#dds_collapsed <- estimateSizeFactors(dds_collapsed)
#sizeFactors(dds_collapsed)
technical_replicates <- c("muestra1_tech1", "muestra1_tech2", "muestra1_tech3", "muestra1_tech4",
                          "muestra2_tech1", "muestra2_tech2", "muestra2_tech3", "muestra2_tech4",
                          "muestra3_tech1", "muestra3_tech2", "muestra3_tech3", "muestra3_tech4",
                          "muestra4_tech1", "muestra4_tech2", "muestra4_tech3", "muestra4_tech4",
                          "muestra5_tech1", "muestra5_tech2", "muestra5_tech3", "muestra5_tech4",
                          "muestra6_tech1", "muestra6_tech2", "muestra6_tech3", "muestra6_tech4")
biological_replicates <- c("muestra1", "muestra1", "muestra1", "muestra1",
                           "muestra2", "muestra2", "muestra2", "muestra2",
                           "muestra3", "muestra3", "muestra3", "muestra3",
                           "muestra4", "muestra4", "muestra4", "muestra4",
                           "muestra5", "muestra5", "muestra5", "muestra5",
                           "muestra6", "muestra6", "muestra6", "muestra6")
strain <- c("BC391", "BC391", "BC391", "BC391",
            "BC391", "BC391", "BC391", "BC391",
            "BC391", "BC391", "BC391", "BC391", 
            "BL323", "BL323", "BL323", "BL323", 
            "BL323", "BL323", "BL323", "BL323", 
            "BL323", "BL323", "BL323", "BL323")


ddsBCyBLcolapsadoRT <- collapseReplicates(NormalizadoBCyBL, groupby = biological_replicates, run = technical_replicates)
dim(ddsBCyBLcolapsadoRT)

ddsBCyBLcolapsadoRT <- DESeq(ddsBCyBLcolapsadoRT)
dim(ddsBCyBLcolapsadoRT)
head(ddsBCyBLcolapsadoRT)


resBCyBLcol <- results(ddsBCyBLcolapsadoRT, contrast = c("cepa", "BL323", "BC391"))
view(resBCyBLcol)

#res_pval <- resBCyBLcol[which(resBCyBLcol$pvalue > 0.95), ] #estudiar sobre p value y adjusted. 
#res_signif <- resBCyBLcol[which(resBCyBLcol$padj > 0.95), ]
res_pval <- resBCyBLcol[which(resBCyBLcol$pvalue < 0.10), ] #estudiar sobre p value y adjusted. 
res_signif <- resBCyBLcol[which(resBCyBLcol$padj < 0.10), ]
view(res_pval)
view(res_signif)

plotMA(resBCyBLcol, main = "MA-plot", ylim = c(-2,2))
colData(vsdBC)
#pasar a ggplot el PCA y cambiar a color por cepa. 
vsdBC <- varianceStabilizingTransformation(ddsBCyBLcolapsadoRT, blind=FALSE)
plotPCA(vsdBC, intgroup= "cepa")#c("cepa", "replicadoB")) #cambiar a ggplot? 
PCAgg <- plotPCA(vsdBC, intgroup= "cepa", returnData = T)#c("cepa", "replicadoB")) #cambiar a ggplot? 

ggplot(PCAgg, aes(x = PC1, y = PC2, color = cepa)) +
  geom_point(size = 5) +  # Adjust size as needed
  labs(title = "PCA de RB de BC y BL al día 14", subtitle = "Color Por Cepa", y = "PC2 28% varianza", x = "PC1 39% varianza") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

top_genes <- head(order(rowVars(assay(vsdBC)), decreasing=TRUE), 50)
mat <- assay(vsdBC)[top_genes, ]
mat <- mat - rowMeans(mat)  # Center rows
str(mat)
dim(mat)
colnames(mat)
colData(mat)
pheatmap(mat) #normalización juega en contra. BC es la que aparece sobreexpresando. Interesante revisar a detalle esas expresiones de todas formas. 

#volcan 
resBCyBLcol
resBCyBLv <- na.omit(resBCyBLcol) #seleccionar solo los 35 primeros o mínimo un zoom sobre eso. 
view(resBCyBLv)
resBCyBLv$significance <- ifelse(resBCyBLv$pvalue < 0.05, "Significant", "Not significant")
resBCyBLv$significancepval <- ifelse(resBCyBLv$significance == "Significant", "Yes", "No" )
#resBCyBLv$ifpval <- 
top_genes <- head(resBCyBLv[order(resBCyBLv$log2FoldChange), ], 10)  # Adjust number of top genes as needed
view(top_genes)
top_genes_BL_pval_0.05 <- top_genes[top_genes$pvalue<0.05]
top_genesmore <- head(resBCyBLv[order(resBCyBLv$log2FoldChange, decreasing = T), ], 10)  # Adjust number of top genes as needed
view(top_genesmore)
top_genes_BC_pval_0.05 <- top_genesmore[top_genesmore$pvalue<0.05]
view(top_genes_BL_pval_0.05)
view(top_genes_BC_pval_0.05) #No hay en BC

view(resBCyBLv)

ggplot(resBCyBLv, aes(x=log2FoldChange, y=pvalue, color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genes, aes(label=row.names(top_genes)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesmore, aes(label=row.names(top_genesmore)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  #geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue")+
  #geom_vline(xintercept = c(-0.807, -1.807, -2.807, -3.807, -4.807), linetype = "dashed", color = "blue")+
  labs(title="Expresión diferencial BL323 vs BC391", x="Log2 Fold Change", y="P-value", subtitle = "Top 10 up/down") +
  theme(legend.position="")

ggplot(resBCyBLv, aes(x=log2FoldChange, y=pvalue, color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genes, aes(label=row.names(top_genes)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesmore, aes(label=row.names(top_genesmore)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")+
  geom_vline(xintercept = c(-0.807, -1.807, -2.807, -3.807, -4.807), linetype = "dashed", color = "blue")+
  annotate("text", x = -2.807, y = max(resBCyBLv$pvalue, na.rm = TRUE), 
           label = "0", color = "blue", hjust = 1, angle = 0)+
  annotate("text", x = -3.807, y = max(resBCyBLv$pvalue, na.rm = TRUE), 
                                                                         label = "-1", color = "blue", hjust = 1, angle = 0)+
  
  annotate("text", x = -4.807, y = max(resBCyBLv$pvalue, na.rm = TRUE), 
                                                                         label = "-2", color = "blue", hjust = 1, angle = 0)+
  
  annotate("text", x = -1.807, y = max(resBCyBLv$pvalue, na.rm = TRUE), 
                                                                         label = "1", color = "blue", hjust = 1, angle = 0)+
  
  annotate("text", x = -0.807, y = max(resBCyBLv$pvalue, na.rm = TRUE), 
                                                                         label = "2", color = "blue", hjust = 1, angle = 0)+
  
  labs(title="Expresión diferencial BL323 vs BC391", x="Log2 Fold Change", y="P-value", subtitle = "Conversión log2(7) = 2.807 = 0") +
  theme(legend.position="") #right

ggplot(resBCyBLv, aes(x=log2FoldChange, y=-log10(pvalue), color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genes, aes(label=row.names(top_genes)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesmore, aes(label=row.names(top_genesmore)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")+
  geom_vline(xintercept = c( -0.807, -1.807, -2.807, -3.807, -4.807), linetype = "dashed", color = "blue")+ #1.193, 0.193,
  annotate("text", x = -2.807, y = 1.5, 
           label = "0", color = "blue", hjust = 1, angle = 0)+
  annotate("text", x = -3.807, y = 1.5, 
           label = "-1", color = "blue", hjust = 1, angle = 0)+
  
  annotate("text", x = -4.807, y = 1.5, 
           label = "-2", color = "blue", hjust = 1, angle = 0)+
  
  annotate("text", x = -1.807, y = 1.5, 
           label = "1", color = "blue", hjust = 1, angle = 0)+
  
  annotate("text", x = -0.807, y = 1.5, 
           label = "2", color = "blue", hjust = 1, angle = 0)+
  
  labs(title="Expresión diferencial BL323 vs BC391", x="Log2 Fold Change", y="-Log10 P-value", subtitle = "Gráfico de Volcán") +
  theme(legend.position="") #right
  #annotate("text", x = .193, y = 1.5, 
  #         label = "3", color = "blue", hjust = 1, angle = 0)+
  
  #annotate("text", x = 1.193, y = 1.5, 
  #         label = "4", color = "blue", hjust = 1, angle = 0)+
  

#de 8 a 30 veces más expresado que en BC. los top 10. 
#genes de 2 a 8 veces más expresados en BC que en BL. 
top_genesBLmas <- head(resBCyBLcol[order(resBCyBLcol$log2FoldChange), ], 20)  
top_genesBCmas <- head(resBCyBLcol[order(resBCyBLcol$log2FoldChange, decreasing = T), ], 20)  
view(top_genesBCmas)
view(top_genesBLmas)
dfExpresiónmasenBLqueBC <- data.frame(top_genesBLmas)
dfExpresiónmasenBCqueBL <- data.frame(top_genesBCmas)
view(dfExpresiónmasenBLqueBC)
view(dfExpresiónmasenBCqueBL)

#Anotaciones