
# Réplica al día 28 BC y BL hasta alineamientos ---------------------------------------------------
#librerías
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

#Librerías de archivos 
BLd28 <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/BLd28"
BCd28 <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/BCd28"
list.files(BLd28) #nombre ejemplo "BLD28 RB1T1.fastq"
list.files(BCd28)

# Alineamiento ------------------------------------------------------------
#Objetivo: usar el código ejemplo y editarlo para funcionalidad. Adaptarlo a generar resultados para BLD14 y BLD28 alineados a referencia BL y analísis entre ellos. 
genome_path <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia"
BC391pathgenome <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BC 391/BC 391.fasta"
BC391index <- paste("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BC 391", "/BC391_index", sep = "")
#BLpath <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BL 323/BL 323.fasta"
#BLindex <- paste("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BL 323", "/BL_index", sep = "")
Musmusculuspath <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/Mus musculus/GCA_000001635.9_GRCm39_genomic.fna"
Musmusculusindex <- paste("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/Mus musculus", "/Mus_index", sep= "")

BCGFF <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BC 391/BC 391.gff3"
#BLGFF <- "C:/Users/david/Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BL 323/BL 323.gff3"

rutaaarchivosBC391D14 <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/BC d14"
rutaaarchivosBC391D28 <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/BCd28"
rutaaarchivosBL323D14 <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/BL d14"
rutaaarchivosBL323D28 <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/BLd28"

rutaoutputBAMBCD14 <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD14 BAM"
rutaoutputBAMBCD28 <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM"
rutaoutputBAMBLD14 <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM"
rutaoutputBAMBLD28 <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM" 

paste(rutaaarchivosBC391D14, "/SRR BR1T1 BC D14.fastq", sep = "")

align( #1 -Mus musculus #Dato ejemplo BCD28
  index = Musmusculusindex,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB1T1.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/MmBCD28_1.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -Mus musculus #Dato ejemplo BLD28
  index = Musmusculusindex,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB1T1.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/MmBLD28_1.bam",
  nthreads = 4,
  input_format = "FASTQ"
)

#Alineamiento de BCD28 a BC index 
#D14
align( #1 -BCD28 1
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB1T1.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_1_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
#D28
align( #1 -BCD28 2
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB1T2.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_2_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BCD28 3
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB1T3.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_3_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BCD28 4
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB1T4.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_4_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BCD28 5
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB2T1.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_5_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BCD28 6
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB2T2.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_6_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BCD28 7
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB2T3.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_7_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BCD28 8
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB2T4.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_8_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BCD28 9
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB3T1.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_9_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BCD28 10
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB3T2.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_10_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BCD28 11
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB3T3.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_11_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BCD28 12
  index = BC391index,
  readfile1 = paste(rutaaarchivosBC391D28, "/BCD28 RB3T4.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_12_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)

#BLD28 a BC. 
align( #1 -BLD28 1
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB1T1.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_1_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 2
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB1T2.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_2_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 3
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB1T3.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_3_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 4
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB1T4.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_4_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 5
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB2T1.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_5_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 6
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB2T2.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_6_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 7
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB2T3.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_7_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 8
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB2T4.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_8_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 9
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB3T1.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_9_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 10
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB3T2.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_10_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 11
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB3T3.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_11_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)
align( #1 -BLD28 12
  index = BC391index,
  readfile1 = paste(rutaaarchivosBL323D28, "/BLD28 RB3T4.fastq", sep = ""),
  #readfile2 = paste(rutaaarchivosBC391D14_2, "/2 SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_12_alineadoaBC.bam",
  nthreads = 4,
  input_format = "FASTQ"
)

#Resumen estadístico de Alineamiento. #comparable estadísticamente? 
ArchivosBCd28 <- c("BCd28_1", "BCd28_2", "BC28_3", "BC28_4",
                  "BC28_5", "BC28_6", "BC28_7", "BC28_8",
                  "BC28_9", "BC28_10", "BC28_11", "BC28_12")
MappedBCd28 <- c(1019, 926, 1122, 1041, 1518, 1415, 1695, 1465, 2212, 2032, 2408, 2216)
TotalReadsBCd28 <- c(316888, 304443, 343236, 315347, 482381, 458231, 525108, 479213, 406085, 386310, 442414, 398867)
alineamientoDFBCD28 <-data.frame(NombreArchivo = ArchivosBCd28, Mapeados = MappedBCd28, LecturasTotales = TotalReadsBCd28)
sumaTRBCd28 <- sum(TotalReadsBCd28)
sumaMRBCd28 <- sum(MappedBCd28)

ArchivosBLd28 <- c("BLd28_1", "BLd28_2", "BLd28_3", "BLd28_4",
                   "BLd28_5", "BLd28_6", "BLd28_7", "BLd28_8",
                   "BLd28_9", "BLd28_10", "BLd28_11", "BLd28_12")
MappedBLd28 <- c(2217, 1860, 2220, 2210, 4337, 3773, 4292, 4146, 4329, 3677, 4350, 4248)
TotalReadsBLd28 <- c(430480, 377710, 421417, 426017, 504917, 435718, 483221, 491119, 563673, 479075, 533634, 544014)
alineamientoDFBLD28 <-data.frame(NombreArchivo = ArchivosBLd28, Mapeados = MappedBLd28, LecturasTotales = TotalReadsBLd28)
sumaTRBLd28 <- sum(TotalReadsBLd28)
sumaMRBLd28 <- sum(MappedBLd28)

#porcentajes
AlineamientoBCD28 <- alineamientoDFBCD28 %>% mutate(Porcentaje = signif(((MappedBCd28/TotalReadsBCd28)*100), digits = 3))
AlineamientoBLD28 <- alineamientoDFBLD28 %>% mutate(Porcentaje = signif(((MappedBLd28/TotalReadsBLd28)*100), digits = 3))

promedioBCD28 <- mean(AlineamientoBCD28$Porcentaje) #0.39 (menos de x2)
promedioBLD28 <- mean(AlineamientoBLD28$Porcentaje) #0.72

#unir ambos
promediolecturasBCD28 <- mean(AlineamientoBCD28$Mapeados[1:12])
promediolecturasBLD28 <- mean(AlineamientoBLD28$Mapeados[1:12])
promediolecturasBLD28/promediolecturasBCD28 #(exactamente 1? raro)
promediolecturasBCD28/promediolecturasBLD28 #(exactamente 1? raro)
resumenBCyBLmapeoD28 <- rbind(alineamientoDFBCD28, alineamientoDFBLD28) #normalizando serían de 155 a 932 lecturas BC vs 1041 a 4633. 
promediolecturastotalesBCD28 <- mean(AlineamientoBCD28$LecturasTotales)
promediolecturastotalesBLD28 <- mean(AlineamientoBLD28$LecturasTotales)
resumenBCD28yBLD28_PorcMapeo <- resumenBCyBLmapeoD28 %>% mutate(PorcentajedeMapeo = c(AlineamientoBCD28$Porcentaje, AlineamientoBLD28$Porcentaje), Resumen = c("Promedios de Alineamiento en Porcentaje", paste0(" BC: ",signif(round(promedioBCD28, digits = 2))), 
                                                                                                                                                            paste0(" BL: ",signif(round(promedioBLD28, digits = 2))),"",
                                                                                                                                                            "Rango de lecturas Mapeadas", "BC: 926 a 2408", "BL: 1860 a 4350", "", "Promedio Lecturas", "BC: 1589", "Promedio Lecturas", "BL: 3471", "2.18x librería alineada", "",
                                                                                                                                                         "Promedio Lecturas Totales", "BC: 404877", "BL: 474250", "14.6% de diferencia en tamaño", rep("", times = 6)))     
View(resumenBCD28yBLD28_PorcMapeo)
promediolecturastotalesBCD28/promediolecturastotalesBLD28


BLGTF <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Genomas referencia/PANDO BL 323/BL 323.gff3"

#Indexar #generar bai ordenado.

#D28 BC
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_1_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_1_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_2_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_2_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_3_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_3_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_4_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_4_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_5_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_5_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_6_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_6_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_7_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_7_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_8_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_8_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_9_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_9_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_10_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_10_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_11_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_11_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/BCD28_12_alineadoaBC.bam", destination = paste(rutaoutputBAMBCD28, "/sortBCD28_12_aBC", sep = ""))

indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_1_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_2_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_3_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_4_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_5_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_6_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_7_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_8_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_9_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_10_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_11_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_12_aBC.bam")

#D28 BL
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_1_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_1_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_2_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_2_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_3_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_3_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_4_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_4_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_5_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_5_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_6_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_6_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_7_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_7_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_8_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_8_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_9_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_9_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_10_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_10_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_11_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_11_aBC", sep = ""))
sortBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/BLD28_12_alineadoaBC.bam", destination = paste(rutaoutputBAMBLD28, "/sortBLD28_12_aBC", sep = ""))

indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_1_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_2_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_3_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_4_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_5_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_6_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_7_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_8_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_9_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_10_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_11_aBC.bam")
indexBam("C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_12_aBC.bam")

# Conteos y Normalización  ------------------------------------------------
#D14 #Might do D14 then D28 BC y D28BL? o algún importe? 
#correrlo en el otro archivo y tenerlo disponible aquí como variable?
#D14 
BL323_1_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD14 BAM/sortBLD1aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
#subsección counts D14
BL323_1_tc_c <-BL323_1_tc$counts 

#D28 BC
BC391d28_1_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_1_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_2_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_2_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_3_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_3_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_4_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_4_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_5_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_5_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_6_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_6_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_7_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_7_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_8_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_8_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_9_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_9_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_10_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_10_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_11_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_11_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BC391d28_12_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BCD28 BAM/sortBCD28_12_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
#D28 BL
BL323d28_1_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_1_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_2_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_2_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_3_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_3_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_4_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_4_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_5_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_5_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_6_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_6_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_7_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_7_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_8_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_8_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_9_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_9_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_10_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_10_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_11_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_11_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
BL323d28_12_tc <- featureCounts(files = "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq/Alineamientos BAM/BLD28 BAM/sortBLD28_12_aBC.bam", #solo se cambia el archivo ej: sortBCD1.bam.bai a sortBCD2.bam.bai
                            annot.ext = BCGFF,
                            isGTFAnnotationFile = TRUE,
                            GTF.featureType = "gene",  # Feature type, e.g., 'exon' for RNA-seq
                            GTF.attrType = "ID",  # Attribute to use for grouping counts (e.g., 'gene_id')
                            useMetaFeatures = TRUE,    # Sum counts over meta-features (e.g., genes)
                            nthreads = 4) 
#subsección counts D28
#D28 BC
BC391d28_1_tc_c <-BC391d28_1_tc$counts 
BC391d28_2_tc_c <-BC391d28_2_tc$counts 
BC391d28_3_tc_c <-BC391d28_3_tc$counts 
BC391d28_4_tc_c <-BC391d28_4_tc$counts 
BC391d28_5_tc_c <-BC391d28_5_tc$counts 
BC391d28_6_tc_c <-BC391d28_6_tc$counts 
BC391d28_7_tc_c <-BC391d28_7_tc$counts 
BC391d28_8_tc_c <-BC391d28_8_tc$counts 
BC391d28_9_tc_c <-BC391d28_9_tc$counts 
BC391d28_10_tc_c <-BC391d28_10_tc$counts 
BC391d28_11_tc_c <-BC391d28_11_tc$counts 
BC391d28_12_tc_c <-BC391d28_12_tc$counts 
#listado BCd29
listBCd28fc_c <- list(BC391d28_1_tc_c,  BC391d28_2_tc_c, BC391d28_3_tc_c,  BC391d28_4_tc_c,
                      BC391d28_5_tc_c,  BC391d28_6_tc_c, BC391d28_7_tc_c,  BC391d28_8_tc_c,
                      BC391d28_9_tc_c, BC391d28_10_tc_c,BC391d28_11_tc_c, BC391d28_12_tc_c)
#D28BL
BL323d28_1_tc_c <- BL323d28_1_tc$counts 
BL323d28_2_tc_c <- BL323d28_2_tc$counts 
BL323d28_3_tc_c <- BL323d28_3_tc$counts 
BL323d28_4_tc_c <- BL323d28_4_tc$counts 
BL323d28_5_tc_c <- BL323d28_5_tc$counts 
BL323d28_6_tc_c <- BL323d28_6_tc$counts 
BL323d28_7_tc_c <- BL323d28_7_tc$counts 
BL323d28_8_tc_c <- BL323d28_8_tc$counts 
BL323d28_9_tc_c <- BL323d28_9_tc$counts 
BL323d28_10_tc_c <-BL323d28_10_tc$counts 
BL323d28_11_tc_c <-BL323d28_11_tc$counts 
BL323d28_12_tc_c <-BL323d28_12_tc$counts 

#listado BLd28
listBLd28fc_c <- list(BL323d28_1_tc_c,  BL323d28_2_tc_c, BL323d28_3_tc_c,  BL323d28_4_tc_c,
                      BL323d28_5_tc_c,  BL323d28_6_tc_c, BL323d28_7_tc_c,  BL323d28_8_tc_c,
                      BL323d28_9_tc_c, BL323d28_10_tc_c,BL323d28_11_tc_c, BL323d28_12_tc_c)

genesBCd28 <- row.names(BC391d28_1_tc_c)
genesBLd28 <- row.names(BL323d28_1_tc_c) #3883 #4002 por referencia a BC
str(genesBLd28)

procFC <- function(x){
  Rn <- row.names(x)
  y <- tibble(x)
  resdf <- as.data.frame(y) %>% mutate(genename = Rn) %>% select(genename, everything())
  return(resdf)
} #En esta función quiero integrar los genes como una variable, cambiar el nombre de los conteos a conteo y especie, filtrar genes no relevantes (será?), 
#para aplicar con todos los datos

#D14
#procesBL <- sapply(listBL1fc_c, procFC)
#
#dfBLd14 <- data.frame(
#  genesBC2 = genesBC, #no se usaron los de BL para normalizar y comparar posteriormente. 
#  BL1T1 = procesBL[,1],
#  BL1T2 = procesBL[,2],
#  BL1T3 = procesBL[,3],
#  BL1T4 = procesBL[,4],
#  BL2T1 = procesBL[,5],
#  BL2T2 = procesBL[,6],
#  BL2T3 = procesBL[,7],
#  BL2T4 = procesBL[,8],
#  BL3T1 = procesBL[,9],
#  BL3T2 = procesBL[,10],
#  BL3T3 = procesBL[,11],
#  BL3T4 = procesBL[,12]
#)
#
#
#dfBLD14_2 <- dfBLd14 %>% select(
#  genesBC2,
#  BL1T1.sortBLD1aBC.bam,
#  BL1T2.sortBLD2aBC.bam,
#  BL1T3.sortBLD3aBC.bam,
#  BL1T4.sortBLD4aBC.bam,
#  BL2T1.sortBLD5aBC.bam,
#  BL2T2.sortBLD6aBC.bam,
#  BL2T3.sortBLD7aBC.bam,
#  BL2T4.sortBLD8aBC.bam,
#  BL3T1.sortBLD9aBC.bam,
#  BL3T2.sortBLD10aBC.bam,
#  BL3T3.sortBLD11aBC.bam,
#  BL3T4.sortBLD12aBC.bam
#)
#view(dfBLD14_2)
#clipr::write_clip(colnames(dfBLd14))
#
#dfBL_c <- dfBLD14_2 %>% rename(
#  Conteo_BL1T1 = BL1T1.sortBLD1aBC.bam,
#  Conteo_BL1T2 = BL1T2.sortBLD2aBC.bam,
#  Conteo_BL1T3 = BL1T3.sortBLD3aBC.bam,
#  Conteo_BL1T4 = BL1T4.sortBLD4aBC.bam,
#  Conteo_BL2T1 = BL2T1.sortBLD5aBC.bam,
#  Conteo_BL2T2 = BL2T2.sortBLD6aBC.bam,
#  Conteo_BL2T3 = BL2T3.sortBLD7aBC.bam,
#  Conteo_BL2T4 = BL2T4.sortBLD8aBC.bam,
#  Conteo_BL3T1 = BL3T1.sortBLD9aBC.bam,
#  Conteo_BL3T2 = BL3T2.sortBLD10aBC.bam,
#  Conteo_BL3T3 = BL3T3.sortBLD11aBC.bam,
#  Conteo_BL3T4 = BL3T4.sortBLD12aBC.bam
#)
#

##D28 BC y BL
procesBCd28 <- sapply(listBCd28fc_c, procFC)
procesBLd28 <- sapply(listBLd28fc_c, procFC)
#BC
dfBCd28 <- data.frame(
  genesBC2 = genesBCd28, #no se usaron los de BL para normalizar y comparar posteriormente. 
  BCd281T1 = procesBCd28[,1],
  BCd281T2 = procesBCd28[,2],
  BCd281T3 = procesBCd28[,3],
  BCd281T4 = procesBCd28[,4],
  BCd282T1 = procesBCd28[,5],
  BCd282T2 = procesBCd28[,6],
  BCd282T3 = procesBCd28[,7],
  BCd282T4 = procesBCd28[,8],
  BCd283T1 = procesBCd28[,9],
  BCd283T2 = procesBCd28[,10],
  BCd283T3 = procesBCd28[,11],
  BCd283T4 = procesBCd28[,12]
)


dfBCD28_2 <- dfBCd28 %>% select(
  genesBC2,
  BCd281T1.sortBCD28_1_aBC.bam,
  BCd281T2.sortBCD28_2_aBC.bam,
  BCd281T3.sortBCD28_3_aBC.bam,
  BCd281T4.sortBCD28_4_aBC.bam,
  BCd282T1.sortBCD28_5_aBC.bam,
  BCd282T2.sortBCD28_6_aBC.bam,
  BCd282T3.sortBCD28_7_aBC.bam,
  BCd282T4.sortBCD28_8_aBC.bam,
  BCd283T1.sortBCD28_9_aBC.bam,
  BCd283T2.sortBCD28_10_aBC.bam,
  BCd283T3.sortBCD28_11_aBC.bam,
  BCd283T4.sortBCD28_12_aBC.bam
)
view(dfBCD28_2)

dfBCd28_c <- dfBCD28_2 %>% rename(Genes1 = genesBC2,
  Conteo_BCd28_1T1 = BCd281T1.sortBCD28_1_aBC.bam,
  Conteo_BCd28_1T2 = BCd281T2.sortBCD28_2_aBC.bam,
  Conteo_BCd28_1T3 = BCd281T3.sortBCD28_3_aBC.bam,
  Conteo_BCd28_1T4 = BCd281T4.sortBCD28_4_aBC.bam,
  Conteo_BCd28_2T1 = BCd282T1.sortBCD28_5_aBC.bam,
  Conteo_BCd28_2T2 = BCd282T2.sortBCD28_6_aBC.bam,
  Conteo_BCd28_2T3 = BCd282T3.sortBCD28_7_aBC.bam,
  Conteo_BCd28_2T4 = BCd282T4.sortBCD28_8_aBC.bam,
  Conteo_BCd28_3T1 = BCd283T1.sortBCD28_9_aBC.bam,
  Conteo_BCd28_3T2 = BCd283T2.sortBCD28_10_aBC.bam,
  Conteo_BCd28_3T3 = BCd283T3.sortBCD28_11_aBC.bam,
  Conteo_BCd28_3T4 = BCd283T4.sortBCD28_12_aBC.bam
)
view(dfBCd28_c)
#dfBLd28 

dfBLd28 <- data.frame(
  genesBC2 = genesBCd28, #no se usaron los de BL para normalizar y comparar posteriormente. 
  BLd281T1 = procesBLd28[,1],
  BLd281T2 = procesBLd28[,2],
  BLd281T3 = procesBLd28[,3],
  BLd281T4 = procesBLd28[,4],
  BLd282T1 = procesBLd28[,5],
  BLd282T2 = procesBLd28[,6],
  BLd282T3 = procesBLd28[,7],
  BLd282T4 = procesBLd28[,8],
  BLd283T1 = procesBLd28[,9],
  BLd283T2 = procesBLd28[,10],
  BLd283T3 = procesBLd28[,11],
  BLd283T4 = procesBLd28[,12]
)

colnames(dfBLd28)

dfBLD28_2 <- dfBLd28 %>% select(
  genesBC2,
  BLd281T1.sortBLD28_1_aBC.bam,
  BLd281T2.sortBLD28_2_aBC.bam,
  BLd281T3.sortBLD28_3_aBC.bam,
  BLd281T4.sortBLD28_4_aBC.bam,
  BLd282T1.sortBLD28_5_aBC.bam,
  BLd282T2.sortBLD28_6_aBC.bam,
  BLd282T3.sortBLD28_7_aBC.bam,
  BLd282T4.sortBLD28_8_aBC.bam,
  BLd283T1.sortBLD28_9_aBC.bam,
  BLd283T2.sortBLD28_10_aBC.bam,
  BLd283T3.sortBLD28_11_aBC.bam,
  BLd283T4.sortBLD28_12_aBC.bam
)
view(dfBLD14_2)

dfBLd28_c <- dfBLD28_2 %>% rename(Genes2 = genesBC2,
  Conteo_BLd28_1T1 = BLd281T1.sortBLD28_1_aBC.bam,
  Conteo_BLd28_1T2 = BLd281T2.sortBLD28_2_aBC.bam,
  Conteo_BLd28_1T3 = BLd281T3.sortBLD28_3_aBC.bam,
  Conteo_BLd28_1T4 = BLd281T4.sortBLD28_4_aBC.bam,
  Conteo_BLd28_2T1 = BLd282T1.sortBLD28_5_aBC.bam,
  Conteo_BLd28_2T2 = BLd282T2.sortBLD28_6_aBC.bam,
  Conteo_BLd28_2T3 = BLd282T3.sortBLD28_7_aBC.bam,
  Conteo_BLd28_2T4 = BLd282T4.sortBLD28_8_aBC.bam,
  Conteo_BLd28_3T1 = BLd283T1.sortBLD28_9_aBC.bam,
  Conteo_BLd28_3T2 = BLd283T2.sortBLD28_10_aBC.bam,
  Conteo_BLd28_3T3 = BLd283T3.sortBLD28_11_aBC.bam,
  Conteo_BLd28_3T4 = BLd283T4.sortBLD28_12_aBC.bam
)
view(dfBLd28_c)
#Unir los marcos D28. 
dfBCD28yBLD28 <- cbind(dfBCd28_c, dfBLd28_c)
dfBCD28yBLD28conteoslimpios <- dfBCD28yBLD28 %>% select(-Genes2)
dfBCD28yBLD28conteoslimpios$conteofilas <- rowSums(dfBCD28yBLD28conteoslimpios[,-1])
dfBCD28yBLD28filtrado <- dfBCD28yBLD28conteoslimpios[dfBCD28yBLD28conteoslimpios$conteofilas >= 1,]
dfBCD28yBLD28filt <- dfBCD28yBLD28filtrado %>% select(-conteofilas)
view(dfBCD28yBLD28filt)
conteogenesBCD28yBLD28 <- dfBCD28yBLD28filt[,-1]
rownames(conteogenesBCD28yBLD28) <- dfBCD28yBLD28filt[,1]

muestrasD28 = colnames(dfBCD28yBLD28filt[,-1])
view(muestrasD28)
view(conteogenesBCD28yBLD28) #1942 *mitad de los datos.

#Tablas de Metadatos para D28

#D28 BC y BL
metadatosBCD28yBLD28 <- data.frame( #Hacer para BLd14 y BLd28
  IDdemuestras = muestrasD28,
  cepa = rep(c("BC391", "BL323"), each = 12),
  replicadoB = rep(1:3, each = 4, times = 2),
  replicadoT = rep(1:4, times = 6)
)
metadatosBCD28yBLD28$replicadoB <- as.factor(metadatosBCD28yBLD28$replicadoB)
metadatosBCD28yBLD28$replicadoT <- as.factor(metadatosBCD28yBLD28$replicadoT)
#Normalizado
ddsBCD28yBLD28 <- DESeqDataSetFromMatrix(countData = conteogenesBCD28yBLD28,
                                   colData = metadatosBCD28yBLD28,
                                   design = ~ cepa + replicadoB)
NormalizadoBCD28yBLD28 <- DESeq(ddsBCD28yBLD28)
str(NormalizadoBCD28yBLD28)
head(counts(NormalizadoBCD28yBLD28, normalized = T))
conteos_normalizadosBCD28yBLD28 <- counts(NormalizadoBCD28yBLD28, normalized = T)
view(conteos_normalizadosBCD28yBLD28)
BCD28yBLD28vst <- varianceStabilizingTransformation(NormalizadoBCD28yBLD28, blind = F)
head(BCD28yBLD28vst)
#PCA D28
PCA1cRbRTD28 <- plotPCA(BCD28yBLD28vst, intgroup = c("cepa", "replicadoB", "replicadoT"), returnData = T) #ambas cepas
str(PCA1cRbRTD28)
plotPCA(BCD28yBLD28vst, intgroup = c("cepa", "replicadoB")) #ambas cepas; menos varianza que primera comparación. 22%. 
plotPCA(BCD28yBLD28vst, intgroup = c("cepa", "replicadoB", "replicadoT")) #ambas cepas; menos varianza que primera comparación. 22%. 
PCAggcepaD28 <- plotPCA(BCD28yBLD28vst, intgroup = "cepa", returnData = T) #se puede ver una amplia dispersión intracepa en BC391, los datos por BL323 están agrupados. #No se está considerando el tamaño de la librería original sino de la librería secundaria las lecturas totales unidas a Mtv para normalizar. La naturaleza de los datos afecta la normalización. La cantidad de lecturas mapeadas a tb por lecturas totales suguieren una mayor expresión de RNA por BL323 al día 14.  
plotPCA(BCD28yBLD28vst, intgroup = "cepa") #se puede ver una amplia dispersión intracepa en BC391, los datos por BL323 están agrupados. #No se está considerando el tamaño de la librería original sino de la librería secundaria las lecturas totales unidas a Mtv para normalizar. La naturaleza de los datos afecta la normalización. La cantidad de lecturas mapeadas a tb por lecturas totales suguieren una mayor expresión de RNA por BL323 al día 14.  

#PCA por GG D28
PCA1forggD28 <- PCA1cRbRTD28
PCA1forggD28$CombinedGroup <- interaction(PCA1forggD28$cepa, PCA1forggD28$replicadoB, PCA1forggD28$replicadoT)
PCA1forggD28$replicadoB <- as.factor(PCA1forggD28$replicadoB)
PCA1forggD28$cepa <- as.factor(PCA1forggD28$cepa)
PCA1forggD28$replicadoT <- as.factor(PCA1forggD28$replicadoT)
pca_ggplot <- ggplot(PCA1forggD28, aes(x = PC1, y = PC2, color = cepa, shape = replicadoB)) +
  geom_point(size = 3) +  # Adjust size as needed
  labs(title = "PCA de archivos BC y BL al día 28", subtitle = "Color Por Cepa, Forma por Replicado Biológico") +
  theme_minimal()
pca_ggplot0 <- ggplot(PCA1forggD28, aes(x = PC1, y = PC2, color = cepa)) +
  geom_point(size = 3) +  # Adjust size as needed
  labs(title = "PCA de archivos BC y BL al día 28", subtitle = "Color Por Cepa") +
  theme_minimal()
pca_ggplot2 <- ggplot(PCA1forggD28, aes(x = PC1, y = PC2, color = cepa, shape = replicadoT)) +
  geom_point(size = 3) + 
  labs(title = "PCA: Color by Strain, Shape by Technical Replicate") +
  theme_minimal()

print(pca_ggplot0)
print(pca_ggplot)
#Volcan y DGE D28 
#colapsar los replicados técnicos. D28
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
ddsBCD28yBLD28colapsadoRT <- collapseReplicates(NormalizadoBCD28yBLD28, groupby = biological_replicates, run = technical_replicates)
#De aquí sale el volcán
resBCD28yBLD28col <- results(ddsBCD28yBLD28colapsadoRT, contrast = c("cepa", "BL323", "BC391")) 
dim(ddsBCD28yBLD28colapsadoRT)
res_pvalD28signif <- resBCD28yBLD28col[which(resBCD28yBLD28col$pvalue < 0.05), ] #estudiar sobre p value y adjusted. 7 de significancia estadística.
res_padjD28signif <- resBCD28yBLD28col[which(resBCD28yBLD28col$padj < 0.05), ] #2 a través de todos los estudios. 

plotMA(resBCD28yBLD28col, main = "MA-plot", ylim = c(-2,2))
vsdBCD28yBLD28 <- varianceStabilizingTransformation(ddsBCD28yBLD28colapsadoRT, blind=FALSE)
colData(vsdBCD28yBLD28)
plotPCA(vsdBCD28yBLD28, intgroup=c("cepa", "replicadoB")) #cambiar a ggplot? #interesante repB1 de BC se podría sacar del estudio. 

top_genesD28 <- head(order(rowVars(assay(vsdBCD28yBLD28)), decreasing=TRUE), 50)
mat <- assay(vsdBCD28yBLD28)[top_genesD28, ]
mat <- mat - rowMeans(mat)  # Center rows
str(mat)
dim(mat)
colnames(mat)
pheatmap(mat)
resBCD28yBLD28col
resBCD28yBLD28v <- na.omit(resBCD28yBLD28col) #seleccionar solo los 35 primeros o mínimo un zoom sobre eso. 
view(resBCD28yBLD28v)
resBCD28yBLD28v$significance <- ifelse(resBCD28yBLD28v$pvalue < 0.05, "Significant", "Not significant") #& abs(resBCyBLv$log2FoldChange) > 1, 
resBCD28yBLD28v$significancepval <- ifelse(resBCD28yBLD28v$significance == "Significant", "Yes", "No" )
#resBCyBLv$ifpval <- 
top_genesBLD28 <- head(resBCD28yBLD28v[order(resBCD28yBLD28v$log2FoldChange), ], 20)  # Adjust number of top genes as needed
top_genesBCD28 <- head(resBCD28yBLD28v[order(resBCD28yBLD28v$log2FoldChange, decreasing = T), ], 20)  # Adjust number of top genes as needed

#top_genes_BLD28_pval_0.05 <- top_genesBLD28[top_genesBLD28$pvalue<0.05]
#top_genes_BCD28_pval_0.05 <- top_genesBCD28[top_genesBCD28$pvalue<0.05]

view(top_genesBLD28)
view(top_genesBCD28)

#view(top_genes_BCD28_pval_0.05)
#view(top_genes_BLD28_pval_0.05)

#Anotación de todos los pvalue en rojo añadir al gráfico. y cambiar eje x a 0.5 de cambio

ggplot(resBCD28yBLD28v, aes(x=log2FoldChange, y=pvalue, color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genesBLD28, aes(label=row.names(top_genesBLD28)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesBCD28, aes(label=row.names(top_genesBCD28)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")+
  #geom_vline(xintercept = c(-0.124, -1.124, -2.124), linetype = "dashed", color = "blue")+
  labs(title="Expresión diferencial BL323 vs BC391 al día 28", x="Log2 Fold Change", y="P-value", subtitle = "Top 10 up/down") +
  theme(legend.position="")

ggplot(resBCD28yBLD28v, aes(x=log2FoldChange, y=pvalue, color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genesBLD28, aes(label=row.names(top_genesBLD28)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesBCD28, aes(label=row.names(top_genesBCD28)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")+
  geom_vline(xintercept = c(-0.124, -1.124, -2.124), linetype = "dashed", color = "blue")+
  annotate("text", x = -0.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
           label = "1", color = "blue", hjust = 1, angle = 0)+
  annotate("text", x = -1.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
           label = "0", color = "blue", hjust = 1, angle = 0)+
  annotate("text", x = -2.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
           label = "-1", color = "blue", hjust = 1, angle = 0)+
  labs(title="Expresión diferencial BL323 vs BC391 al día 28", x="Log2 Fold Change", y="P-value", subtitle = "Conversión -log2(2.18) = -1.124 = 0") +
  theme(legend.position="") #right

ggplot(resBCD28yBLD28v, aes(x=log2FoldChange, y=-log10(pvalue), color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genesBLD28, aes(label=row.names(top_genesBLD28)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesBCD28, aes(label=row.names(top_genesBCD28)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")+
  geom_vline(xintercept = c(-0.124, -1.124, -2.124), linetype = "dashed", color = "blue")+
  annotate("text", x = -0.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
           label = "1", color = "blue", hjust = 1, angle = 0)+
  annotate("text", x = -1.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
           label = "0", color = "blue", hjust = 1, angle = 0)+
  annotate("text", x = -2.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
           label = "-1", color = "blue", hjust = 1, angle = 0)+
  labs(title="Expresión diferencial BL323 vs BC391 al día 28", x="Log2 Fold Change", y="-Log10 P-value", subtitle = "Gráfico de Volcán") +
  theme(legend.position="") #right 

#de 8 a 30 veces más expresado que en BC. los top 10. 
#genes de 2 a 8 veces más expresados en BC que en BL. 
top_genesBLD28mas <- head(resBCD28yBLD28col[order(resBCD28yBLD28col$log2FoldChange), ], 20)  
top_genesBCD28mas <- head(resBCD28yBLD28col[order(resBCD28yBLD28col$log2FoldChange, decreasing = T), ], 20)  
view(top_genesBCD28mas)
view(top_genesBLD28mas)
dfExpresiónmasenBLqueBCD28 <- data.frame(top_genesBLD28mas)
dfExpresiónmasenBCqueBLD28 <- data.frame(top_genesBCD28mas)
view(dfExpresiónmasenBLqueBC)
view(dfExpresiónmasenBCqueBL)





##D14 vs D28 BL
dfBL_c #no tiene datos; curioso
view(dfBL_c)
dfBLd28_c
view(dfBLd28_c)
#unirdatos en una tabla
dfBLD14yBLD28 <- cbind(dfBL_c, dfBLd28_c)
view(dfBLD14yBLD28) #genesBC2
dfBLD14yBLD28conteoslimpios <- dfBLD14yBLD28 %>% select(-Genes2)
dfBLD14yBLD28conteoslimpios$conteofilas <- rowSums(dfBLD14yBLD28conteoslimpios[,-1])
dfBLD14yBLD28filtrado <- dfBLD14yBLD28conteoslimpios[dfBLD14yBLD28conteoslimpios$conteofilas >= 1,]
dfBLD14yBLD28filt <- dfBLD14yBLD28filtrado %>% select(-conteofilas)
muestrasD14yD28 = colnames(dfBLD14yBLD28filt[,-1])
conteogenesBLD14yBLD28 <- dfBLD14yBLD28filt[,-1] #esta es la variable buena para seguir
rownames(conteogenesBLD14yBLD28) <- dfBLD14yBLD28filt[,1]
str(conteogenesBLD14yBLD28)
view(conteogenesBLD14yBLD28)

#Metadatos BLD14 y BLD28 
metadatosBLD14yBLD28 <- data.frame(
  IDdemuestras = muestrasD14yD28,
  Día = rep(c("D14", "D28"), each = 12), #cambiar cepa por Día como variable principal
  replicadoB = rep(1:3, each = 4, times = 2),
  replicadoT = rep(1:4, times = 6)
)
#Normalizado BLD14 y BLD28
ddsBLD14yBLD28 <- DESeqDataSetFromMatrix(countData = conteogenesBLD14yBLD28,
                                   colData = metadatosBLD14yBLD28,
                                   design = ~ Día + replicadoB) #no lo hace porque hay 0s, remover.
NormalizadoBLD14yBLD28 <- DESeq(ddsBLD14yBLD28)
#head(counts(NormalizadoBLD14yBLD28, normalized = T))
conteos_normalizadosBLD14yBLD28 <- counts(NormalizadoBLD14yBLD28, normalized = T)
BLD14yBLD28vst <- varianceStabilizingTransformation(NormalizadoBLD14yBLD28, blind = F)

#PCA D14 y D28
PCA1cRbRTD14yD28 <- plotPCA(BLD14yBLD28vst, intgroup = c("Día", "replicadoB", "replicadoT"), returnData = T) #ambas cepas
str(PCA1cRbRTD28)
plotPCA(BLD14yBLD28vst, intgroup = c("Día", "replicadoB")) # 
plotPCA(BLD14yBLD28vst, intgroup = c("Día", "replicadoB", "replicadoT")) #grupo 2 al día 14 es el que más varía.  

PCAggcepaD14yD28 <- plotPCA(BLD14yBLD28vst, intgroup = "Día", returnData = T) #  
plotPCA(BLD14yBLD28vst, intgroup = "Día") 

#PCA por GG D14 y D28
PCA1forggD14yD28 <- PCA1cRbRTD14yD28
PCA1forggD14yD28$CombinedGroup <- interaction(PCA1forggD14yD28$Día, PCA1forggD14yD28$replicadoB, PCA1forggD14yD28$replicadoT)
PCA1forggD14yD28$replicadoB <- as.factor(PCA1forggD14yD28$replicadoB)
PCA1forggD14yD28$Día <- as.factor(PCA1forggD14yD28$Día)
PCA1forggD14yD28$replicadoT <- as.factor(PCA1forggD14yD28$replicadoT)
pca_ggplot <- ggplot(PCA1forggD14yD28, aes(x = PC1, y = PC2, color = Día, shape = replicadoB)) +
  geom_point(size = 3) +  # Adjust size as needed
  labs(title = "PCA de archivos BL al día 14 y BL al día 28", subtitle = "Color Por Día, Forma por Replicado Biológico") +
  theme_minimal()
pca_ggplot0 <- ggplot(PCA1forggD14yD28, aes(x = PC1, y = PC2, color = Día)) +
  geom_point(size = 3) +  # Adjust size as needed
  labs(title = "PCA de archivos BL día 14 y BL al día 28", subtitle = "Color Por Día") +
  theme_minimal()
pca_ggplot2 <- ggplot(PCA1forggD14yD28, aes(x = PC1, y = PC2, color = cepa, shape = replicadoT)) +
  geom_point(size = 3) + 
  labs(title = "PCA: Color by Strain, Shape by Technical Replicate") +
  theme_minimal()

print(pca_ggplot0)
print(pca_ggplot)

#DGE y volcán D14 y D28.
NormalizadoBLD14yBLD28
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
#strainBLD14yBLD28 <- c("BL323", "BL323", "BL323", "BL323",
#                       "BL323", "BL323", "BL323", "BL323",
#                       "BL323", "BL323", "BL323", "BL323",
#                       "BL323", "BL323", "BL323", "BL323",
#                       "BL323", "BL323", "BL323", "BL323",
#                       "BL323", "BL323", "BL323", "BL323")
Día <- c(rep("D14", times = 12), rep("D28", times = 12))

# DGE/Análisis Transcriptómico BLd28 vs BLd14 ------------------------------------------------------
NormalizadoBLD14yBLD28
str(NormalizadoBLD14yBLD28)

#colapsar los replicados técnicos. 

ddsBLD14yBLD28colapsadoRT <- collapseReplicates(NormalizadoBLD14yBLD28, groupby = biological_replicates, run = technical_replicates)
#De aquí sale el volcán
resBLD14yBLD28col <- results(ddsBLD14yBLD28colapsadoRT, contrast = c("Día", "D14", "D28")) #relativo al 14. a la izquierda lo que se expresa menos y a la derecha más.

dim(ddsBLD14yBLD28colapsadoRT)
res_pvalBLD14vsD28signif <- resBLD14yBLD28col[which(resBLD14yBLD28col$pvalue < 0.05), ] #estudiar sobre p value y adjusted. 7 de significancia estadística.
res_padjBLD14vsD28signif <- resBLD14yBLD28col[which(resBLD14yBLD28col$padj < 0.05), ] #2 a través de todos los estudios. 

plotMA(resBLD14yBLD28col, main = "MA-plot", ylim = c(-2,2))
vsdBLD14yBLD28 <- varianceStabilizingTransformation(ddsBLD14yBLD28colapsadoRT, blind=FALSE)
colData(vsdBLD14yBLD28)
plotPCA(vsdBLD14yBLD28, intgroup=c("Día", "replicadoB")) #cambiar a ggplot? #interesante repB1 de BC se podría sacar del estudio. 

top_genesBLD14yD28 <- head(order(rowVars(assay(vsdBLD14yBLD28)), decreasing=TRUE), 50)
mat <- assay(vsdBLD14yBLD28)[top_genesBLD14yD28, ]
mat <- mat - rowMeans(mat)  # Center rows
str(mat)
dim(mat)
colnames(mat)
pheatmap(mat)

#Volcán BCD28 vs BLD28
resBLD14yBLD28col
resBLD14yBLD28v <- na.omit(resBLD14yBLD28col) #seleccionar solo los 35 primeros o mínimo un zoom sobre eso. 
view(resBLD14yBLD28v)
resBLD14yBLD28v$significance <- ifelse(resBLD14yBLD28v$pvalue < 0.05, "Significant", "Not significant") #& abs(resBCyBLv$log2FoldChange) > 1, 
resBLD14yBLD28v$significancepval <- ifelse(resBLD14yBLD28v$significance == "Significant", "Yes", "No" )
#resBCyBLv$ifpval <- 
top_genesBLD14v28 <- head(resBLD14yBLD28v[order(resBLD14yBLD28v$log2FoldChange), ], 10)  # Adjust number of top genes as needed
top_genesBLD28v14 <- head(resBLD14yBLD28v[order(resBLD14yBLD28v$log2FoldChange, decreasing = T), ], 10)  # Adjust number of top genes as needed
view(top_genesBLD14v28)
view(top_genesBLD28v14)

#Anotación de todos los pvalue en rojo añadir al gráfico. cambiar el eje x por 0.5 de cambios para ver mejor.

ggplot(resBLD14yBLD28v, aes(x=log2FoldChange, y=pvalue, color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genesBLD14v28, aes(label=row.names(top_genesBLD14v28)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesBLD28v14, aes(label=row.names(top_genesBLD28v14)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")+
  #geom_vline(xintercept = c(-0.124, -1.124, -2.124), linetype = "dashed", color = "blue")+
  labs(title="Expresión diferencial BL323 al día 14 vs BL323 al día 28", x="Log2 Fold Change", y="P-value", subtitle = "Top 10 up/down") +
  theme(legend.position="")
#factor normalización. BLD14 es .808 de D28. D28 es 1.237 de BLD14. Entonces multiplica BLD14 por 1.237

#gráfico no relevante normalización adecuada
ggplot(resBLD14yBLD28v, aes(x=log2FoldChange, y=pvalue, color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genesBLD14v28, aes(label=row.names(top_genesBLD14v28)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesBLD28v14, aes(label=row.names(top_genesBLD28v14)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")+
  #geom_vline(xintercept = c(-1.21, -0.21, .79), linetype = "dashed", color = "blue")+
  #annotate("text", x = -1.21, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
  #         label = "-1", color = "blue", hjust = 1, angle = 0)+
  #annotate("text", x = -.21, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
  #         label = "0", color = "blue", hjust = 1, angle = 0)+
  #annotate("text", x = .79, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
  #         label = "1", color = "blue", hjust = 1, angle = 0)+
  labs(title="Expresión diferencial BL323 al día 14 vs BL323 al día 28", x="Log2 Fold Change", y="P-value", subtitle = "Conversión log2(1.237) = 0.21 = 0") +
  theme(legend.position="") #right

ggplot(resBLD14yBLD28v, aes(x=log2FoldChange, y=-log10(pvalue), color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genesBLD14v28, aes(label=row.names(top_genesBLD14v28)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesBLD28v14, aes(label=row.names(top_genesBLD28v14)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")+
  #geom_vline(xintercept = c(-0.124, -1.124, -2.124), linetype = "dashed", color = "blue")+
  #annotate("text", x = -0.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
  #         label = "1", color = "blue", hjust = 1, angle = 0)+
  #annotate("text", x = -1.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
  #         label = "0", color = "blue", hjust = 1, angle = 0)+
  #annotate("text", x = -2.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
  #         label = "-1", color = "blue", hjust = 1, angle = 0)+
  labs(title="Expresión diferencial BL323 al día 14 vs BL323 al día 28", x="Log2 Fold Change", y="-Log10 P-value", subtitle = "Gráfico de Volcán") +
  theme(legend.position="") #right


#de 8 a 30 veces más expresado que en BC. los top 10. 
#genes de 2 a 8 veces más expresados en BC que en BL. 
top_genesBLD14masD28 <- head(resBCD28yBLD28col[order(resBCD28yBLD28col$log2FoldChange), ], 20)  
top_genesBLD28masD14 <- head(resBCD28yBLD28col[order(resBCD28yBLD28col$log2FoldChange, decreasing = T), ], 20)  
view(top_genesBLD14masD28)
view(top_genesBLD28masD14)
dfExpresiónmasenBL14queBL28 <- data.frame(top_genesBLD14masD28)
dfExpresiónmasenBL28queBL14 <- data.frame(top_genesBLD28masD14)
view(dfExpresiónmasenBL14queBL28)
view(dfExpresiónmasenBL28queBL14)

#D14 y D28. BC.
dfBC_c #D14
dfBCd28_c
dfBCD14yBCD28 <- cbind(dfBC_c, dfBCd28_c)
view(dfBCD14yBCD28) #genesBC2
dfBCD14yBCD28conteoslimpios <- dfBCD14yBCD28 %>% select(-Genes1)
dfBCD14yBCD28conteoslimpios$conteofilas <- rowSums(dfBCD14yBCD28conteoslimpios[,-1])
dfBCD14yBCD28filtrado <- dfBCD14yBCD28conteoslimpios[dfBCD14yBCD28conteoslimpios$conteofilas >= 1,]
dfBCD14yBCD28filt <- dfBCD14yBCD28filtrado %>% select(-conteofilas)
muestrasBCD14yD28 = colnames(dfBCD14yBCD28filt[,-1])
conteogenesBCD14yBCD28 <- dfBCD14yBCD28filt[,-1] #esta es la variable buena para seguir
rownames(conteogenesBCD14yBCD28) <- dfBCD14yBCD28filt[,1]
#Metadatos
metadatosBCD14yBCD28 <- data.frame(
  IDdemuestras = muestrasBCD14yD28,
  Día = rep(c("D14", "D28"), each = 12), #cambiar cepa por Día como variable principal
  replicadoB = rep(1:3, each = 4, times = 2),
  replicadoT = rep(1:4, times = 6)
)
#Normalizado BCD14 y BCD28
ddsBCD14yBCD28 <- DESeqDataSetFromMatrix(countData = conteogenesBCD14yBCD28,
                                         colData = metadatosBCD14yBCD28,
                                         design = ~ Día + replicadoB) #no lo hace porque hay 0s, remover.
NormalizadoBCD14yBCD28 <- DESeq(ddsBCD14yBCD28)
#head(counts(NormalizadoBCD14yBCD28, normalized = T))
conteos_normalizadosBCD14yBCD28 <- counts(NormalizadoBCD14yBCD28, normalized = T)
BCD14yBCD28vst <- varianceStabilizingTransformation(NormalizadoBCD14yBCD28, blind = F)

#PCA BC D14 y D28
PCA1cRbRTBCD14yD28 <- plotPCA(BCD14yBCD28vst, intgroup = c("Día", "replicadoB", "replicadoT"), returnData = T) #ambas cepas
str(PCA1cRbRTBCD14D28)
plotPCA(BCD14yBCD28vst, intgroup = c("Día", "replicadoB")) # 
plotPCA(BCD14yBCD28vst, intgroup = c("Día", "replicadoB", "replicadoT")) #grupo 2 al día 14 es el que más varía.  

PCAggcepaBCD14yD28 <- plotPCA(BCD14yBCD28vst, intgroup = "Día", returnData = T) #  
plotPCA(BCD14yBCD28vst, intgroup = "Día") 

#PCA por GG BC D14 y D28
PCA1forggBCD14yD28 <- PCA1cRbRTBCD14yD28
PCA1forggBCD14yD28$CombinedGroup <- interaction(PCA1forggBCD14yD28$Día, PCA1forggBCD14yD28$replicadoB, PCA1forggBCD14yD28$replicadoT)
PCA1forggBCD14yD28$replicadoB <- as.factor(PCA1forggBCD14yD28$replicadoB)
PCA1forggBCD14yD28$Día <- as.factor(PCA1forggBCD14yD28$Día)
PCA1forggBCD14yD28$replicadoT <- as.factor(PCA1forggBCD14yD28$replicadoT)
pca_ggplot <- ggplot(PCA1forggBCD14yD28, aes(x = PC1, y = PC2, color = Día, shape = replicadoB)) +
  geom_point(size = 3) +  # Adjust size as needed
  labs(title = "PCA de archivos BC al día 14 y BC al día 28", subtitle = "Color Por Día, Forma por Replicado Biológico") +
  theme_minimal()
pca_ggplot0 <- ggplot(PCA1forggBCD14yD28, aes(x = PC1, y = PC2, color = Día)) +
  geom_point(size = 3) +  # Adjust size as needed
  labs(title = "PCA de archivos BC día 14 y BC al día 28", subtitle = "Color Por Día") +
  theme_minimal()
pca_ggplot2 <- ggplot(PCA1forggBCD14yD28, aes(x = PC1, y = PC2, color = cepa, shape = replicadoT)) +
  geom_point(size = 3) + 
  labs(title = "PCA: Color by Strain, Shape by Technical Replicate") +
  theme_minimal()

print(pca_ggplot0)
print(pca_ggplot)
#obtener factor de normalización
#DGE volcán BCD14 vs BCD28
NormalizadoBCD14yBCD28
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
#strainBLD14yBLD28 <- c("BL323", "BL323", "BL323", "BL323",
#                       "BL323", "BL323", "BL323", "BL323",
#                       "BL323", "BL323", "BL323", "BL323",
#                       "BL323", "BL323", "BL323", "BL323",
#                       "BL323", "BL323", "BL323", "BL323",
#                       "BL323", "BL323", "BL323", "BL323")
Día <- c(rep("D14", times = 12), rep("D28", times = 12))

# DGE/Análisis Transcriptómico BCd28 vs BCd14 ------------------------------------------------------
NormalizadoBCD14yBCD28
str(NormalizadoBCD14yBCD28)

#colapsar los replicados técnicos. 

ddsBCD14yBCD28colapsadoRT <- collapseReplicates(NormalizadoBCD14yBCD28, groupby = biological_replicates, run = technical_replicates)
#De aquí sale el volcán
resBCD14yBCD28col <- results(ddsBCD14yBCD28colapsadoRT, contrast = c("Día", "D14", "D28")) #relativo al 14. a la izquierda lo que se expresa menos y a la derecha más.

dim(ddsBCD14yBCD28colapsadoRT)
res_pvalBCD14vsD28signif <- resBCD14yBCD28col[which(resBCD14yBCD28col$pvalue < 0.05), ] #estudiar sobre p value y adjusted. 7 de significancia estadística.
res_padjBCD14vsD28signif <- resBCD14yBCD28col[which(resBCD14yBCD28col$padj < 0.05), ] #2 a través de todos los estudios. 

plotMA(resBCD14yBCD28col, main = "MA-plot", ylim = c(-2,2))
vsdBCD14yBCD28 <- varianceStabilizingTransformation(ddsBCD14yBCD28colapsadoRT, blind=FALSE)
colData(vsdBCD14yBCD28)
plotPCA(vsdBCD14yBCD28, intgroup=c("Día", "replicadoB")) #cambiar a ggplot? #interesante repB1 de BC se podría sacar del estudio. 

top_genesBCD14yD28 <- head(order(rowVars(assay(vsdBCD14yBCD28)), decreasing=TRUE), 50)
mat <- assay(vsdBCD14yBCD28)[top_genesBCD14yD28, ]
mat <- mat - rowMeans(mat)  # Center rows
str(mat)
dim(mat)
colnames(mat)
pheatmap(mat)

#Volcán BCD14 vs BCD28
resBCD14yBCD28col
resBCD14yBCD28v <- na.omit(resBCD14yBCD28col) #seleccionar solo los 35 primeros o mínimo un zoom sobre eso. 
view(resBCD14yBCD28v)
resBCD14yBCD28v$significance <- ifelse(resBCD14yBCD28v$pvalue < 0.05, "Significant", "Not significant") #& abs(resBCyBLv$log2FoldChange) > 1, 
resBCD14yBCD28v$significancepval <- ifelse(resBCD14yBCD28v$significance == "Significant", "Yes", "No" )
#resBCyBLv$ifpval <- 
top_genesBCD14v28 <- head(resBCD14yBCD28v[order(resBCD14yBCD28v$log2FoldChange), ], 10)  # Adjust number of top genes as needed
top_genesBCD28v14 <- head(resBCD14yBCD28v[order(resBCD14yBCD28v$log2FoldChange, decreasing = T), ], 10)  # Adjust number of top genes as needed
view(top_genesBCD14v28)
view(top_genesBCD28v14)

#Anotación de todos los pvalue en rojo añadir al gráfico. cambiar el eje x por 0.5 de cambios para ver mejor.

ggplot(resBCD14yBCD28v, aes(x=log2FoldChange, y=pvalue, color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genesBCD14v28, aes(label=row.names(top_genesBCD14v28)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesBCD28v14, aes(label=row.names(top_genesBCD28v14)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")+
  #geom_vline(xintercept = c(-0.124, -1.124, -2.124), linetype = "dashed", color = "blue")+
  labs(title="Expresión diferencial BC391 al día 14 vs BC391 al día 28", x="Log2 Fold Change", y="P-value", subtitle = "Top 10 up/down") +
  theme(legend.position="")
#factor normalización. BLD14 es .808 de D28. D28 es 1.237 de BLD14. Entonces multiplica BLD14 por 1.237

ggplot(resBCD14yBCD28v, aes(x=log2FoldChange, y=-log10(pvalue), color = significancepval)) + #y=-log10(pvalue)
  geom_point(alpha=0.5, size=2) +
  geom_text(data=top_genesBCD14v28, aes(label=row.names(top_genesBCD14v28)), vjust=1, hjust=1, size=3) +
  geom_text(data = top_genesBCD28v14, aes(label=row.names(top_genesBCD28v14)), vjust=1, hjust=1, size=3) +
  scale_color_manual(values = c("Yes" = "red", "No"= "black")) +  # Red for significant, black for non-significant
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")+
  #geom_vline(xintercept = c(-0.124, -1.124, -2.124), linetype = "dashed", color = "blue")+
  #annotate("text", x = -0.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
  #         label = "1", color = "blue", hjust = 1, angle = 0)+
  #annotate("text", x = -1.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
  #         label = "0", color = "blue", hjust = 1, angle = 0)+
  #annotate("text", x = -2.124, y = max(resBCD28yBLD28v$pvalue, na.rm = TRUE), 
  #         label = "-1", color = "blue", hjust = 1, angle = 0)+
  labs(title="Expresión diferencial BC391 al día 14 vs BC391 al día 28", x="Log2 Fold Change", y="-Log10 P-value", subtitle = "Gráfico de Volcán") +
  theme(legend.position="") #right


#Listas DGE BCD14 vs BCD28 
top_genesBCD14masD28 <- head(resBCD14yBCD28col[order(resBCD14yBCD28col$log2FoldChange), ], 20)  
top_genesBCD28masD14 <- head(resBCD14yBCD28col[order(resBCD14yBCD28col$log2FoldChange, decreasing = T), ], 20)  
view(top_genesBCD14masD28)
view(top_genesBCD28masD14)
dfExpresiónmasenBC14queBC28 <- data.frame(top_genesBLD14masD28)
dfExpresiónmasenBC28queBC14 <- data.frame(top_genesBLD28masD14)
view(dfExpresiónmasenBC14queBC28)
view(dfExpresiónmasenBC28queBC14)

#by time 14

#BC391
BCvsBL_D14 <- dfExpresiónmasenBLqueBC #las sobrerepresentadas de BC
#BL323
BLvsBC_D14 <- dfExpresiónmasenBCqueBL #las sobrerepresentadas de BL

#by time 28
#BC391
BCvsBL_D28 <- dfExpresiónmasenBCqueBLD28 #sobrerepresentados en BCD28 vs BLD28
View(BCvsBL_D28)
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


