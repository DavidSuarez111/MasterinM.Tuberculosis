
# Importar librerías ------------------------------------------------------

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


# Direcciones a genoma referencia y archivos ------------------------------
directorio_general <- "C:/Users/david/1.Archivos-Directorios/1. Maestría/0Tesis Maestría/1. AT en R/1. Datos limpios fastq"

#Genoma referencia
genoma_referencia_balbc <- paste(directorio_general, "/Genomas referencia/balbc/balbc.fasta", sep="")
anotaciones_balbc <- paste(directorio_general, "/Genomas referencia/balbc/balbc_annotations.gff3", sep = "")

#Archivos RNAseq o .fastq. Cada variable tiene 12 archivos revisables por list.files()
transcriptomas_murinos_infeccion_bc_d14 <- paste(directorio_general,"/BC d14", sep = "")
transcriptomas_murinos_infeccion_bc_d28 <-paste(directorio_general,"/BC d28", sep = "")
transcriptomas_murinos_infeccion_bl_d14 <-paste(directorio_general,"/BL d14", sep = "")
transcriptomas_murinos_infeccion_bl_d28 <-paste(directorio_general,"/BL d28", sep = "")
#list.files(), list.dirs() #funciones para revisar los archivos en cada directorio. 

#Directorios donde se pondrán los archivos de indexados
balbc_index <- paste(directorio_general, "/Genomas referencia/balbc", "/balbc_index", sep = "")

#Directorios donde se pondrán los archivos alineados (.bam)
directorio_bc_d14_bam <- paste(directorio_general, "/Alineamientos BAM/Mice analyses/BC d14 BAM", sep = "")
directorio_bc_d28_bam <- paste(directorio_general, "/Alineamientos BAM/Mice analyses/BC d28 BAM", sep = "")
directorio_bl_d14_bam <- paste(directorio_general, "/Alineamientos BAM/Mice analyses/BL d14 BAM", sep = "")
directorio_bl_d28_bam <- paste(directorio_general, "/Alineamientos BAM/Mice analyses/BL d28 BAM", sep = "")

#Generar el indexado del genoma murino balbc
buildindex( #balbc
  basename = balbc_index, 
  reference = genoma_referencia_balbc)


# Alineamientos a ratón ---------------------------------------------------
#.fastq a .bam
#Se ocupa el genoma indexado y las anotaciones. 

#12 alineamiento BC d14. Cada archivo fastq se alínea al genoma indexado de balbc para poder pasar a contar la expresión genética.
align( #1 -Mus musculus balbc BC d14 1er archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR1T1 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_1.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d14 2do archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR1T2 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_2.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d14 3er archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR1T3 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_3.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d14 4to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR1T4 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_4.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BC d14 5to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR2T1 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_5.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d14 6to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR2T2 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_6.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BC d14 7o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR2T3 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcd14_7.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d14 8vo archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR2T4 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_8.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BC d14 9no archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR3T1 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_9.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d14 10o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR3T2 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_10.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BC d14 11vo archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR3T3 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_11.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #12 -Mus musculus balbc BC d14 12o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d14, "/SRR BR3T4 BC D14.fastq", sep = ""),
  output_file = paste(directorio_bc_d14_bam, "/balbcbcd14_12.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

#12 alineamientos BL d14

align( #1 -Mus musculus balbc BL d14 1er archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR1T1 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_1.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d14 2do archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR1T2 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_2.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d14 3er archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR1T3 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_3.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d14 4to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR1T4 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_4.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BL d14 5to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR2T1 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_5.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d14 6to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR2T2 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_6.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BL d14 7o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR2T3 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbld14_7.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d14 8vo archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR2T4 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_8.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BL d14 9no archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR3T1 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_9.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d14 10o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR3T2 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_10.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BL d14 11vo archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR3T3 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_11.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #12 -Mus musculus balbc BL d14 12o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR3T4 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d14_bam, "/balbcbld14_12.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

#12 alineamientos BC d28

align( #1 -Mus musculus balbc BC d28 1er archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR1T1 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_1.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d28 2do archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR1T2 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_2.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d28 3er archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR1T3 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_3.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d28 4to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR1T4 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_4.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BC d28 5to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR2T1 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_5.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d28 6to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR2T2 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_6.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BC d28 7o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR2T3 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcd28_7.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d28 8vo archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR2T4 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_8.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BC d28 9no archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR3T1 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_9.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BC d28 10o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR3T2 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_10.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BC d28 11vo archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28 , "/SRR BR3T3 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_11.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #12 -Mus musculus balbc BC d28 12o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bc_d28, "/SRR BR3T4 BC D28.fastq", sep = ""),
  output_file = paste(directorio_bc_d28_bam, "/balbcbcd28_12.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)
#12 alineamientos BL d28

align( #1 -Mus musculus balbc BL d28 1er archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d28, "/SRR BR1T1 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_1.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d28 2do archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d28, "/SRR BR1T2 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_2.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d28 3er archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d28, "/SRR BR1T3 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_3.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d28 4to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d28, "/SRR BR1T4 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_4.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BL d28 5to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d28, "/SRR BR2T1 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_5.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d14 6to archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR2T2 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_6.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BL d28 7o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d28, "/SRR BR2T3 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbld28_7.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d28 8vo archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d28, "/SRR BR2T4 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_8.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BL d28 9no archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d28, "/SRR BR3T1 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_9.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #1 -Mus musculus balbc BL d28 10o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d14, "/SRR BR3T2 BL D14.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_10.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)

align( #1 -Mus musculus balbc BL d14 11vo archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d28, "/SRR BR3T3 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_11.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
) 

align( #12 -Mus musculus balbc BL d14 12o archivo
  index = balbc_index,
  readfile1 = paste(transcriptomas_murinos_infeccion_bl_d28, "/SRR BR3T4 BL D28.fastq", sep = ""),
  output_file = paste(directorio_bl_d28_bam, "/balbcbld28_12.bam", sep = ""),
  nthreads = 4,
  input_format = "FASTQ"
)


# Resumen estadístico de alineamientos ------------------------------------
#4 tablas

#1 BC D14

#2 BC D14

#3 BL D14

#4 BL D28

# Ordenar los alineamientos .bam -------------------------------
#12 archivos BC d14
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_1.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_1.bam", sep = ""), "/sortBCD1", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_2.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_2.bam", sep = ""), "/sortBCD2", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_3.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_3.bam", sep = ""), "/sortBCD3", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_4.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_4.bam", sep = ""), "/sortBCD4", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_5.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_5.bam", sep = ""), "/sortBCD5", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_6.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_6.bam", sep = ""), "/sortBCD6", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_7.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_7.bam", sep = ""), "/sortBCD7", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_8.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_8.bam", sep = ""), "/sortBCD8", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_9.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_9.bam", sep = ""), "/sortBCD9", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_10.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_10.bam", sep = ""), "/sortBCD10", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_11.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_11.bam", sep = ""), "/sortBCD11", sep = ""))
sortBam(paste(directorio_bc_d14_bam, "/balbcbcd14_11.bam", sep = ""), destination = paste(paste(directorio_bc_d14_bam, "/balbcbcd14_12.bam", sep = ""), "/sortBCD12", sep = ""))

#12 archivos BC d28
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_1.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_1.bam", sep = ""), "/sortBCD1", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_2.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_2.bam", sep = ""), "/sortBCD2", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_3.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_3.bam", sep = ""), "/sortBCD3", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_4.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_4.bam", sep = ""), "/sortBCD4", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_5.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_5.bam", sep = ""), "/sortBCD5", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_6.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_6.bam", sep = ""), "/sortBCD6", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_7.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_7.bam", sep = ""), "/sortBCD7", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_8.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_8.bam", sep = ""), "/sortBCD8", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_9.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_9.bam", sep = ""), "/sortBCD9", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_10.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_10.bam", sep = ""), "/sortBCD10", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_11.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_11.bam", sep = ""), "/sortBCD11", sep = ""))
sortBam(paste(directorio_bc_d28_bam, "/balbcbcd28_12.bam", sep = ""), destination = paste(paste(directorio_bc_d28_bam, "/balbcbcd28_12.bam", sep = ""), "/sortBCD12", sep = ""))
#12 archivos BL d14
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_1.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_1.bam", sep = ""), "/sortBLD1", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_2.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_2.bam", sep = ""), "/sortBLD2", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_3.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_3.bam", sep = ""), "/sortBLD3", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_4.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_4.bam", sep = ""), "/sortBLD4", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_5.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_5.bam", sep = ""), "/sortBLD5", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_6.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_6.bam", sep = ""), "/sortBLD6", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_7.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_7.bam", sep = ""), "/sortBLD7", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_8.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_8.bam", sep = ""), "/sortBLD8", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_9.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_9.bam", sep = ""), "/sortBLD9", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_10.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_10.bam", sep = ""), "/sortBLD10", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_11.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_11.bam", sep = ""), "/sortBLD11", sep = ""))
sortBam(paste(directorio_bl_d14_bam, "/balbcbld14_11.bam", sep = ""), destination = paste(paste(directorio_bl_d14_bam, "/balbcbld14_12.bam", sep = ""), "/sortBLD12", sep = ""))
#12 archivos BL d28
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_1.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_1.bam", sep = ""), "/sortBLD1", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_2.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_2.bam", sep = ""), "/sortBLD2", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_3.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_3.bam", sep = ""), "/sortBLD3", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_4.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_4.bam", sep = ""), "/sortBLD4", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_5.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_5.bam", sep = ""), "/sortBLD5", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_6.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_6.bam", sep = ""), "/sortBLD6", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_7.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_7.bam", sep = ""), "/sortBLD7", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_8.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_8.bam", sep = ""), "/sortBLD8", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_9.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_9.bam", sep = ""), "/sortBLD9", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_10.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_10.bam", sep = ""), "/sortBLD10", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_11.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_11.bam", sep = ""), "/sortBLD11", sep = ""))
sortBam(paste(directorio_bl_d28_bam, "/balbcbld28_12.bam", sep = ""), destination = paste(paste(directorio_bl_d28_bam, "/balbcbld28_12.bam", sep = ""), "/sortBLD12", sep = ""))

# Indexación de archivos .bam a .bam.bai ----------------------------------
#12 archivos BC d14
indexBam(paste(directorio_bc_d14_bam, "/sortBCD1",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD2",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD3",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD4",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD5",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD6",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD7",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD8",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD9",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD10",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD11",sep = ""))
indexBam(paste(directorio_bc_d14_bam, "/sortBCD12",sep = ""))
#12 archivos BC d28
indexBam(paste(directorio_bc_d28_bam, "/sortBCD1",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD2",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD3",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD4",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD5",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD6",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD7",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD8",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD9",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD10",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD11",sep = ""))
indexBam(paste(directorio_bc_d28_bam, "/sortBCD12",sep = ""))
#12 archivos BL d14
indexBam(paste(directorio_bl_d14_bam, "/sortBLD1",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD2",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD3",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD4",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD5",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD6",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD7",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD8",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD9",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD10",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD11",sep = ""))
indexBam(paste(directorio_bl_d14_bam, "/sortBLD12",sep = ""))
#12 archivos BL d28
indexBam(paste(directorio_bl_d28_bam, "/sortBLD1",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD2",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD3",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD4",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD5",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD6",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD7",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD8",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD9",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD10",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD11",sep = ""))
indexBam(paste(directorio_bl_d28_bam, "/sortBLD12",sep = ""))


# Conteo de transcritos -featureCounts() ------------------------------------
#podría definir una función con feature counts adentro y que le pase los valores por defecto de todo lo dicho y solo que el primero de la dirección sea asignado y usarlo en lapply

#12 archivos BC d14


#12 archivos BC d28

#12 archivos BL d14

#12 archivos BL d28

# Unir los marcos de datos ------------------------------------------------


# Normalización -----------------------------------------------------------


# Volcan ------------------------------------------------------------------


# Filtrado ----------------------------------------------------------------


# Exportación a excel -----------------------------------------------------








