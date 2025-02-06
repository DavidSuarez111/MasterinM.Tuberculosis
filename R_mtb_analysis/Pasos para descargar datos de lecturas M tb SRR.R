##Documento con el propósito de descargar los datos usando SRAtoolkit en terminal.
#El documento debe generar la unión de la función a ser ejecutada junto a el nombre de cada archivo para que al correrlo en la terminal se descarguen los archivos de interés (de la cepa Beijing Clásica 391 y Tipo Beijing 323 al día 14) en formato fastq.

##Pasos para descargar datos de NCBI

#1. Variables con nombre de los archivos por cepa por día


numaccesNCBIBCd14 <- c("SRR18348642",
                       "SRR18348641",
                       "SRR18348640",
                       "SRR18348639",
                       "SRR18348638",
                       "SRR18348637",
                       "SRR18348636",
                       "SRR18348666",
                       "SRR18348665",
                       "SRR18348633",
                       "SRR18348672",
                       "SRR18348671")

numaccesNCBIBLd14 <- c("SRR18348652",
                       "SRR18348623",
                       "SRR18348622",
                       "SRR18348651",
                       "SRR18348650",
                       "SRR18348649",
                       "SRR18348628",
                       "SRR18348627",
                       "SRR18348626",
                       "SRR18348625",
                       "SRR18348624",
                       "SRR18348674")

numaccesNCBIBCd28 <- c("SRR18348670",
                       "SRR18348669",
                       "SRR18348668",
                       "SRR18348667",
                       "SRR18348605",
                       "SRR18348664",
                       "SRR18348663",
                       "SRR18348662",
                       "SRR18348612",
                       "SRR18348611",
                       "SRR18348610",
                       "SRR18348609")

numaccesNCBIBLd28 <- c("SRR18348673",
                       "SRR18348620",
                       "SRR18348621",
                       "SRR18348680",
                       "SRR18348679",
                       "SRR18348678",
                       "SRR18348677",
                       "SRR18348676",
                       "SRR18348675",
                       "SRR18348648",
                       "SRR18348646",
                       "SRR18348645")
#

#Variables de las listas de los números de acceso a los archivos transcriptómicos por cepa y día
#Cepa BC y BL al día 14
numaccesNCBIBCd14

numaccesNCBIBLd14

#Cepa BC y BL al día 28
numaccesNCBIBCd28

numaccesNCBIBLd28


##Generar uniones de las funciones con los nombres de los archivos 

n1 <- c()
n2 <- c()

#func_y_doc_BCd14 
for (n in numaccesNCBIBCd28) {
  n1 <- c(paste0("fastq-dump --split-3 ", n))
  print(n1)
}

#func_y_doc_BLd14
for(n in numaccesNCBIBLd28) {
  n2 <- c(paste0("fastq-dump --split-3 ", n))
  print(n2)
}

print(n1) #BCd14 (listo, descargado por fastq-dump)
print(n2) #BLd14


#Pegar los valores de n1 y n2 en la terminal 1 por 1. 
#Se descargaron con éxito los archivos usando terminal en apoyo del programa SRAtoolkit usando las lineas impresas por n1 a n2. 

#n3 <- c()
#n4 <- c()
#func_y_doc_BCd28
#for (n in numaccesNCBIBCd28) {
#  n3 <- c(paste0("fastq-dump --split-3 ", n))
#}
#func_y_doc_BLd28
#for (n in numaccesNCBIBLd28) {
#  n4 <- c(paste0("fastq-dump --split-3 ", n))
#}
#print(n3) #BCd28
#print(n4) #BLd28

#Posibilidades para enriquecer el documento
#Loop through SRA accession numbers ##No se automatizó/no se hizo
#Execute fastq-dump command (assuming system calls are allowed)
#system(fastq_dump)

