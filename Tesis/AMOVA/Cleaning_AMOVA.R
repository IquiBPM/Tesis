library(readr)

SCPY<-readr::read_tsv("C:/Users/user/Downloads/AMOVA.tsv")

######
#Normalizar los ID

SCPY2<-SCPY

S<-SCPY[c(1:167),1]
SID<-data.frame(S,
                "n"=c(1:167))
ID<-strsplit(x=SID$Sample_id,split = "")
ID2<-lapply(ID, '[', c(5,6,7,8,9))
ID3<-lapply(ID2, FUN=paste, collapse="")
B1<-unlist(ID3)

S2<-SCPY[c(168:265),1]
SID2<-data.frame(S2,
                "n"=c(168:265))
IID<-strsplit(x=SID2$Sample_id,split = "")
IID2<-lapply(IID, '[', c(3,4,5,6,7))
IID3<-lapply(IID2, FUN=paste, collapse="")
B2<-unlist(IID3)

B3<-c(B1,B2)

SCPY$Sample_id<-unlist(B3)

AMOVA<-SCPY[,-2]
AMOVA<-AMOVA[,-2]

Labeles<-readr::read_tsv("C:/Users/user/Downloads/Labeles_B.txt")

new_file<-merge(y=AMOVA,x=Labeles,by="Sample_id") 

write.table(new_file,file = "Secuencias_Pakistan_AMOVA.tsv",sep = "\t",row.names = FALSE,col.names = TRUE)

new_file2[c(1,2,3),c(1,2,3,4,5)]
AMOVA_B[c(1,2,3),c(1,2,3,4)]

II<-AMOVA_B[c(1,2,3),c(1,2,3,4,5)]

AMOVA_B[which(AMOVA_B[,2]=="Burusho"),3]<-c(rep("Gilgit-Baltistan",19))
AMOVA_B[which(AMOVA_B[,2]=="Kalash"),3]<-c(rep("Khyber Pakhtunkhwa",16))

######
#Datos que no se encuentran en la tabla para el AMOVA

setwd("C:/Users/user/Downloads")
AMOVA_Y<-read.csv("C:/Users/user/Downloads/metadatos_AMOVA_Y.csv")
Est_desc<-read.table("tablaY_regiones.txt", sep = "\t", skip = 1)

Eliminados<-Est_desc[c(which(Est_desc$V1 %in% AMOVA_Y$Sample_id == FALSE)),]

nombre<-c("ID","Grupo","Haplogrupo","Haplogrupo (simplificado)", "Origen")
Elimina2<-rbind(nombre, Eliminados)

colnames(Eliminados) <-nombre

write.csv(Eliminados, file = "datos_eliminados", row.names = FALSE)
