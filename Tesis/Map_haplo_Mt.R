setwd("C:/Users/IquiB/OneDrive/Escritorio/Iqui/CCG/Análisis de haplogrupos/Mapa Pakistán")
library(tidyverse)
library(dplyr)
library(ggplot2)
library("mxmaps")
library(raster)
library(broom)
library(reshape2)
library(viridis)
library(rgdal)
library(RColorBrewer)
library(ggnewscale)
library(scatterpie)
library(matrixStats)

######
#Obtención de los datos de Punjabi

hapPJ <- read.csv("1KG_para_PJL.csv")
hapPJ2 <- read.csv("Punjabi.csv", header = F)

hapPJ$SampleID[hapPJ$SampleID %in% hapPJ2$V1]

P1 <- hapPJ[which(hapPJ$SampleID==1583),]
P2 <- hapPJ[which(hapPJ$SampleID==1586),]
P3 <- hapPJ[which(hapPJ$SampleID==1589),]
P4 <- hapPJ[which(hapPJ$SampleID==2490),]
P5 <- hapPJ[which(hapPJ$SampleID==2493),]
P6 <- hapPJ[which(hapPJ$SampleID==2597),]
P7 <- hapPJ[which(hapPJ$SampleID==2600),]
P8 <- hapPJ[which(hapPJ$SampleID==2603),]
P9 <- hapPJ[which(hapPJ$SampleID==2648),]
P10 <- hapPJ[which(hapPJ$SampleID==2651),]
P11 <- hapPJ[which(hapPJ$SampleID==2654),]
P12 <- hapPJ[which(hapPJ$SampleID==2657),]
P13 <- hapPJ[which(hapPJ$SampleID==2660),]
P14 <- hapPJ[which(hapPJ$SampleID==2681),]
P15 <- hapPJ[which(hapPJ$SampleID==2684),]
P16 <- hapPJ[which(hapPJ$SampleID==2687),]
P17 <- hapPJ[which(hapPJ$SampleID==2690),]
P18 <- hapPJ[which(hapPJ$SampleID==2696),]
P19 <- hapPJ[which(hapPJ$SampleID==2724),]
P20 <- hapPJ[which(hapPJ$SampleID==2727),]
P21 <- hapPJ[which(hapPJ$SampleID==2733),]
P22 <- hapPJ[which(hapPJ$SampleID==2736),]
P23<- hapPJ[which(hapPJ$SampleID==2774),]
P24 <- hapPJ[which(hapPJ$SampleID==2780),]
P25<- hapPJ[which(hapPJ$SampleID==2783),]
P26<- hapPJ[which(hapPJ$SampleID==2786),]
P27<- hapPJ[which(hapPJ$SampleID==2789),]
P28<- hapPJ[which(hapPJ$SampleID==2792),]
P29<- hapPJ[which(hapPJ$SampleID==3015),]
P30<- hapPJ[which(hapPJ$SampleID==3018),]
P31 <- hapPJ[which(hapPJ$SampleID==2699),]
P32 <- hapPJ[which(hapPJ$SampleID==3021),]
P33<- hapPJ[which(hapPJ$SampleID==3228),]
P34<- hapPJ[which(hapPJ$SampleID==3234),]
P35<- hapPJ[which(hapPJ$SampleID==3237),]
P36<- hapPJ[which(hapPJ$SampleID==3490),]
P37<- hapPJ[which(hapPJ$SampleID==3624),]
P38 <- hapPJ[which(hapPJ$SampleID==3629),]
P39 <- hapPJ[which(hapPJ$SampleID==3636),]
P40<- hapPJ[which(hapPJ$SampleID==3649),]
P41<- hapPJ[which(hapPJ$SampleID==3652),]
P42<- hapPJ[which(hapPJ$SampleID==3660),]
P43<- hapPJ[which(hapPJ$SampleID==3663),]
P44<- hapPJ[which(hapPJ$SampleID==3667),]
P45<- hapPJ[which(hapPJ$SampleID==3702),]
P46<- hapPJ[which(hapPJ$SampleID==3705),]
P47 <- hapPJ[which(hapPJ$SampleID==3708),]
P48 <- hapPJ[which(hapPJ$SampleID==3767),]

popP<-rbind(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,
           P21,P22,P23,P24,P25,P26,P27,P28,P29,P30,P31,P32,P33,P34,P35,P36,P37,P38,P39,P40,
           P41,P42,P43,P44,P45,P46,P47,P48)

popP$Population<-"Punjabi"
popP$Region<-""

MtHC <- read.csv("hoja completa2 de HMt.csv")

MtHC$Haplogroup <- MtHC$MtDNA.haplogroup
MtHC$SampleID <- MtHC$Sample.ID
MtHC$Region <- MtHC$Regional.group

Mt<-MtHC[,c("SampleID","Haplogroup","Population","Region")]

MtF<-rbind(Mt,popP)

write.table(MtF, file="haploMt_resumen.txt",quote=FALSE, row.names=FALSE,sep=",")

######
#Creación del data frame con los datos dela ascendencia materna

MtF <- read.csv("hoja completa2 de HMt.csv")

pakHazM <- MtF%>%
  filter(Population == "Hazara")
pakPunM <-MtF %>%
  filter(Population == "Punjabi")
pakKalM <-MtF[which(MtF[,3]=="Kalash"),]
pakMakM <-MtF[which(MtF[,3]=="Makrani"),]
pakBraM <-MtF[which(MtF[,3]=="Brahui"),]
pakPatM <-MtF[which(MtF[,3]=="Pathan"),]
pakBurM <-MtF[which(MtF[,3]=="Burusho"),]
pakBalM <-MtF[which(MtF[,3]=="Balochi"),]
pakSinM <-MtF[which(MtF[,3]=="Sindhi"),]

popPMt<-rbind(pakSinM,pakBalM,pakBurM,pakPatM,pakBraM,pakMakM,pakKalM,pakHazM,pakPunM)

write.table(popPMt, file="hapMt_Pak.txt",quote=FALSE, row.names=FALSE,sep=",")
######
#Creación del mapa de Pakistán con gráficas de haplogrupos mitocondriales con un carácter

popPMt<-read.csv("hapMt_Pak.txt")

f4<-popPMt[,c("Population","Haplogroup")]

M1<-strsplit(x=f4$Haplogroup,split = "")
M2<-lapply(M1, '[', 1)
M3<-lapply(M2, FUN=paste, collapse="")
f4$Haplogroup1<-unlist(M3)

f_summary4 <- dcast(melt(f4),  Population ~ Haplogroup1)
ncol(f_summary4)
f_summary4$samples <-  apply(f_summary4[,2:19], 1, sum)
f_summary4[,2:19] <- f_summary4[,2:19]/f_summary4[,"samples"]

shapefile_Pa <- readOGR(dsn = "PAK_adm1.shp")

data_Pa <- tidy(shapefile_Pa)

data_Pa$value4 <- NA
get_hap <- function(hap)
{
  for (comunidad in f_summary4$Population){
    data_Pa[data_Pa$Population==comunidad,]$value4 <- f_summary4[f_summary4$f==comunidad,hap]
    
  }
  return(data_Pa)
}

f_summary4$long <- NA
f_summary4$lat <- NA

f_summary4[f_summary4$Population=="Kalash", ]$lat <- 35.73
f_summary4[f_summary4$Population=="Kalash", ]$long <- 71.76
f_summary4[f_summary4$Population=="Brahui", ]$lat <- 27.77
f_summary4[f_summary4$Population=="Brahui", ]$long <- 66.61
f_summary4[f_summary4$Population=="Pathan", ]$lat <- 30.95
f_summary4[f_summary4$Population=="Pathan", ]$long <- 70.95
f_summary4[f_summary4$Population=="Burusho", ]$lat <- 36.52
f_summary4[f_summary4$Population=="Burusho", ]$long <- 73.85
f_summary4[f_summary4$Population=="Balochi", ]$lat <- 27.71
f_summary4[f_summary4$Population=="Balochi", ]$long <- 64.81
f_summary4[f_summary4$Population=="Sindhi", ]$lat <- 26.04
f_summary4[f_summary4$Population=="Sindhi", ]$long <- 68.97
f_summary4[f_summary4$Population=="Makrani", ]$lat <- 24.60
f_summary4[f_summary4$Population=="Makrani", ]$long <- 68.07
f_summary4[f_summary4$Population=="Hazara", ]$lat <- 34.94
f_summary4[f_summary4$Population=="Hazara", ]$long <- 72.41
f_summary4[f_summary4$Population=="Punjabi", ]$lat <- 30
f_summary4[f_summary4$Population=="Punjabi", ]$long <- 69.5

f_summary4$Population <- as.factor(f_summary4$Population)
f_summary4$long <- as.numeric(f_summary4$long)
f_summary4$lat <- as.numeric(f_summary4$lat)

f_summary4$samples2 <-  apply(f_summary4[,2:19], 1, sum)

unique(f4$Haplogroup1)

f_summary4$radio<- (f_summary4$samples)/40

Pn4 <- ggplot(data_Pa)  +
  labs(title = "Haplogroups of indigenous populations from Pakistan") +
  geom_polygon(aes( x= long, y = lat, group = group),
               fill = "white", color="black" , linewidth=0.5) +
  theme_void() +
  theme(panel.background = element_rect(size= 0.5, color = "white", fill = "white"))

sp_piecharts4 <- Pn4 + geom_scatterpie(aes(x=long, y=lat, group = Population, r=radio),
                                       data=f_summary4, cols=colnames(f_summary4)[2:19], alpha=0.8) + 
  labs(fill="Haplogroup", title = "Haplogroups Mt of indigenous populations from Pakistan") +
  theme(plot.title = element_text(family = "serif",
                                  face = "bold",
                                  size = 25,
                                  hjust = .8,
                                  vjust = 1)) +
  scale_fill_manual(values = c("H"="#1E90FF", "M"="#104E8B", "T"="#E066FF", "J"="#7A378B", "L"="#FF3030", 
                               "W"="#8B1A1A", "U"="#FFA500", "K"="#8B5A00", "R"="#C0FF3E", "N"="#698B22", 
                               "D"="#F0FFF0", "I"="#FFFF00", "X"="#CDCD00", "C"="#FFE4C4", "Z"="blue", "F"="black",
                               "B"="red", "G"="green"))



