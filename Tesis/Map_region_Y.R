setwd("C:/Users/user/Downloads")
library(tidyverse)
library(dplyr)
library(ggplot2)
library(raster)
library(broom)
library(reshape2)
library(viridis)
library(rgdal)
library(RColorBrewer)
library(ggnewscale)
library(scatterpie)
library(matrixStats)
library(paletteer)

haploy<-read.csv("hapY_Pak.txt")

######
#Clasificación de los haplogrupos (tres carácteres) por región

unique(haploy$Haplogroup)
M1<-strsplit(x=haploy$Haplogroup,split = "")
M1<-lapply(M1, '[', c(1,2,3))
M1<-lapply(M1, FUN=na.omit)
M1<-lapply(M1, FUN=paste, collapse="")
haploy$Haplogroup1<-unlist(M1)
HY<-unique(haploy$Haplogroup1)

AF<-rbind(haploy[which(haploy$Haplogroup1=="E1b"),])
AF$Region<-"África"

MI<-rbind(haploy[which(haploy$Haplogroup1=="G2a"),],haploy[which(haploy$Haplogroup1=="J1a"),],
          haploy[which(haploy$Haplogroup1=="L1a"),])
MI$Region<-"Medio Oriente"

EU<-rbind(haploy[which(haploy$Haplogroup1=="I2a"),],
          haploy[which(haploy$Haplogroup1=="R1b"),])
EU$Region<-"Europa"

EA<-rbind(haploy[which(haploy$Haplogroup1=="C2b"),],haploy[which(haploy$Haplogroup1=="C2e"),],
          haploy[which(haploy$Haplogroup1=="O2a"),])
EA$Region<-"Asia Oriental"

CS<-rbind(haploy[which(haploy$Haplogroup1=="C5"),],haploy[which(haploy$Haplogroup1=="G1a"),],
          haploy[which(haploy$Haplogroup1=="G2b"),],haploy[which(haploy$Haplogroup1=="H0"),],
          haploy[which(haploy$Haplogroup1=="H1"),],haploy[which(haploy$Haplogroup1=="H1a"),],
          haploy[which(haploy$Haplogroup1=="H1b"),],haploy[which(haploy$Haplogroup1=="J1"),],
          haploy[which(haploy$Haplogroup1=="J2a"),],haploy[which(haploy$Haplogroup1=="J2b"),],
          haploy[which(haploy$Haplogroup1=="L1"),],haploy[which(haploy$Haplogroup1=="L1b"),],
          haploy[which(haploy$Haplogroup1=="Q1b"),],haploy[which(haploy$Haplogroup1=="R1a"),],
          haploy[which(haploy$Haplogroup1=="R2a"),],haploy[which(haploy$Haplogroup1=="Q1a"),],
          haploy[which(haploy$Haplogroup1=="C1b"),])
CS$Region<-"Asia Central y del Sur"

world<-rbind(AF,CS,EA,EU,MI)
unique(world$Haplogroup1)
#write.csv(world, file = "tablaY_anexo", col.names = TRUE, sep = ",",row.names = FALSE)

######
#Tablas de las regiones de origen

f_summary<- dcast(melt(world),  Population ~ Region)
ncol(f_summary)
f_summary$Samples <-  apply(f_summary[,2:6], 1, sum)

f_summary2<- dcast(melt(world), Region ~ Population)
ncol(f_summary2)
f_summary2$Samples <-  apply(f_summary2[,2:10], 1, sum)
#write.table(f_summary2,file="HaploY_region_origin3(la_buenax3).csv",quote=FALSE,col.name=TRUE,row.names=FALSE,sep=",")

######
#Graficación de los porcentajes de las regiones sobre el mapa de Pakistán

f_summary[,2:6] <- f_summary[,2:6]/f_summary[,"Samples"]
setwd("C:/Users/IquiB/OneDrive/Escritorio/Iqui/CCG/Análisis de haplogrupos/Mapa Pakistán")
shapefile_Pa <- readOGR(dsn = "PAK_adm1.shp")
data_Pa <- tidy(shapefile_Pa)

data_Pa$value <- NA
get_hap <- function(hap)
{
  for (comunidad in f_summary$Population){
    data_Pa[data_Pa$Population==comunidad,]$value <- f_summary[f_summary$f==comunidad,hap]
    
  }
  return(data_Pa)
}

f_summary$long <- NA
f_summary$lat <- NA
f_summary[f_summary$Population=="Kalash", ]$lat <- 35.73
f_summary[f_summary$Population=="Kalash", ]$long <- 71.76
f_summary[f_summary$Population=="Brahui", ]$lat <- 27.5
f_summary[f_summary$Population=="Brahui", ]$long <- 66.61
f_summary[f_summary$Population=="Pathan", ]$lat <- 30.95
f_summary[f_summary$Population=="Pathan", ]$long <- 70.5
f_summary[f_summary$Population=="Burusho", ]$lat <- 36.52
f_summary[f_summary$Population=="Burusho", ]$long <- 74
f_summary[f_summary$Population=="Balochi", ]$lat <- 28
f_summary[f_summary$Population=="Balochi", ]$long <- 64.81
f_summary[f_summary$Population=="Sindhi", ]$lat <- 27
f_summary[f_summary$Population=="Sindhi", ]$long <- 68.5
f_summary[f_summary$Population=="Makrani", ]$lat <- 25
f_summary[f_summary$Population=="Makrani", ]$long <- 68.07
f_summary[f_summary$Population=="Hazara", ]$lat <- 34.94
f_summary[f_summary$Population=="Hazara", ]$long <- 73
f_summary[f_summary$Population=="Punjabi", ]$lat <- 31
f_summary[f_summary$Population=="Punjabi", ]$long <- 72.5

f_summary$Population <- as.factor(f_summary$Population)
f_summary$long <- as.numeric(f_summary$long)
f_summary$lat <- as.numeric(f_summary$lat)

f_summary$samples2 <-  apply(f_summary[,2:6], 1, sum)

f_summary$radio<- (f_summary$Samples)/40

Pn <- ggplot(data_Pa)  +
  labs(title = "Region of origion of haplogrups Y from Pakistan") +
  geom_polygon(aes( x= long, y = lat, group = group),
               fill = "white", color="black" , size=0.5) +
  theme_void() +
  theme(panel.background = element_rect(size= 0.5, color = "white", fill = "white"))

sp_piecharts<- Pn + geom_scatterpie(data=f_summary, aes(x=long, y=lat, group = Population, r=radio),
                                    cols=colnames(f_summary)[2:6], alpha=1) +
  labs(fill="Región", title = "Región de origen de los haplogrupos del cromosoma Y") +
  theme(legend.text = element_text(size = 12),legend.position=c(.85,.25),legend.title = element_text(size=12),
        plot.title.position="plot",plot.title = element_text(family = "serif",
                                                             face = "bold",
                                                             size = 20,
                                                             hjust = .8,
                                                             vjust = 1)) +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C")) +
  geom_scatterpie_legend(f_summary$radio, x=65, y=34, n=3, labeller = function(x) 40*x)

######
#Deducción de las regiones de origen de los haplogrupos presentes en Pakistán
#con base a su distrubución mundial

YF <- read.csv("Haplo_Y_World.csv")

M3<-strsplit(x=YF$Haplogroup,split = "")
M3<-lapply(M3, '[', c(1,2,3))
M3<-lapply(M3, FUN=na.omit)
M3<-lapply(M3, FUN=paste, collapse="")
YF$Haplogroup1<-unlist(M3)

HP<-YF %>%
  filter(Population == "Punjabi" | Population == "Brahui" | Population == "Kalash"|
           Population == "Makrani" | Population == "Hazara" | Population == "Pathan"|
           Population == "Burusho" | Population == "Balochi" | Population == "Sindhi")
haplo<-unique(HP$Haplogroup1)

haplo[order(haplo)]

YF2<-YF
YFHP<-rbind(YF2[which(YF2$Haplogroup1=="C1b"),],YF2[which(YF2$Haplogroup1=="C2b"),],
            YF2[which(YF2$Haplogroup1=="C2e"),],YF2[which(YF2$Haplogroup1=="C5"),],
            YF2[which(YF2$Haplogroup1=="E1b"),],YF2[which(YF2$Haplogroup1=="G1a"),],
            YF2[which(YF2$Haplogroup1=="G2a"),],YF2[which(YF2$Haplogroup1=="G2b"),],
            YF2[which(YF2$Haplogroup1=="H0"),],YF2[which(YF2$Haplogroup1=="H1"),],
            YF2[which(YF2$Haplogroup1=="H1a"),],YF2[which(YF2$Haplogroup1=="H1b"),],
            YF2[which(YF2$Haplogroup1=="I2a"),])

unique(YFHP$Region)

#
CSA<-YFHP %>%
  filter(Population == "Hazara" | Region == "CENTRAL_SOUTH_ASIA")

unique(CSA$Haplogroup1)
NCSA <- c(count(CSA$Haplogroup1=="C1b"),count(CSA$Haplogroup1=="C2b"),
          count(CSA$Haplogroup1=="C2e"),count(CSA$Haplogroup1=="C5"),
          count(CSA$Haplogroup1=="E1b"),count(CSA$Haplogroup1=="G1a"),
          count(CSA$Haplogroup1=="G2a"),
          count(CSA$Haplogroup1=="G2b"),count(CSA$Haplogroup1=="H0"),
          count(CSA$Haplogroup1=="H1"),
          count(CSA$Haplogroup1=="H1a"),count(CSA$Haplogroup1=="H1b"),
          count(CSA$Haplogroup1=="I2a"))
sum(NCSA)
CSA<-data.frame(
  "Haplogroup"=c("C1b","C2b","C2e","C5","E1b","G1a","G2a","G2b","H0",
                 "H1","H1a","H1b","I2a"),
  "ab_rel"=(NCSA/sum(NCSA))*100,
  "Absolute_Abundance"=NCSA,
  "Region"= "Asia Central y del Sur")

#
EU<-YFHP %>%
  filter(Region == "EUROPE")

unique(EU$Haplogroup1)
NEU <- c(count(EU$Haplogroup1=="E1b"),count(EU$Haplogroup1=="G2a"),
         count(EU$Haplogroup1=="I2a"))

sum(NEU)
EU<-data.frame(
  "Haplogroup"=c("E1b","G2a","I2a"),
  "ab_rel"=(NEU/sum(NEU))*100,
  "Absolute_Abundance"=NEU,
  "Region"= "Europa")

#
ME<-YFHP %>%
  filter(Region == "MIDDLE_EAST")

unique(ME$Haplogroup1)
NME <- c(count(ME$Haplogroup1=="E1b"),count(ME$Haplogroup1=="G2a"))

sum(NME)
ME<-data.frame(
  "Haplogroup"=c("E1b","G2a"),
  "ab_rel"=(NME/sum(NME))*100,
  "Absolute_Abundance"=NME,
  "Region"= "Medio Oriente")

#
EA<-YFHP %>%
  filter(Region == "EAST_ASIA" | Region == "mixed" & Population != "Hazara")

unique(EA$Haplogroup1)
NEA <- c(count(EA$Haplogroup1=="C1b"),count(EA$Haplogroup1=="C2b"),
         count(EA$Haplogroup1=="C2e"))

sum(NEA)
EA<-data.frame(
  "Haplogroup"=c("C1b","C2b","C2e"),
  "ab_rel"=(NEA/sum(NEA))*100,
  "Absolute_Abundance"=NEA,
  "Region"= "Asia Oriental")

#
AM<-YFHP %>%
  filter(Region == "AMERICA")

unique(AM$Haplogroup1)
NAM <- c(count(AM$Haplogroup1=="G2a"))

sum(NAM)
AM<-data.frame(
  "Haplogroup"=c("G2a"),
  "ab_rel"=(NAM/sum(NAM))*100,
  "Absolute_Abundance"=NAM,
  "Region"= "América")

#
AF<-YFHP %>%
  filter(Region == "AFRICA")

unique(AF$Haplogroup1)
NAF <- count(AF$Haplogroup1=="E1b")

sum(NAF)
AF<-data.frame(
  "Haplogroup"=c("E1b"),
  "ab_rel"=(NAF/sum(NAF))*100,
  "Absolute_Abundance"=NAF,
  "Region"= "África")

#
OC<-YFHP %>%
  filter(Region == "OCEANIA")

unique(OC$Haplogroup1)
NOC <-count(OC$Haplogroup1=="C1b")

sum(NOC)
OC<-data.frame(
  "Haplogroup"=c("C1b"),
  "ab_rel"=(NOC/sum(NOC))*100,
  "Absolute_Abundance"=NOC,
  "Region"= "Oceania")

######
#Graficación de la distribución mundial de los haplogrupos presentes en Pakistán

unique(YFHP$Region)
YFHP2<-rbind(EU,ME,CSA,EA,AM,AF,OC)
unique(YFHP2$Haplogroup)

colnames(YFHP2)<-c("Haplogrupo","Abundancia Relativa","Abundancia_Absoluta","Región")

ggplot(YFHP2, aes(fill=Haplogrupo, x=Región, y=Abundancia_Absoluta)) + 
  labs(title = "Abundancia de los haplogrupos del cromosoma Y en el mundo") +
  geom_bar(stack="dodge", stat="identity") +
  scale_fill_manual(values=paletteer_d("ggsci::default_igv")) +
  theme(axis.text.x = element_text(size = 8, hjust = 1, angle = 40),
        plot.title.position="plot",plot.title = element_text(family = "serif",
                                                             face = "bold",
                                                             size = 15,
                                                             hjust = .8,
                                                             vjust = 1))+
  scale_y_continuous(breaks = seq(0,60,10))