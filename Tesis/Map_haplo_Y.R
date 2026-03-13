setwd("C:/Users/IquiB/OneDrive/Escritorio/Iqui/CCG/Análisis de haplogrupos/Mapa Pakistán")
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

######
#Creación del data frame con los datos dela ascendencia paterna

haploy2 <- read.csv("Haplogrupos Y P.csv")

pakHaz <-haploy2 %>%
  filter(Population == "Hazara")
pakPun <-haploy2 %>%
  filter(Population == "Punjabi")
pakKal <-haploy2[which(haploy2[,3]=="Kalash"),]
pakMak <-haploy2[which(haploy2[,3]=="Makrani"),]
pakBra <-haploy2[which(haploy2[,3]=="Brahui"),]
pakPat <-haploy2[which(haploy2[,3]=="Pathan"),]
pakBur <-haploy2[which(haploy2[,3]=="Burusho"),]
pakBal <-haploy2[which(haploy2[,3]=="Balochi"),]
pakSin <-haploy2[which(haploy2[,3]=="Sindhi"),]

populations<-rbind(pakSin,pakBal,pakBur,pakPat,pakBra,pakMak,pakKal,pakHaz,pakPun)

write.table(f2, file="hapY_Pak.txt",quote=FALSE, row.names=FALSE,sep=",")

######
#Creación del mapa de Pakistán con gráficas de haplogrupos del cromosoma Y con un carácter

shapefile_Pa <- readOGR(dsn = "PAK_adm1.shp")
data_Pa <- tidy(shapefile_Pa)
f3<-populations[,c("Population","Haplogroup")]

HY1<-strsplit(x=f3$Haplogroup,split = "")
HY3<-lapply(HY1, '[', 1)
HYF2<-lapply(HY3, FUN=paste, collapse="")
f3$Haplogroup3<-unlist(HYF2)

f_summary3 <- dcast(melt(f3),  Population ~ Haplogroup3)
ncol(f_summary3)
f_summary3$samples <-  apply(f_summary3[,2:11], 1, sum)
f_summary3[,2:11] <- f_summary3[,2:11]/f_summary3[,"samples"]

data_Pa$value3 <- NA
get_hap <- function(hap)
{
  for (comunidad in f_summary3$Population){
    data_Pa[data_Pa$Population==comunidad,]$value3 <- f_summary3[f_summary3$f==comunidad,hap]
    
  }
  return(data_Pa)
}

f_summary3$long <- NA
f_summary3$lat <- NA

f_summary3[f_summary3$Population=="Kalash", ]$lat <- 35.73
f_summary3[f_summary3$Population=="Kalash", ]$long <- 71.76
f_summary3[f_summary3$Population=="Brahui", ]$lat <- 27.5
f_summary3[f_summary3$Population=="Brahui", ]$long <- 66.61
f_summary3[f_summary3$Population=="Pathan", ]$lat <- 30.95
f_summary3[f_summary3$Population=="Pathan", ]$long <- 70.5
f_summary3[f_summary3$Population=="Burusho", ]$lat <- 36.52
f_summary3[f_summary3$Population=="Burusho", ]$long <- 74
f_summary3[f_summary3$Population=="Balochi", ]$lat <- 28
f_summary3[f_summary3$Population=="Balochi", ]$long <- 64.81
f_summary3[f_summary3$Population=="Sindhi", ]$lat <- 27
f_summary3[f_summary3$Population=="Sindhi", ]$long <- 68.5
f_summary3[f_summary3$Population=="Makrani", ]$lat <- 25
f_summary3[f_summary3$Population=="Makrani", ]$long <- 68.07
f_summary3[f_summary3$Population=="Hazara", ]$lat <- 34.94
f_summary3[f_summary3$Population=="Hazara", ]$long <- 73
f_summary3[f_summary3$Population=="Punjabi", ]$lat <- 31
f_summary3[f_summary3$Population=="Punjabi", ]$long <- 72.5

f_summary3$Population <- as.factor(f_summary3$Population)
f_summary3$long <- as.numeric(f_summary3$long)
f_summary3$lat <- as.numeric(f_summary3$lat)

f_summary3$samples2 <-  apply(f_summary3[,2:11], 1, sum)

write.table(f_summary3, file="haploY_para_mapa.txt",quote=FALSE, row.names=FALSE,sep=",")

f_summary3$radio<- as.numeric((f_summary3$samples)/40)


Pn3 <- ggplot(data_Pa)  +
  labs(title = "Haplogroups of indigenous populations from Pakistan") +
  geom_polygon(aes( x= long, y = lat, group = group),
               fill = "white", color="black" , linewidth=0.5) +
  theme_void() +
  theme(panel.background = element_rect(size= 0.5, color = "white", fill = "white"))

sp_piecharts3 <- Pn3 + geom_scatterpie(aes(x=long, y=lat, group = Population, r=radio),
                                       data=f_summary3, cols=colnames(f_summary3)[2:11], alpha=0.8) +
  labs(fill="Haplogrupo", title = "Haplogrupos del cromosoma Y en grupos culturales de Pakistán") +
  theme(legend.text = element_text(size = 12),legend.title =element_text(size=12),
        plot.title.position="plot",plot.title = element_text(family = "serif",
                                                             face = "bold",
                                                             size = 19,
                                                             hjust = .8,
                                                             vjust = 1)) +
  scale_fill_brewer(palette="Paired") + 
  geom_scatterpie_legend(f_summary3$radio, x=65, y=34, n=3, labeller = function(x) 40*x)
sp_piecharts3


######
##Creación del mapa de Pakistán con gráficas de haplogrupos del cromosoma Y con dos carácteres

shapefile_Pa <- readOGR(dsn = "PAK_adm1.shp")
data_Pa <- tidy(shapefile_Pa)

f2<-read.csv("hapY_Pak.txt")


HY<-strsplit(x=f2$Haplogroup,split = "")
HY2<-lapply(HY, '[', 1:2)
HYF<-lapply(HY2, FUN=paste, collapse="")
f2$Haplogroup2<-unlist(HYF)

is(f2$Haplogroup)
unique(f2$Haplogroup2)

f_summary2 <- dcast(melt(f2),  Population ~ Haplogroup2)
ncol(f_summary2)
f_summary2$samples <-  apply(f_summary2[,2:17], 1, sum)
f_summary2[,2:17] <- f_summary2[,2:17]/f_summary2[,"samples"]

data_Pa$value2 <- NA
get_hap <- function(hap)
{
  for (comunidad in f_summary2$Population){
    data_Pa[data_Pa$Population==comunidad,]$value2 <- f_summary2[f_summary2$f==comunidad,hap]
    
  }
  return(data_Pa)
}

f_summary2$long <- NA
f_summary2$lat <- NA

f_summary2[f_summary2$Population=="Kalash", ]$lat <- 35.73
f_summary2[f_summary2$Population=="Kalash", ]$long <- 71.76
f_summary2[f_summary2$Population=="Brahui", ]$lat <- 27.77
f_summary2[f_summary2$Population=="Brahui", ]$long <- 66.61
f_summary2[f_summary2$Population=="Pathan", ]$lat <- 30.95
f_summary2[f_summary2$Population=="Pathan", ]$long <- 70.95
f_summary2[f_summary2$Population=="Burusho", ]$lat <- 36.52
f_summary2[f_summary2$Population=="Burusho", ]$long <- 73.85
f_summary2[f_summary2$Population=="Balochi", ]$lat <- 27.71
f_summary2[f_summary2$Population=="Balochi", ]$long <- 64.81
f_summary2[f_summary2$Population=="Sindhi", ]$lat <- 26.04
f_summary2[f_summary2$Population=="Sindhi", ]$long <- 68.97
f_summary2[f_summary2$Population=="Makrani", ]$lat <- 24.60
f_summary2[f_summary2$Population=="Makrani", ]$long <- 68.07
f_summary2[f_summary2$Population=="Hazara", ]$lat <- 34.94
f_summary2[f_summary2$Population=="Hazara", ]$long <- 72.41
f_summary2[f_summary2$Population=="Punjabi", ]$lat <- 30
f_summary2[f_summary2$Population=="Punjabi", ]$long <- 69.5

f_summary2$Population <- as.factor(f_summary2$Population)
f_summary2$long <- as.numeric(f_summary2$long)
f_summary2$lat <- as.numeric(f_summary2$lat)

f_summary2$samples2 <-  apply(f_summary2[,2:17], 1, sum)

f_summary2$radio<- (f_summary2$samples)/40

Pn2 <- ggplot(data_Pa)  +
  labs(title = "Haplogroups of indigenous populations from Pakistan") +
  geom_polygon(aes( x= long, y = lat, group = group),
               fill = "white", color="black" , size=0.5) +
  theme_void() +
  theme(panel.background = element_rect(size= 0.5, color = "white", fill = "white"))

sp_piecharts <- Pn2 + geom_scatterpie(aes(x=long, y=lat, group = Population, r=radio),
                                     data=f_summary2, cols=colnames(f_summary2)[2:17], alpha=0.8) + 
  labs(fill="Haplogroup", title = "Haplogroups Y of indigenous populations from Pakistan") +
  theme(plot.title = element_text(family = "serif",
                                  face = "bold",
                                  size = 25,
                                  hjust = .8,
                                  vjust = 1)) +
  scale_fill_manual(values = c("R2"="#1E90FF", "R1"="#104E8B", "J2"="#E066FF", "J1"="#7A378B", "H1"="#FF3030", 
                               "H0"="#8B1A1A", "G2"="#FFA500", "G1"="#8B5A00", "C2"="#C0FF3E", "C1"="#698B22", 
                               "C5"="#F0FFF0", "I2"="#FFFF00", "Q1"="#CDCD00", "L1"="#FFE4C4", "E1"="blue", "O2"="black"))

