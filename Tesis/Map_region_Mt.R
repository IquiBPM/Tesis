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

haplom<-read.csv("hapMt_Pak.txt")

######
#Clasificación de los haplogrupos (tres carácteres) por región

unique(haplom$Haplogroup)
Mt1<-strsplit(x=haplom$Haplogroup,split = "")
Mt1<-lapply(Mt1, '[', c(1,2,3))
Mt1[[61]]<-gsub("[+]",NA,Mt1[[61]])
Mt1[[71]]<-gsub("[+]",NA,Mt1[[71]])
Mt1[[126]]<-gsub("[+]",NA,Mt1[[126]])
Mt1[[161]]<-gsub("[+]",NA,Mt1[[161]])
Mt1[[163]]<-gsub("[+]",NA,Mt1[[163]])
Mt1<-lapply(Mt1, FUN=na.omit)
Mt1<-lapply(Mt1, FUN=paste, collapse="")
haplom$Haplogroup1<-unlist(Mt1)
HM<-unique(haplom$Haplogroup1)

AF<-rbind(haplom[which(haplom$Haplogroup1=="L0a"),],haplom[which(haplom$Haplogroup1=="L1b"),],
          haplom[which(haplom$Haplogroup1=="L1c"),],haplom[which(haplom$Haplogroup1=="L2a"),],
          haplom[which(haplom$Haplogroup1=="L3d"),])
AF$Region<-"África"

EU<-rbind(haplom[which(haplom$Haplogroup1=="T1a"),],haplom[which(haplom$Haplogroup1=="T2b"),],
            haplom[which(haplom$Haplogroup1=="H5a"),],haplom[which(haplom$Haplogroup1=="U5a"),]
          )
EU$Region<-"Europa"

EA<-rbind(haplom[which(haplom$Haplogroup1=="B4c"),],haplom[which(haplom$Haplogroup1=="C4"),],
          haplom[which(haplom$Haplogroup1=="C4a"),],haplom[which(haplom$Haplogroup1=="D4j"),],
          haplom[which(haplom$Haplogroup1=="D5a"),],haplom[which(haplom$Haplogroup1=="F1"),],
          haplom[which(haplom$Haplogroup1=="F1a"),],haplom[which(haplom$Haplogroup1=="F1b"),],
          haplom[which(haplom$Haplogroup1=="F2"),],
          haplom[which(haplom$Haplogroup1=="M4"),],haplom[which(haplom$Haplogroup1=="M9a"),])
EA$Region<-"Asia Oriental"

CS<-rbind(haplom[which(haplom$Haplogroup1=="G2a"),],haplom[which(haplom$Haplogroup1=="G3b"),],
          haplom[which(haplom$Haplogroup1=="H1"),],
          haplom[which(haplom$Haplogroup1=="H13"),],
          haplom[which(haplom$Haplogroup1=="H14"),],haplom[which(haplom$Haplogroup1=="H2a"),],
          haplom[which(haplom$Haplogroup1=="H2b"),],haplom[which(haplom$Haplogroup1=="H4"),],
          haplom[which(haplom$Haplogroup1=="H6a"),],
          haplom[which(haplom$Haplogroup1=="HV"),],
          haplom[which(haplom$Haplogroup1=="HV2"),],haplom[which(haplom$Haplogroup1=="I1b"),],
          haplom[which(haplom$Haplogroup1=="I4"),],haplom[which(haplom$Haplogroup1=="J1b"),],
          haplom[which(haplom$Haplogroup1=="J1d"),],
          haplom[which(haplom$Haplogroup1=="J2b"),],haplom[which(haplom$Haplogroup1=="K2a"),],

          haplom[which(haplom$Haplogroup1=="M18"),],
          haplom[which(haplom$Haplogroup1=="M2a"),],haplom[which(haplom$Haplogroup1=="M33"),],
          haplom[which(haplom$Haplogroup1=="M38"),],haplom[which(haplom$Haplogroup1=="M39"),],
          haplom[which(haplom$Haplogroup1=="M3a"),],haplom[which(haplom$Haplogroup1=="M3c"),],
          haplom[which(haplom$Haplogroup1=="M52"),],
          haplom[which(haplom$Haplogroup1=="M5a"),],haplom[which(haplom$Haplogroup1=="M6"),],
          haplom[which(haplom$Haplogroup1=="M65"),],haplom[which(haplom$Haplogroup1=="M6a"),],
          haplom[which(haplom$Haplogroup1=="N1d"),],
          haplom[which(haplom$Haplogroup1=="N1e"),],
          haplom[which(haplom$Haplogroup1=="R"),],
          haplom[which(haplom$Haplogroup1=="R0a"),],haplom[which(haplom$Haplogroup1=="R2"),],
          haplom[which(haplom$Haplogroup1=="R5a"),],haplom[which(haplom$Haplogroup1=="R6a"),],
          haplom[which(haplom$Haplogroup1=="T2"),],
          haplom[which(haplom$Haplogroup1=="T2a"),],haplom[which(haplom$Haplogroup1=="T2c"),],
          haplom[which(haplom$Haplogroup1=="U2"),],haplom[which(haplom$Haplogroup1=="U2b"),],
          haplom[which(haplom$Haplogroup1=="U2a"),],
          haplom[which(haplom$Haplogroup1=="U2e"),],haplom[which(haplom$Haplogroup1=="U4a"),],
          haplom[which(haplom$Haplogroup1=="U4b"),],
          haplom[which(haplom$Haplogroup1=="U7"),],haplom[which(haplom$Haplogroup1=="U7a"),],
          haplom[which(haplom$Haplogroup1=="U7b"),],haplom[which(haplom$Haplogroup1=="U9b"),],
          haplom[which(haplom$Haplogroup1=="W"),],haplom[which(haplom$Haplogroup1=="W1c"),],
          haplom[which(haplom$Haplogroup1=="W3"),],haplom[which(haplom$Haplogroup1=="W3a"),],
          haplom[which(haplom$Haplogroup1=="W4"),],haplom[which(haplom$Haplogroup1=="X2"),],
          haplom[which(haplom$Haplogroup1=="X2p"),],haplom[which(haplom$Haplogroup1=="Z3"),])
CS$Region<-"Asia Central y del Sur"

ME<-rbind(haplom[which(haplom$Haplogroup1=="H"),],haplom[which(haplom$Haplogroup1=="H7b"),],
          haplom[which(haplom$Haplogroup1=="J2a"),],
          haplom[which(haplom$Haplogroup1=="M30"),],haplom[which(haplom$Haplogroup1=="U1a"),])
ME$Region<-"Medio Oriente"

world<-rbind(AF,CS,EA,EU,ME)
unique(world$Haplogroup1)
#write.csv(world, file = "tablaMt_anexo", col.names = TRUE, row.names = FALSE, sep = ",")

######
#Tablas de las regiones de origen
f_summary<- dcast(melt(world),  Population ~ Region)
ncol(f_summary)
f_summary$Samples <-  apply(f_summary[,2:6], 1, sum)

f_summary2<- dcast(melt(world), Region ~ Population)
ncol(f_summary2)
f_summary2$Samples <-  apply(f_summary2[,2:10], 1, sum)
##write.table(f_summary2,file="HaploMT_region_origin3(la_buena).csv",quote=FALSE,col.name=TRUE,row.names=FALSE,sep=",")

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
  labs(fill="Región", title = "Región de origen de los haplogrupos del DNA mitocondrial") +
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

MF <- read.csv("Haplo_Mt_World.txt")
M3<-strsplit(x=MF$Haplogroup,split = "")
M3<-lapply(M3, '[', c(1,2,3))
M3[[64]]<-gsub("[+]",NA,M3[[64]])
M3[[66]]<-gsub("[+]",NA,M3[[66]])
M3[[89]]<-gsub("[+]",NA,M3[[89]])
M3[[118]]<-gsub("[+]",NA,M3[[118]])
M3[[165]]<-gsub("[+]",NA,M3[[165]])
M3[[239]]<-gsub("[+]",NA,M3[[239]])
M3[[244]]<-gsub("[+]",NA,M3[[244]])
M3[[326]]<-gsub("[+]",NA,M3[[326]])
M3[[327]]<-gsub("[+]",NA,M3[[327]])
M3[[333]]<-gsub("[+]",NA,M3[[333]])
M3[[336]]<-gsub("[+]",NA,M3[[336]])
M3[[453]]<-gsub("[+]",NA,M3[[453]])
M3[[468]]<-gsub("[+]",NA,M3[[468]])
M3[[536]]<-gsub("[+]",NA,M3[[536]])
M3[[542]]<-gsub("[+]",NA,M3[[542]])
M3<-lapply(M3, FUN=na.omit)
M3<-lapply(M3, FUN=paste, collapse="")
MF$Haplogroup1<-unlist(M3)

#MF2<-MF[which(MF$Population!="Punjabi"),]
#MF2<-MF2[which(MF2$Population!="Hazara"),]
#MF2<-MF2[which(MF2$Population!="Kalash"),]
#MF2<-MF2[which(MF2$Population!="Makrani"),]
#MF2<-MF2[which(MF2$Population!="Brahui"),]
#MF2<-MF2[which(MF2$Population!="Pathan"),]
#MF2<-MF2[which(MF2$Population!="Burusho"),]
#MF2<-MF2[which(MF2$Population!="Balochi"),]
#MF2<-MF2[which(MF2$Population!="Sindhi"),]

HMW<-unique(MF2$Haplogroup1)
HMW

MFHP<-rbind(MF2[which(MF2$Haplogroup1=="H2a"),],MF2[which(MF2$Haplogroup1=="M30"),],
            MF2[which(MF2$Haplogroup1=="M2a"),],MF2[which(MF2$Haplogroup1=="T2b"),],
            MF2[which(MF2$Haplogroup1=="J1b"),],MF2[which(MF2$Haplogroup1=="L1c"),],
            MF2[which(MF2$Haplogroup1=="L0a"),],MF2[which(MF2$Haplogroup1=="M5a"),],
            MF2[which(MF2$Haplogroup1=="M18"),],MF2[which(MF2$Haplogroup1=="W3a"),],
            MF2[which(MF2$Haplogroup1=="M3a"),],MF2[which(MF2$Haplogroup1=="U2a"),],
            MF2[which(MF2$Haplogroup1=="W4"),],MF2[which(MF2$Haplogroup1=="K2a"),],
            MF2[which(MF2$Haplogroup1=="M65"),],MF2[which(MF2$Haplogroup1=="H13"),],
            MF2[which(MF2$Haplogroup1=="U4a"),],MF2[which(MF2$Haplogroup1=="R2"),],
            MF2[which(MF2$Haplogroup1=="N1e"),],MF2[which(MF2$Haplogroup1=="U2b"),],
            MF2[which(MF2$Haplogroup1=="N1d"),],MF2[which(MF2$Haplogroup1=="HV"),],
            MF2[which(MF2$Haplogroup1=="HV2"),],MF2[which(MF2$Haplogroup1=="U4b"),],
            MF2[which(MF2$Haplogroup1=="D5a"),],MF2[which(MF2$Haplogroup1=="W3"),],
            
            MF2[which(MF2$Haplogroup1=="I4"),],MF2[which(MF2$Haplogroup1=="U1a"),],
            MF2[which(MF2$Haplogroup1=="H"),],MF2[which(MF2$Haplogroup1=="R5a"),],
            MF2[which(MF2$Haplogroup1=="I1b"),],MF2[which(MF2$Haplogroup1=="X2"),],
            MF2[which(MF2$Haplogroup1=="C4"),],MF2[which(MF2$Haplogroup1=="T2c"),],
            MF2[which(MF2$Haplogroup1=="U9b"),],MF2[which(MF2$Haplogroup1=="J1d"),],
            MF2[which(MF2$Haplogroup1=="T1a"),],MF2[which(MF2$Haplogroup1=="H2b"),],
            MF2[which(MF2$Haplogroup1=="U2"),],MF2[which(MF2$Haplogroup1=="R0a"),],
            MF2[which(MF2$Haplogroup1=="U7"),],MF2[which(MF2$Haplogroup1=="J2b"),],
            MF2[which(MF2$Haplogroup1=="D4j"),],MF2[which(MF2$Haplogroup1=="M3c"),],
            MF2[which(MF2$Haplogroup1=="H1"),],MF2[which(MF2$Haplogroup1=="H4"),],
            MF2[which(MF2$Haplogroup1=="T2"),],MF2[which(MF2$Haplogroup1=="H6a"),],
            MF2[which(MF2$Haplogroup1=="L1b"),],MF2[which(MF2$Haplogroup1=="L2a"),],
            
            MF2[which(MF2$Haplogroup1=="L3d"),],MF2[which(MF2$Haplogroup1=="H14"),],
            MF2[which(MF2$Haplogroup1=="W"),],MF2[which(MF2$Haplogroup1=="T2a"),],
            MF2[which(MF2$Haplogroup1=="U2e"),],MF2[which(MF2$Haplogroup1=="U5a"),],
            MF2[which(MF2$Haplogroup1=="Z3"),],MF2[which(MF2$Haplogroup1=="M4"),],
            MF2[which(MF2$Haplogroup1=="M9a"),],MF2[which(MF2$Haplogroup1=="F2"),],
            MF2[which(MF2$Haplogroup1=="H7b"),],MF2[which(MF2$Haplogroup1=="B4c"),],
            MF2[which(MF2$Haplogroup1=="C4a"),],MF2[which(MF2$Haplogroup1=="H5a"),],
            MF2[which(MF2$Haplogroup1=="F1"),],MF2[which(MF2$Haplogroup1=="G2a"),],
            MF2[which(MF2$Haplogroup1=="M6a"),],MF2[which(MF2$Haplogroup1=="M39"),],
            MF2[which(MF2$Haplogroup1=="R6a"),],MF2[which(MF2$Haplogroup1=="U7a"),],
            MF2[which(MF2$Haplogroup1=="M33"),],MF2[which(MF2$Haplogroup1=="F1a"),],
            MF2[which(MF2$Haplogroup1=="U7b"),],MF2[which(MF2$Haplogroup1=="X2p"),],
            
            MF2[which(MF2$Haplogroup1=="G3b"),],MF2[which(MF2$Haplogroup1=="R"),],
            MF2[which(MF2$Haplogroup1=="W1c"),],MF2[which(MF2$Haplogroup1=="M52"),],
            MF2[which(MF2$Haplogroup1=="F1b"),],MF2[which(MF2$Haplogroup1=="M38"),],
            MF2[which(MF2$Haplogroup1=="M6"),])

unique(MFHP$Region)



#
EU<-MFHP %>%
  filter(Region == "Europe")

unique(EU$Haplogroup1)
NEU <- c(count(EU$Haplogroup1=="H2a"),count(EU$Haplogroup1=="U1a"),
         count(EU$Haplogroup1=="T2b"),count(EU$Haplogroup1=="H"),
         count(EU$Haplogroup1=="H13"),count(EU$Haplogroup1=="T1a"),
         count(EU$Haplogroup1=="U4a"),count(EU$Haplogroup1=="H2b"),
         count(EU$Haplogroup1=="J2b"),count(EU$Haplogroup1=="U5a"),
         count(EU$Haplogroup1=="H1"),count(EU$Haplogroup1=="H5a"),
         count(EU$Haplogroup1=="U2e"))

sum(NEU)
EU<-data.frame(
  "Haplogroup"=c("H2a","U1a","T2b","H","H13","T1a","U4a","H2b","J2b","U5a","H1","H5a","U2e"),
  "ab_rel"=(NEU/sum(NEU))*100,
  "Absolute_Abundance"=NEU,
  "Region"= "Europa")

#
ME<-MFHP %>%
  filter(Region == "MiddleEast/NorthAfrica")

unique(ME$Haplogroup1)
NME <- c(count(ME$Haplogroup1=="H2a"),count(ME$Haplogroup1=="M30"),
         count(ME$Haplogroup1=="T2b"),count(ME$Haplogroup1=="L0a"),
         count(ME$Haplogroup1=="H13"),count(ME$Haplogroup1=="HV"),
         count(ME$Haplogroup1=="U1a"),
         count(ME$Haplogroup1=="H"),count(ME$Haplogroup1=="R5a"),
         count(ME$Haplogroup1=="T1a"),count(ME$Haplogroup1=="R0a"),
         count(ME$Haplogroup1=="U7"),count(ME$Haplogroup1=="H1"),
         count(ME$Haplogroup1=="L1b"),
         count(ME$Haplogroup1=="L2a"),count(ME$Haplogroup1=="U5a"),
         count(ME$Haplogroup1=="H7b"),count(ME$Haplogroup1=="H5a"))

sum(NME)
ME<-data.frame(
  "Haplogroup"=c("H2a","M30","T2b","L0a","H13","HV","U1a","H","R5a","T1a","R0a",
                 "U7","H1","L1b","L2a","U5a","H7b","H5a"),
  "ab_rel"=(NME/sum(NME))*100,
  "Absolute_Abundance"=NME,
  "Region"= "Medio Oriente")

#
EA<-MFHP %>%
  filter(Region == "EastAsia")

m<-MFHP %>%
         filter(Region == "mixed")

EA<-rbind(EA,m)

unique(EA$Haplogroup1)

NEA <- c(count(EA$Haplogroup1=="U4a"),count(EA$Haplogroup1=="D5a"),
         count(EA$Haplogroup1=="C4"),count(EA$Haplogroup1=="D4j"),
         count(EA$Haplogroup1=="Z3"),count(EA$Haplogroup1=="M9a"),
         count(EA$Haplogroup1=="F2"),
         count(EA$Haplogroup1=="B4c"),count(EA$Haplogroup1=="C4a"),
         count(EA$Haplogroup1=="F1"),count(EA$Haplogroup1=="G2a"),
         count(EA$Haplogroup1=="F1a"),count(EA$Haplogroup1=="F1b"),
         count(EA$Haplogroup1=="H2a"),
         count(EA$Haplogroup1=="M18"),count(EA$Haplogroup1=="R2"),
         count(EA$Haplogroup1=="U1a"),count(EA$Haplogroup1=="H"),
         count(EA$Haplogroup1=="U5a"))

sum(NEA)
EA<-data.frame(
  "Haplogroup"=c("U4a","D5a","C4","D4j","Z3","M9a","F2","B4c","C4a","F1","G2a",
                 "F1a","F1b","H2a","M18","R2","U1a","H","U5a"),
  "ab_rel"=(NEA/sum(NEA))*100,
  "Absolute_Abundance"=NEA,
  "Region"= "Asia Oriental")


#
AF<-MFHP %>%
  filter(Region == "Africa")

unique(AF$Haplogroup1)
NAF <- c(count(AF$Haplogroup1=="L1c"),count(AF$Haplogroup1=="L0a"),
             count(AF$Haplogroup1=="L1b"),count(AF$Haplogroup1=="L2a"),
             count(AF$Haplogroup1=="L3d"))

sum(NAF)
AF<-data.frame(
  "Haplogroup"=c("L1c","L0a","L1b","L2a","L3d"),
  "ab_rel"=(NAF/sum(NAF))*100,
  "Absolute_Abundance"=NAF,
  "Region"= "África")

######
#Graficación de la distribución mundial de los haplogrupos presentes en Pakistán

library(paletteer)

unique(MFHP$Region)
MFHP2<-rbind(EU,ME,EA,AF)
unique(MFHP2$Haplogroup)

colnames(MFHP2)<-c("Haplogroup","Relative_Abundance","Absolute_Abundance","Region")

ggplot(MFHP2, aes(fill=Haplogroup, x=Region, y=Relative_Abundance)) + 
  labs(title = "Relative abundance of mitichondrial haplogroups in the world.") +
  geom_bar(stack="dodge", stat="identity") +
  scale_fill_manual(values=paletteer_d("ggsci::default_igv"))+
  theme(axis.text.x = element_text(size = 8, hjust = 1, angle = 40),
        plot.title.position="plot",plot.title = element_text(family = "serif",
                                                             face = "bold",
                                                             size = 15,
                                                             hjust = .8,
                                                             vjust = 1))+
  scale_y_continuous(breaks = seq(0,200,20))

ggplot(MFHP2, aes(fill=Haplogroup, x=Region, y=Absolute_Abundance)) + 
  labs(title = "Abundance of Y chromosome haplogroups in the world.") +
  geom_bar(stack="dodge", stat="identity") 
