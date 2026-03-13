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
library(knitr)

setwd("C:/Users/IquiB/OneDrive/Escritorio/Iqui/CCG/Análisis de haplogrupos/Mapa Pakistán/Enviar")
haplom<-read.csv("hapMt_Pak.txt")
RegionMt<-read.csv("HaploMT_region_origin3(la_buena).csv")
setwd("C:/Users/IquiB/OneDrive/Escritorio/Iqui/CCG/Análisis de haplogrupos/Mapa Pakistán/Las buenas")
haploy<-read.csv("hapY_Pak.txt")
RegionY<-read.csv("HaploY_region_origin3(la_buenax3).csv")

N <-data.frame(
  "Población"=c("Sindhi","Balochi","Burusho","Pathan","Brahui","Makrani","Kalash","Hazara","Punjabi"),
  
  "Mt"=c(nrow(haplom[which(haplom$Population=="Sindhi"),]),nrow(haplom[which(haplom$Population=="Balochi"),]),
        nrow(haplom[which(haplom$Population=="Burusho"),]),nrow(haplom[which(haplom$Population=="Pathan"),]),
        nrow(haplom[which(haplom$Population=="Brahui"),]),nrow(haplom[which(haplom$Population=="Makrani"),]),
        nrow(haplom[which(haplom$Population=="Kalash"),]),nrow(haplom[which(haplom$Population=="Hazara"),]),
        nrow(haplom[which(haplom$Population=="Punjabi"),])),
  
  "Y"=c(nrow(haploy[which(haploy$Population=="Sindhi"),]),nrow(haploy[which(haploy$Population=="Balochi"),]),
        nrow(haploy[which(haploy$Population=="Burusho"),]),nrow(haploy[which(haploy$Population=="Pathan"),]),
        nrow(haploy[which(haploy$Population=="Brahui"),]),nrow(haploy[which(haploy$Population=="Makrani"),]),
        nrow(haploy[which(haploy$Population=="Kalash"),]),nrow(haploy[which(haploy$Population=="Hazara"),]),
        nrow(haploy[which(haploy$Population=="Punjabi"),])))

write.csv(N, file="número de muestras.csv",quote=FALSE)

###########################

help(wilcox.test)
rank(RegionMt[])

prub<-data.frame("Middle_east" = c(RegionY[2,5],RegionMt[2,5]),
                 "Sex"= c("Y","Mt"))

wilcox.test(Middle_east ~ Sex, data=prub)

D<-(5*9)

prub2<-data.frame("Middle_east" = c(RegionY[5,2],RegionY[5,3],
                                   RegionY[5,4],RegionY[5,5],
                                   RegionY[5,6],RegionY[5,7],
                                   RegionY[5,8],RegionY[5,9],
                                   RegionY[5,10]
                                   ,RegionMt[5,2],RegionMt[5,3]
                                   ,RegionMt[5,4],RegionMt[5,5],
                                   RegionMt[5,6],RegionMt[5,7]
                                   ,RegionMt[5,8],RegionMt[5,9],
                                   RegionMt[5,10]),
                 "Sex"= c(rep("Y",9),rep("Mt",9)))

wilcox.test(Middle_east ~ Sex, data=prub2)


nMt+(-1*nY)/nMt+nY

##"nMt= número de muesrtas mitocondriales de un 
#grupo correspondoientes a una región.
#nY = némro de muestras del cromosoma Y de un 
#grupo correspondientes a una región.
#Tn = muestras totales de un grupo correspondoientes a una región.

Y<-c(RegionY[1:5,2],RegionY[1:5,3])

RegionY[1:5,2]+RegionMt[1:5,2]

#########################
FBalMt<-rep(0, 24)
FBalY<-c(1,1,rep(0,22))
WFBal<-data.frame(FBalY, FBalMt)
Bal_Af<-wilcox.test(FBalY, FBalMt, exact=F, data=WFBal, alternative = "two.side",
            paired=F)
rm(FBalMt,FBalY,WFBal)

FMaMt<-c(rep(1, 7),rep(0,13))
FMaY<-c(1,rep(0,19))
WFMa<-data.frame(FMaY, FMaMt)
Mak_Af<-wilcox.test(FMaY, FMaMt, exact=F, data=WFMa, alternative = c("two.side"))
rm(FMaMt,FMaY,WFMa)

FSiMt<-c(rep(1, 3),rep(0,17))
FSiY<-rep(0,20)
WFSi<-data.frame(FSiY, FSiMt)
Sin_Af<-wilcox.test(FSiY, FSiMt, exact=F, data=WFSi, alternative = c("two.side"))
rm(FSiMt,WFSi,FSiY)


CSABalMt<-c(rep(1, 20),rep(0,4))
CSABalY<-c(rep(0,10),rep(1,14))
WCSABal<-data.frame(CSABalY, CSABalMt)
Bal_CSA<-wilcox.test(CSABalY, CSABalMt, exact=F, data=WCSABal)
rm(CSABalMt,CSABalY,WCSABal)

CSABrMt<-c(rep(1, 16),rep(0,9))
CSABrY<-c(rep(0,9),rep(1,16))
WCSABr<-data.frame(CSABrY, CSABrMt)
Bra_CSA<-wilcox.test(CSABrY, CSABrMt, exact=F, data=WCSABr, alternative = "greater")
rm(CSABrMt,CSABrY,WCSABr)

CSABuMt<-c(rep(1, 14),rep(0,6))
CSABuY<-c(rep(0,5),rep(1,15))
WCSABu<-data.frame(CSABuY, CSABuMt)
Bur_CSA<-wilcox.test(CSABuY, CSABuMt, exact=F, data=WCSABu, alternative = c("two.side"))
rm(CSABuMt,CSABuY,WCSABu)

CSAHaMt<-c(rep(1, 12),rep(0,10))
CSAHaY<-c(rep(0,16),rep(1,3),rep(NA,3))
WCSAHa<-data.frame(CSAHaY, CSAHaMt)
Haz_CSA<-wilcox.test(CSAHaY, CSAHaMt, exact=F, data=WCSAHa, alternative = c("two.side"))
rm(CSAHaMt,CSAHaY,WCSAHa)

CSAKaMt<-c(rep(1, 18),rep(0,0))
CSAKaY<-c(rep(0,7),rep(1,10),rep(NA,1))
WCSAKa<-data.frame(CSAKaY, CSAKaMt)
Kal_CSA<-wilcox.test(CSAKaY, CSAKaMt, exact=F, data=WCSAKa, alternative = c("two.side"))
rm(CSAKaMt,CSAKaY,WCSAKa)

CSAMaMt<-c(rep(1, 9),rep(0,11))
CSAMaY<-c(rep(0,4),rep(1,16),rep(NA,0))
WCSAMa<-data.frame(CSAMaY, CSAMaMt)
Mak_CSA<-wilcox.test(CSAMaY, CSAMaMt, exact=F, data=WCSAMa, alternative = c("two.side"))
rm(CSAMaMt,CSAMaY,WCSAMa)

CSAPaMt<-c(rep(1, 15),rep(0,4))
CSAPaY<-c(rep(0,5),rep(1,14),rep(NA,0))
WCSAKa<-data.frame(CSAPaY, CSAPaMt)
Pat_CSA<-wilcox.test(CSAPaY, CSAPaMt, exact=F, data=WCSAPa, alternative = c("two.side"))
rm(CSAPaMt,CSAPaY,WCSAKa)

CSAPuMt<-c(rep(1, 43),rep(0,5))
CSAPuY<-c(rep(0,6),rep(1,42),rep(NA,0))
WCSAPu<-data.frame(CSAPuY, CSAPuMt)
Pun_CSA<-wilcox.test(CSAPuY, CSAPuMt, exact=F, data=WCSAPu, alternative = c("two.side"))
rm(CSAPuMt,CSAPuY,WCSAPu)

CSASiMt<-c(rep(0,15),rep(1, 5))
CSASiY<-c(rep(NA,0),rep(0,2),rep(1,18))
WCSASi<-data.frame(CSASiY, CSASiMt)
Sin_CSA<-wilcox.test(CSASiY, CSASiMt, exact=F, data=WCSASi, alternative = c("two.side"))
rm(CSASiMt,CSASiY,WCSASi)



EABrMt<-c(rep(1,4),rep(0,21))
EABrY<-c(rep(0,24),rep(1,1))
WEABr<-data.frame(EABrY, EABrMt)
Bra_EA<-wilcox.test(EABrY, EABrMt, exact=F, data=WEABr, alternative = c("two.side"))
rm(EABrMt,EABrY,WEABr)

EABuMt<-c(rep(1, 3),rep(0,17))
EABuY<-c(rep(0,18),rep(1,2))
WEABu<-data.frame(EABuY, EABuMt)
Bur_EA<-wilcox.test(EABuY, EABuMt, exact=F, data=WEABu, alternative = c("two.side"))
rm(EABuMt,EABuY,WEABu)

EAHaMt<-c(rep(1, 5),rep(0,17))
EAHaY<-c(rep(0,12),rep(1,7),rep(NA,3))
WEAHa<-data.frame(EAHaY, EAHaMt)
Haz_EA<-wilcox.test(EAHaY, EAHaMt, exact=F, data=WEAHa, alternative = c("two.side"))
rm(EAHaMt,EAHaY,WEAHa)

EAPaMt<-c(rep(1, 1),rep(0,19))
EAPaY<-c(rep(0,20),rep(1,0),rep(NA,0))
WEAKa<-data.frame(EAPaY, EAPaMt)
Pat_EA<-wilcox.test(EAPaY, EAPaMt, exact=F, data=WEAPa, alternative = c("two.side"))
rm(EAPaMt,EAPaY,WEAKa)

EAPuMt<-c(rep(1, 2),rep(0,46))
EAPuY<-c(rep(0,48),rep(1,0),rep(NA,0))
WEAPu<-data.frame(EAPuY, EAPuMt)
Pun_EA<-wilcox.test(EAPuY, EAPuMt, exact=F, data=WEAPu, alternative = c("two.side"))
rm(EAPuMt,EAPuY,WEAPu)


EBalMt<-rep(0, 24)
EBalY<-c(1,1,rep(0,22))
WEBal<-data.frame(EBalY, EBalMt)
Bal_E<-wilcox.test(EBalY, EBalMt, exact=F, data=WEBal, alternative = c("two.side"))
rm(EBalMt,EBalY,WEBal)

EBrMt<-c(1,1,1,1,rep(0, 21))
EBrY<-rep(0,25)
WEBr<-data.frame(EBrY, EBrMt)
Bra_E<-wilcox.test(EBrY, EBrMt, exact=F, data=WEBr, alternative = c("two.side"))
rm(EBrMt,EBrY,WEBr)

EBuMt<-c(1,1,rep(0, 18))
EBuY<-rep(0,20)
WEBu<-data.frame(EBuY, EBuMt)
Bur_E<-wilcox.test(EBuY, EBuMt, exact=F, data=WEBu, alternative = c("two.side"))
rm(EBuMt,EBuY,WEBu)

EHaMt<-c(1,1,1,rep(0, 19))
EHaY<-c(rep(1,9),rep(0,10),rep(NA,3))
WEHa<-data.frame(EHaY, EHaMt)
Haz_E<-wilcox.test(EHaY, EHaMt, exact=F, data=WEHa, alternative = c("two.side"))
rm(EHaMt,EHaY,WEHa)

EMaMt<-c(1,1,1,rep(0, 17))
EMaY<-c(1,rep(0,19))
WEMa<-data.frame(EMaY, EMaMt)
Mak_E<-wilcox.test(EMaY, EMaMt, exact=F, data=WEMa, alternative = c("two.side"))
rm(EMaMt,EMaY,WEMa)

ESiMt<-c(1,rep(0, 19))
ESiY<-rep(0,20)
WESi<-data.frame(ESiY, ESiMt)
Sin_E<-wilcox.test(ESiY, ESiMt, exact=F, data=WESi, alternative = c("two.side"))
rm(ESiMt,ESiY,WESi)


MEBalMt<-c(1,1,1,1,rep(0, 20))
MEBalY<-c(rep(1,6),rep(0,18))
WMEBal<-data.frame(MEBalY, MEBalMt)
Bal_ME<-wilcox.test(MEBalY, MEBalMt, exact=F, data=WMEBal, alternative = c("two.side"))
rm(MEBalMt,MEBalY,WMEBal)

MEBrMt<-c(1,1,1,1,rep(0, 21))
MEBrY<-c(rep(1,8),rep(0,17))
WMEBr<-data.frame(MEBrY, MEBrMt)
Bra_ME<-wilcox.test(MEBrY, MEBrMt, exact=F, data=WMEBr, alternative = c("two.side"))
rm(MEBrMt,MEBrY,WMEBr)

MEBuMt<-c(1,rep(0, 19))
MEBuY<-c(1,1,1,rep(0,17))
WMEBu<-data.frame(MEBuY, MEBuMt)
Bur_ME<-wilcox.test(MEBuY, MEBuMt, exact=F, data=WMEBu, alternative = c("two.side"))
rm(MEBuMt,MEBuY,WMEBu)

MEHaMt<-c(1,1,rep(0, 20))
MEHaY<-c(rep(0,19),rep(NA,3))
WMEHa<-data.frame(MEHaY, MEHaMt)
Haz_ME<-wilcox.test(MEHaY, MEHaMt, exact=F, data=WMEHa, alternative = c("two.side"))
rm(MEHaMt,MEHaY,WMEHa)

MEKaMt<-rep(0, 18)
MEKaY<-c(rep(1,7),rep(0,10),NA)
WMEKa<-data.frame(MEKaY, MEKaMt)
Kal_ME<-wilcox.test(MEKaY, MEKaMt, exact=F, data=WMEKa, alternative = c("two.side"))
rm(MEKaMt,MEKaY,WMEKa)

MEMaMt<-c(1,rep(0, 19))
MEMaY<-c(1,1,rep(0,18))
WMEMa<-data.frame(MEMaY, MEMaMt)
Mak_ME<-wilcox.test(MEMaY, MEMaMt, exact=F, data=WMEMa, alternative = c("two.side"))
rm(MEMaMt,MEMaY,WMEMa)

MEPaMt<-c(1,1,1,rep(0, 16))
MEPaY<-c(rep(1,3),rep(0,16))
WMEPa<-data.frame(MEPaY, MEPaMt)
Pat_ME<-wilcox.test(MEPaY, MEPaMt, exact=F, data=WMEPa, alternative = c("two.side"))
rm(MEPaMt,MEPaY,WMEPa)

MEPuMt<-c(1,1,1,rep(0, 45))
MEPuY<-c(rep(1,6),rep(0,42))
WMEPu<-data.frame(MEPuY, MEPuMt)
Pun_ME<-wilcox.test(MEPuY, MEPuMt, exact=F, data=WMEPu, alternative = c("two.side"))
rm(MEPuMt,MEPuY,WMEPu)

MESiMt<-c(1,rep(0, 19))
MESiY<-c(1,1,rep(0,18))
WMESi<-data.frame(MESiY, MESiMt)
Sin_ME<-wilcox.test(MESiY, MESiMt, exact=F, data=WMESi, alternative = c("two.side"))
rm(MESiMt,MESiY,WMESi)

############################
library(dplyr)
library(knitr)
library(ggplot2)

p_val<-(c(Bal_Af$p.value, Bal_CSA$p.value,Bal_E$p.value, Bal_ME$p.value, Bra_CSA$p.value, Bra_E$p.value,
        Bra_EA$p.value, Bra_ME$p.value, Bur_CSA$p.value, Bur_E$p.value, Bur_EA$p.value, Bur_ME$p.value,
        Haz_CSA$p.value, Haz_E$p.value, Haz_EA$p.value, Haz_ME$p.value, Kal_CSA$p.value, Kal_ME$p.value,
        Mak_Af$p.value, Mak_CSA$p.value, Mak_E$p.value, Mak_ME$p.value, Pat_CSA$p.value, Pat_EA$p.value,
        Pat_ME$p.value, Sin_Af$p.value, Sin_CSA$p.value, Sin_E$p.value, Sin_ME$p.value, Pun_CSA$p.value,
        Pun_EA$p.value, Pun_ME$p.value))

test<-(c("Bal_Af", "Bal_CSA", "Bal_E", "Bal_ME", "Bra_CSA", "Bra_E", "Bra_EA", "Bra_ME", "Bur_CSA", 
          "Bur_E", "Bur_EA", "Bur_ME", "Haz_CSA", "Haz_E", "Haz_EA", "Haz_ME", "Kal_CSA", "Kal_ME", 
          "Mak_Af", "Mak_CSA", "Mak_E", "Mak_ME", "Pat_CSA", "Pat_EA", "Pat_ME", "Sin_Af", "Sin_CSA", 
          "Sin_E", "Sin_ME", "Pun_CSA", "Pun_EA", "Pun_ME"))

resultados <- data.frame(test, p_val)
resultados <- arrange(resultados, p_val)
resultados <- mutate(resultados, indice = 1:length(p_val))
resultados <- mutate(resultados, `d*i/n` = 0.05*(indice/length(p_val)))
resultados <- mutate(resultados, significancia = p_val <= `d*i/n`)
kable(resultados, align = "c")

BH <-p.adjust(p = p_val, method = "BH") <= 0.05
B <- p.adjust(p = p_val, method = "bonferroni") <= 0.05
write.csv(data.frame("Test" = test, "p_value"= p_val,"Bonferroni" = B, "Benjamini_Hochberg" = BH ),
          file="sesgo sexual")

pi_0 <- sum(resultados$p_val > 0.5)/(length(resultados$p_val) * 0.5)
ggplot(data = resultados, aes(resultados$p_val)) + geom_histogram(aes(y = ..density..), 
               fill = "grey", color = "black") + labs(title = "Distribución de los m p-values", 
               x = "p value") + geom_hline(yintercept = 0.625, linetype = "dashed", color = "red", 
               size = 1.5) + geom_segment(aes(x = 0.5, y = 0, xend = 0.5, yend = 0.67), colour = "red", 
               size = 1.5) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))

##############################

sex_bias<-data.frame("Región" = c("África", "Asia Central y del Sur",
                                  "Asia Oriental","Europa",
                                  "Medio Oriente","XYZ"),
                     "Balochi" = c(-1,0.25,NA,-1,-0.2,1),
                     "Brahui" = c(NA,0,0.6,1,-0.5,-1),
                     "Burusho" = c(NA,-0.04,0.25,1,-0.25,-1),
                     "Hazara" = c(NA,0.6,0,-0.5,1,-1),
                     "Kalash" = c(NA,0.28,NA,NA,-1,1),
                     "Makrani" = c(0.29,-0.23,NA,0.5,-0.5,1),
                     "Pathan" = c(NA,0.03,1,0,0,1),
                     "Punjabi" = c(NA,0.01,1,NA,-0.7,-1),
                     "Sindhi" = c(1,-0.56,NA,1,-1,1))

ggplot(sex_bias, aes(fill=sex_bias[,1], y=sex_bias[,2], x=sex_bias[,1])) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Sesgo sexual en el grupo cultural Balochi", y="Sesgo",x="Región") +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C"))+
  theme(axis.text.x = element_text(size = 9, hjust = 1, angle = 40),
        legend.position=c(50,50),
        plot.title = element_text(family = "serif",face = "bold",
                                  size = 15
        ))+
  scale_y_continuous(breaks = seq(-2,1,.2))


ggplot(sex_bias, aes(fill=sex_bias[,1], y=sex_bias[,3], x=sex_bias[,1],)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Sesgo sexual en el grupo cultural Brahui", y="Sesgo",x="Región") +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C"))+
  theme(axis.text.x = element_text(size = 9, hjust = 1, angle = 40),
        legend.position=c(50,50),
        plot.title = element_text(family = "serif",face = "bold",
                                  size = 15
        ))+
  scale_y_continuous(breaks = seq(-2,1,.2))

ggplot(sex_bias, aes(fill=sex_bias[,1], y=sex_bias[,4], x=sex_bias[,1],)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Sesgo sexual en el grupo cultural Burusho", y="Sesgo",x="Región") +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C"))+
  theme(axis.text.x = element_text(size = 9, hjust = 1, angle = 40),
        legend.position=c(50,50),
        plot.title = element_text(family = "serif",face = "bold",
                                  size = 15
        ))+
  scale_y_continuous(breaks = seq(-2,1,.2))

ggplot(sex_bias, aes(fill=sex_bias[,1], y=sex_bias[,5], x=sex_bias[,1],)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Sesgo sexual en el grupo cultural Hazara", y="Sesgo",x="Región") +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C"))+
  theme(axis.text.x = element_text(size = 9, hjust = 1, angle = 40),
        legend.position=c(50,50),
        plot.title = element_text(family = "serif",face = "bold",
                                  size = 15
        ))+
  scale_y_continuous(breaks = seq(-2,1,.2))

ggplot(sex_bias, aes(fill=sex_bias[,1], y=sex_bias[,6], x=sex_bias[,1],)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Sesgo sexual en el grupo cultural Kalash", y="Sesgo",x="Región") +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C"))+
  theme(axis.text.x = element_text(size = 9, hjust = 1, angle = 40),
        legend.position=c(50,50),
        plot.title = element_text(family = "serif",face = "bold",
                                  size = 15
        ))+
  scale_y_continuous(breaks = seq(-2,1,.2))

ggplot(sex_bias, aes(fill=sex_bias[,1], y=sex_bias[,7], x=sex_bias[,1],)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Sesgo sexual en el grupo cultural Makrani", y="Sesgo",x="Región") +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C"))+
  theme(axis.text.x = element_text(size = 9, hjust = 1, angle = 40),
        legend.position=c(50,50),
        plot.title = element_text(family = "serif",face = "bold",
                                  size = 15
        ))+
  scale_y_continuous(breaks = seq(-2,1,.2))

ggplot(sex_bias, aes(fill=sex_bias[,1], y=sex_bias[,8], x=sex_bias[,1],)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Sesgo secual en el grupo cultural Pathan", y="Sesgo",x="Región") +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C"))+
  theme(axis.text.x = element_text(size = 9, hjust = 1, angle = 40),
        legend.position=c(50,50),
        plot.title = element_text(family = "serif",face = "bold",
                                  size = 15
        ))+
  scale_y_continuous(breaks = seq(-2,1,.2))

ggplot(sex_bias, aes(fill=sex_bias[,1], y=sex_bias[,9], x=sex_bias[,1],)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Sesgo secual en el grupo cultural Punjabi", y="Sesgo",x="Región") +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C"))+
  theme(axis.text.x = element_text(size = 9, hjust = 1, angle = 40),
        legend.position=c(50,50),
        plot.title = element_text(family = "serif",face = "bold",
                                  size = 15
        ))+
  scale_y_continuous(breaks = seq(-2,1,.2))

ggplot(sex_bias, aes(fill=sex_bias[,1], y=sex_bias[,10], x=sex_bias[,1],)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Sesgo sexual en el grupo sexual Sindhi", y="Sesgo",x="Región") +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C"))+
  theme(axis.text.x = element_text(size = 9, hjust = 1, angle = 40),
        legend.position=c(50,50),
        plot.title = element_text(family = "serif",face = "bold",
                                  size = 15
        ))+
  scale_y_continuous(breaks = seq(-2,1,.2))
