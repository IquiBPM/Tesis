setwd("C:/Users/user/Downloads")
library("ggplot2")

######
#Gráfica del porcentaje de aparición de las regiones de origen
#de los haplogrupos del cromosoma Y en Pakistán

popPY<-read.csv("tablaY_regiones.txt", sep="\t")
popPY[c(which(popPY[,5]=="fr")),5]<-"África"

HY<-unique(popPY[,5])

nY <- c()
for (i in 1:length(HY)) {
  nY[i] <- length(which(popPY$Reg==HY[i]))
}
FY<-data.frame("Región"= HY,
                     "n" = nY)
rm(i, nY)

FY2<-FY[order(FY[,1]),]
FY2$per<-(FY2$n/sum(FY2$n))*100
FY2$per<-as.character(FY2$per)
per<-strsplit(x=FY2$per,split = "")
per<-lapply(per, '[', c(1:4))
per<-lapply(per, FUN=paste, collapse="")
FY2$perc<-unlist(per)
FY2$ymax<-cumsum(FY2$per)
FY2$ymin = c(0, head(FY2$ymax, n=-1))
FY2$labelPosition <- ((FY2$ymax + FY2$ymin) / 2)
FY2$label <- paste0(FY2$perc,"%", "\n", FY2$n)
ggplot(FY2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Región)) +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C")) +
  geom_rect(linewidth=7) +
  geom_text( x=4.3, aes(y=labelPosition, label=label), size=3.25) + 
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  ggtitle("Haplogrupos del cromosoma Y con mayor porcentaje de aparición")

#Tabla de porcentajes Y
HY<-strsplit(x=popPY$Haplogrupo,split = "")
HY2<-lapply(HY, '[', 1)
HYF<-lapply(HY2, FUN=paste, collapse="")

popPY$Haplogroup2<-unlist(HYF)

f_summary2 <- dcast(melt(popPY),  popPY[,2] ~ popPY[,6])
ncol(f_summary2)
f_summary2$samples <-  apply(f_summary2[,2:11], 1, sum)
f_summary2[,2:11] <- f_summary2[,2:11]/f_summary2[,"samples"]
f_summary2[,2:11]<-f_summary2[,2:11]*100
write.csv(f_summary2,file="porcentaje de haplogrupos Y por cada grupo", row.names = F)

######
#Gráfica del porcentaje de aparición de las regiones de origen
#de los haplogrupos del DNA mitocondriál en Pakistán

popPY<-read.csv("tablaMt_regiones.txt", sep="\t")
popPY[c(which(popPY[,5]=="fr")),5]<-"África"

HY<-unique(popPY[,5])

nY <- c()
for (i in 1:length(HY)) {
  nY[i] <- length(which(popPY$Reg==HY[i]))
}
FY<-data.frame("Región"= HY,
               "n" = nY)
rm(i, nY)

FY2<-FY[order(FY[,1]),]
FY2$per<-(FY2$n/sum(FY2$n))*100
FY2$per<-as.character(FY2$per)
per<-strsplit(x=FY2$per,split = "")
per<-lapply(per, '[', c(1:4))
per<-lapply(per, FUN=paste, collapse="")
FY2$perc<-unlist(per)
FY2[2,4]<-c("75")
FY2$ymax<-cumsum(FY2$per)
FY2$ymin = c(0, head(FY2$ymax, n=-1))
FY2$labelPosition <- ((FY2$ymax + FY2$ymin) / 2)
FY2$label <- paste0(FY2$perc,"%", "\n", FY2$n)
ggplot(FY2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Región)) +
  scale_fill_manual(values = c("África"="#FB9A99", "Asia Central y del Sur"="#B2DF8A", 
                               "Asia Oriental"="#33A02C", "Europa"="#A6CEE3", "Medio Oriente"="#E31A1C")) +
  geom_rect(linewidth=7) +
  geom_text( x=4.3, aes(y=labelPosition, label=label), size=3.25) + 
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  ggtitle("Haplogrupos del DNA Mt con mayor porcentaje de aparición")

#Tabla de porcentajes Mt
HY<-strsplit(x=popPY$Haplogrupo,split = "")
HY2<-lapply(HY, '[', 1)
HYF<-lapply(HY2, FUN=paste, collapse="")

popPY$Haplogroup2<-unlist(HYF)

f_summary2 <- dcast(melt(popPY),  popPY[,2] ~ popPY[,6])
ncol(f_summary2)
f_summary2$samples <-  apply(f_summary2[,2:19], 1, sum)
f_summary2[,2:19] <- f_summary2[,2:19]/f_summary2[,"samples"]
f_summary2[,2:19]<-f_summary2[,2:19]*100
write.csv(f_summary2,file="porcentaje de haplogrupos Mt por cada grupo", row.names = F)
