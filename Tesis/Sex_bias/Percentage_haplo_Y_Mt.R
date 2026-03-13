setwd("C:/Users/user/Downloads")
library("ggplot2")
popPY<-read.csv("hapY_Pak.txt")

popPMt<-read.csv("hapMt_Pak.txt")

##hapMt2<-popPMt[,c("Sample.ID","Population","MtDNA.haplogroup")]

##write.csv(hapMt2,"hapMt2_Pak.txt",col.names = T)

length(which(popPY$Haplogroup=="R1a1a1b2a2a"))

##lapply(HY, FUN=length(which(popPY$Haplogroup=="HY")))

HY<-unique(popPY[,3])
nY <- c()
for (i in 1:length(HY)) {
  nY[i] <- length(which(popPY$Haplogroup==HY[i]))
}
n_por_hY<-data.frame("Haplogrupo"= HY,
           "n" = nY)

otros<-data.frame("Haplogrupo" = "Otros",
           "n" = sum(c(n_por_hY[order(n_por_hY$n,decreasing = T)[10:length(n_por_hY$n)],2])))
FY<-rbind(n_por_hY[c(order(n_por_hY$n,decreasing = TRUE)[1:9]),],
          otros)
rm(otros, i, nY)

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
ggplot(FY2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Haplogrupo)) +
  scale_fill_brewer(palette="Paired") +
  geom_rect(linewidth=7) +
  geom_text( x=4.3, aes(y=labelPosition, label=label), size=3.25) + 
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  ggtitle("Haplogrupos del cromosoma Y con mayor porcentaje de aparición")



M1<-strsplit(x=popPY$Haplogroup,split = "")
M1<-lapply(M1, '[', c(1))
M1<-lapply(M1, FUN=na.omit)
M1<-lapply(M1, FUN=paste, collapse="")
popPY$Haplogroup1<-unlist(M1)
HY<-unique(popPY[,4])

nY <- c()
for (i in 1:length(HY)) {
  nY[i] <- length(which(popPY$Haplogroup1==HY[i]))
}
n_por_hY<-data.frame("Haplogrupo"= HY,
                     "n" = nY)
otros<-data.frame("Haplogrupo" = "Otros",
                  "n" = sum(c(n_por_hY[order(n_por_hY$n,decreasing = T)[8:length(n_por_hY$n)],2])))
FY<-rbind(n_por_hY[c(order(n_por_hY$n,decreasing = TRUE)[1:7]),],
          otros)
rm(otros, i, nY)

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
ggplot(FY2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Haplogrupo)) +
  scale_fill_brewer(palette="Paired") +
  geom_rect(linewidth=7) +
  geom_text( x=4.3, aes(y=labelPosition, label=label), size=3.25) + 
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  ggtitle("Haplogrupos del cromosoma Y simplificados con mayor porcentaje de aparición")

############################################3
HMt<-unique(popPMt[,2])
nMt <- c()
for (i in 1:length(HMt)) {
  nMt[i] <- length(which(popPMt$Haplogroup==HMt[i]))
}
n_por_hMt<-data.frame("Haplogrupo"= HMt,
                     "n" = nMt)
otras<-data.frame("Haplogrupo" = "Otros",
                  "n" = sum(c(n_por_hMt[order(n_por_hMt$n,decreasing = T)[6:length(n_por_hMt$n)],2])))
FMt<-rbind(n_por_hMt[c(order(n_por_hMt$n,decreasing = TRUE)[1:5]),],
          otras)
rm(otras, i, nMt)

FMt2<-FMt[order(FMt[,1]),]
FMt2$per<-(FMt2$n/sum(FMt2$n))*100
FMt2$per<-as.character(FMt2$per)
per<-strsplit(x=FMt2$per,split = "")
per<-lapply(per, '[', c(1:4))
per<-lapply(per, FUN=paste, collapse="")
FMt2$perc<-unlist(per)
FMt2$ymax<-cumsum(FMt2$per)
FMt2$ymin = c(0, head(FMt2$ymax, n=-1))
FMt2$labelPosition <- ((FMt2$ymax + FMt2$ymin) / 2)
FMt2$label <- paste0(FMt2$perc,"%", "\n", FMt2$n)
ggplot(FMt2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Haplogrupo)) +
  geom_rect(linewidth=7) +
  scale_fill_brewer(palette="Paired") +
  geom_text( x=4.3, aes(y=labelPosition, label=label), size=3.25) + 
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  ggtitle("Haplogrupos del DNA mitocondrial con mayor porcentaje de aparición")


Mt1<-strsplit(x=popPMt$Haplogroup,split = "")
Mt1<-lapply(Mt1, '[', c(1))
Mt1[[61]]<-gsub("[+]",NA,Mt1[[61]])
Mt1[[71]]<-gsub("[+]",NA,Mt1[[71]])
Mt1[[126]]<-gsub("[+]",NA,Mt1[[126]])
Mt1[[161]]<-gsub("[+]",NA,Mt1[[161]])
Mt1[[163]]<-gsub("[+]",NA,Mt1[[163]])
Mt1<-lapply(Mt1, FUN=na.omit)
Mt1<-lapply(Mt1, FUN=paste, collapse="")
popPMt$Haplogroup1<-unlist(Mt1)

HMt<-unique(popPMt[,4])
nMt <- c()
for (i in 1:length(HMt)) {
  nMt[i] <- length(which(popPMt$Haplogroup1==HMt[i]))
}
n_por_hMt<-data.frame("Haplogrupo"= HMt,
                      "n" = nMt)
otras<-data.frame("Haplogrupo" = "Otros",
                  "n" = sum(c(n_por_hMt[order(n_por_hMt$n,decreasing = T)[9:length(n_por_hMt$n)],2])))
FMt<-rbind(n_por_hMt[c(order(n_por_hMt$n,decreasing = TRUE)[1:8]),],
           otras)
rm(otras, i, nMt)

FMt2<-FMt[order(FMt[,1]),]
FMt2$per<-(FMt2$n/sum(FMt2$n))*100
FMt2$per<-as.character(FMt2$per)
per<-strsplit(x=FMt2$per,split = "")
per<-lapply(per, '[', c(1:4))
per<-lapply(per, FUN=paste, collapse="")
FMt2$perc<-unlist(per)
FMt2$ymax<-cumsum(FMt2$per)
FMt2$ymin = c(0, head(FMt2$ymax, n=-1))
FMt2$labelPosition <- ((FMt2$ymax + FMt2$ymin) / 2)
FMt2$label <- paste0(FMt2$perc,"%", "\n", FMt2$n)
ggplot(FMt2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Haplogrupo)) +
  geom_rect(linewidth=7) +
  scale_fill_brewer(palette="Paired") +
  geom_text( x=4.3, aes(y=labelPosition, label=label), size=3.25) + 
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  ggtitle("Haplogrupos del DNA mitocondrial simplificados con mayor porcentaje de aparición")

