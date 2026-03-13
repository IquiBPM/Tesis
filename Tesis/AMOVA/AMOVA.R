

AMOVA_B <- readr::read_tsv("C:/Users/user/Downloads/Secuencias_Pakistan_AMOVA.tsv",
                           col_names = TRUE, col_types = c(rep("c",664419)))

#########Preparación de los datos
######
#Agragar la columna de provinicas
Balochistan<-AMOVA_B[which(AMOVA_B[,3]=="Balochistan"),]
Khyber<-AMOVA_B[which(AMOVA_B[,3]=="Khyber Pakhtunkhwa"),]
Punjab<-AMOVA_B[which(AMOVA_B[,3]=="Punjab"),]
Sindh<-AMOVA_B[which(AMOVA_B[,3]=="Sindh"),]
Gilgit<-AMOVA_B[which(AMOVA_B[,3]=="Gilgit-Baltistan"),]
Gilgit2<-rbind(Gilgit,Khyber)

######
#Agregar la colimna de regiones
unique(AMOVA_B[,5])

Af<-AMOVA_B[which(AMOVA_B[,5]=="Africa"),]
EA<-AMOVA_B[which(AMOVA_B[,5]=="East_Asia"),]
Eu<-AMOVA_B[which(AMOVA_B[,5]=="Europe"),]
CSA<-AMOVA_B[which(AMOVA_B[,5]=="Central_South_Asia"),]
ME<-AMOVA_B[which(AMOVA_B[,5]=="Middle_East"),]

######
#Separación de los grupos (un grupo vs todos los demás)
unique(AMOVA_B[,2])

Brahui<-AMOVA_B[which(AMOVA_B[,2]=="Brahui"),]
Brahui$Vs<-c(rep("Brahui",23))
NoBrahui<-AMOVA_B[which(AMOVA_B[,2]!="Brahui"),]
NoBrahui$Vs<-c(rep("No_Brahui",170))
Brahui<-rbind(Brahui,NoBrahui)
rm(NoBrahui)

Balochi<-AMOVA_B[which(AMOVA_B[,2]=="Balochi"),]
Balochi$Vs<-c(rep("Balochi",22))
NoBalochi<-AMOVA_B[which(AMOVA_B[,2]!="Balochi"),]
NoBalochi$Vs<-c(rep("No_Balochi",171))
Balochi<-rbind(Balochi,NoBalochi)
rm(NoBalochi)

Hazara<-AMOVA_B[which(AMOVA_B[,2]=="Hazara"),]
Hazara$Vs<-c(rep("Hazara",13))
NoHazara<-AMOVA_B[which(AMOVA_B[,2]!="Hazara"),]
NoHazara$Vs<-c(rep("No_Hazara",180))
Hazara<-rbind(Hazara,NoHazara)
rm(NoHazara)

Makrani<-AMOVA_B[which(AMOVA_B[,2]=="Makrani"),]
Makrani$Vs<-c(rep("Makrani",19))
NoMakrani<-AMOVA_B[which(AMOVA_B[,2]!="Makrani"),]
NoMakrani$Vs<-c(rep("No_Makrani",174))
Makrani<-rbind(Makrani,NoMakrani)
rm(NoMakrani)

Sindhi<-AMOVA_B[which(AMOVA_B[,2]=="Sindhi"),]
Sindhi$Vs<-c(rep("Sindhi",18))
NoSindhi<-AMOVA_B[which(AMOVA_B[,2]!="Sindhi"),]
NoSindhi$Vs<-c(rep("No_Sindhi",175))
Sindhi<-rbind(Sindhi,NoSindhi)
rm(NoSindhi)

Pathan<-AMOVA_B[which(AMOVA_B[,2]=="Pathan"),]
Pathan$Vs<-c(rep("Pathan",18))
NoPathan<-AMOVA_B[which(AMOVA_B[,2]!="Pathan"),]
NoPathan$Vs<-c(rep("No_Pathan",175))
Pathan<-rbind(Pathan,NoPathan)
rm(NoPathan)

Kalash<-AMOVA_B[which(AMOVA_B[,2]=="Kalash"),]
Kalash$Vs<-c(rep("Kalash",16))
NoKalash<-AMOVA_B[which(AMOVA_B[,2]!="Kalash"),]
NoKalash$Vs<-c(rep("No_Kalash",177))
Kalash<-rbind(Kalash,NoKalash)
rm(NoKalash,Kalash3)

Burusho<-AMOVA_B[which(AMOVA_B[,2]=="Burusho"),]
Burusho$Vs<-c(rep("Burusho",19))
NoBurusho<-AMOVA_B[which(AMOVA_B[,2]!="Burusho"),]
NoBurusho$Vs<-c(rep("No_Burusho",174))
Burusho<-rbind(Burusho,NoBurusho)
rm(NoBurusho)

Punjabi<-AMOVA_B[which(AMOVA_B[,2]=="Punjabi"),]
Punjabi$Vs<-c(rep("Punjabi",45))
NoPunjabi<-AMOVA_B[which(AMOVA_B[,2]!="Punjabi"),]
NoPunjabi$Vs<-c(rep("No_Punjabi",148))
Punjabi<-rbind(Punjabi,NoPunjabi)
rm(NoPunjabi)



library("poppr")
set.seed(1999)

#########Análisis AMOVA
######
#Balochistan
Balochistan2<-Balochistan[, -c(1,2,3,4,5)]
colnames(Balochistan2)<-c(1:664416)
Balochistan3<-cbind(subset(Balochistan, select = c(Sample_id,Pop,Province,Haplogroup,Region)),Balochistan2)

amova_Y <- df2genind(Balochistan3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(Balochistan3$Pop))
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

######
#Khyber
Khyber2<-Khyber[, -c(1,2,3,4,5)]
colnames(Khyber2)<-c(1:664416)
Khyber3<-cbind(subset(Khyber, select = c(Sample_id,Pop,Province,Haplogroup,Region)),Khyber2)
rm(Khybe,Khyber2,Khybe3)
amova_Y <- df2genind(Khyber3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(Khyber3$Pop))
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

######
#Punjab
Punjab2<-Punjab[, -c(1,2,3,4,5)]
colnames(Punjab2)<-c(1:664416)
Punjab3<-cbind(subset(Punjab, select = c(Sample_id,Pop,Province,Haplogroup,Region)),Punjab2)

amova_Y <- df2genind(Punjab3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(Punjab3$Pop))
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

######
#Sindh
Sindh2<-Sindh[, -c(1,2,3,4,5)]
colnames(Sindh2)<-c(1:664416)
Sindh3<-cbind(subset(Sindh, select = c(Sample_id,Pop,Province,Haplogroup,Region)),Sindh2)

amova_Y <- df2genind(Sindh3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(Sindh3$Pop))
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

######
#Gilgit
Gilgit3<-Gilgit2[, -c(1,2,3,4,5)]
colnames(Gilgit3)<-c(1:664416)
Gilgit4<-cbind(subset(Gilgit2, select = c(Sample_id,Pop,Province,Haplogroup,Region)),Gilgit3)
rm(Gilgit2,Gilgit3)
amova_Y <- df2genind(Gilgit4[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(Gilgit4$Pop))
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)


######
#Middle East
ME2<-ME[, -c(1,2,3,4,5)]
colnames(ME2)<-c(1:664416)
ME3<-cbind(subset(ME, select = c(Sample_id,Pop,Province,Haplogroup,Region)),ME2)
rm(ME,ME2)
amova_Y <- df2genind(ME3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(ME3$Pop), strata = ME3[,c("Province","Pop")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Province/Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

######
#East Asia
EA2<-EA[, -c(1,2,3,4,5)]
colnames(EA2)<-c(1:664416)
EA3<-cbind(subset(EA, select = c(Sample_id,Pop,Province,Haplogroup,Region)),EA2)
rm(EA,EA2)
amova_Y <- df2genind(EA3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(EA3$Pop), strata = EA3[,c("Province","Pop")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Province/Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

######
#Europe
Eu2<-Eu[, -c(1,2,3,4,5)]
colnames(Eu2)<-c(1:664416)
Eu3<-cbind(subset(Eu, select = c(Sample_id,Pop,Province,Haplogroup,Region)),Eu2)
rm(Eu,Eu2)
amova_Y <- df2genind(Eu3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(Eu3$Pop), strata = Eu3[,c("Province","Pop")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Province/Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

amova_Y2 <- df2genind(Eu3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(Eu3$Pop))
test_amova_y <- poppr.amova(amova_Y, ~Pop)

######
#Central_S_Asia
CSA2<-CSA[, -c(1,2,3,4,5)]
colnames(CSA2)<-c(1:664416)
CSA3<-cbind(subset(CSA, select = c(Sample_id,Pop,Province,Haplogroup,Region)),CSA2)
rm(CSA,CSA2)
amova_Y <- df2genind(CSA3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(CSA3$Pop), strata = CSA3[,c("Province","Pop")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Province/Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

######
#Todas las ascendencias
AMOVA_B2<-AMOVA_B[, -c(1,2,3,4,5)]
colnames(AMOVA_B2)<-c(1:664416)
AMOVA_B3<-cbind(subset(AMOVA_B, select = c(Sample_id,Pop,Province,Haplogroup,Region)),AMOVA_B2)

amova_Y <- df2genind(AMOVA_B3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(AMOVA_B3$Pop), strata = AMOVA_B3[,c("Province","Pop")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Province/Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

######
#Comparación de todos los grupos
AMOVA_B2<-AMOVA_B[, -c(1,2,3,4,5)]
colnames(AMOVA_B2)<-c(1:664416)
AMOVA_B3<-cbind(subset(AMOVA_B, select = c(Sample_id,Pop,Province,Haplogroup,Region)),AMOVA_B2)
rm(AMOVA_B2)
amova_Y <- df2genind(AMOVA_B3[,-c(1,2,3,4,5)], ploidy=1, pop = as.factor(AMOVA_B3$Pop))
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

######
#Un solo grupo contra todos los demás

#Brahui
Brahui2<-Brahui[, -c(1,2,3,4,5,664420)]
colnames(Brahui2)<-c(1:664416)
Brahui3<-cbind(subset(Brahui, select = c(Sample_id,Pop,Province,Haplogroup,Region,Vs)),Brahui2)
rm(Brahui,Brahui2)
amova_Y <- df2genind(Brahui3[,-c(1,2,3,4,5,664420)], ploidy=1, pop = as.factor(Brahui3$Pop), strata = Brahui3[,c("Pop","Vs")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Vs/Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

#Balochi
Balochi2<-Balochi[, -c(1,2,3,4,5,664420)]
colnames(Balochi2)<-c(1:664416)
Balochi3<-cbind(subset(Balochi, select = c(Sample_id,Pop,Province,Haplogroup,Region,Vs)),Balochi2)
rm(Balochi,Balochi2)
amova_Y <- df2genind(Balochi3[,-c(1,2,3,4,5,664420)], ploidy=1, pop = as.factor(Balochi3$Pop), strata = Balochi3[,c("Pop","Vs")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Vs/Pop)
test_amova_y_sig   <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

#Hazara
Hazara2<-Hazara[, -c(1,2,3,4,5,664420)]
colnames(Hazara2)<-c(1:664416)
Hazara3<-cbind(subset(Hazara, select = c(Sample_id,Pop,Province,Haplogroup,Region,Vs)),Hazara2)
rm(Hazara,Hazara2)
amova_Y <- df2genind(Hazara3[,-c(1,2,3,4,5,664420)], ploidy=1, pop = as.factor(Hazara3$Pop), strata = Hazara3[,c("Pop","Vs")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Vs/Pop)
test_amova_y_sig <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

#Makrani
Makrani2<-Makrani[, -c(1,2,3,4,5,664420)]
colnames(Makrani2)<-c(1:664416)
Makrani3<-cbind(subset(Makrani, select = c(Sample_id,Pop,Province,Haplogroup,Region,Vs)),Makrani2)
rm(Makrani,Makrani2)
amova_Y <- df2genind(Makrani3[,-c(1,2,3,4,5,664420)], ploidy=1, pop = as.factor(Makrani3$Pop), strata = Makrani3[,c("Pop","Vs")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Vs/Pop)
test_amova_y_sig <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

#Sindhi
Sindhi2<-Sindhi[, -c(1,2,3,4,5,664420)]
colnames(Sindhi2)<-c(1:664416)
Sindhi3<-cbind(subset(Sindhi, select = c(Sample_id,Pop,Province,Haplogroup,Region,Vs)),Sindhi2)
rm(Sindhi,Sindhi2)
amova_Y <- df2genind(Sindhi3[,-c(1,2,3,4,5,664420)], ploidy=1, pop = as.factor(Sindhi3$Pop), strata = Sindhi3[,c("Pop","Vs")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Vs/Pop)
test_amova_y_sig <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

#Pathan
Pathan2<-Pathan[, -c(1,2,3,4,5,664420)]
colnames(Pathan2)<-c(1:664416)
Pathan3<-cbind(subset(Pathan, select = c(Sample_id,Pop,Province,Haplogroup,Region,Vs)),Pathan2)
rm(Pathan,Pathan2)
amova_Y <- df2genind(Pathan3[,-c(1,2,3,4,5,664420)], ploidy=1, pop = as.factor(Pathan3$Pop), strata = Pathan3[,c("Pop","Vs")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Vs/Pop)
test_amova_y_sig <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)
rm(amova_Y,Pathan3,test_amova_y,test_amova_y_sig)

##Kalash
Kalash2<-Kalash[, -c(1,2,3,4,5,664420)]
colnames(Kalash2)<-c(1:664416)
Kalash3<-cbind(subset(Kalash, select = c(Sample_id,Pop,Province,Haplogroup,Region,Vs)),Kalash2)
rm(Kalash,Kalash2)
amova_Y <- df2genind(Kalash3[,-c(1,2,3,4,5,664420)], ploidy=1, pop = as.factor(Kalash3$Pop), strata = Kalash3[,c("Pop","Vs")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Vs/Pop)
test_amova_y_sig <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)
rm(amova_Y,Pathan3,test_amova_y,test_amova_y_sig)

#Burusho
Burusho2<-Burusho[, -c(1,2,3,4,5,664420)]
colnames(Burusho2)<-c(1:664416)
Burusho3<-cbind(subset(Burusho, select = c(Sample_id,Pop,Province,Haplogroup,Region,Vs)),Burusho2)
rm(Burusho,Burusho2)
amova_Y <- df2genind(Burusho3[,-c(1,2,3,4,5,664420)], ploidy=1, pop = as.factor(Burusho3$Pop), strata = Burusho3[,c("Pop","Vs")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Vs/Pop)
test_amova_y_sig <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)

rm(amova_Y,Pathan3,test_amova_y,test_amova_y_sig)

amova_Y2 <- df2genind(Burusho3[,-c(1,2,3,4,5,664420)], ploidy=1, pop = as.factor(Burusho3$Vs))
amova_Y2 <- as.genclone(amova_Y2)
amova_Y2
test_amova_y2 <- poppr.amova(amova_Y2, ~Vs)
test_amova_y_sig2 <- randtest(test_amova_y2, nrepet = 1000)
plot(test_amova_y_sig2)
rm(amova_Y,Pathan3,test_amova_y,test_amova_y_sig)

##Punjabi
Punjabi2<-Punjabi[, -c(1,2,3,4,5,664420)]
colnames(Punjabi2)<-c(1:664416)
Punjabi3<-cbind(subset(Punjabi, select = c(Sample_id,Pop,Province,Haplogroup,Region,Vs)),Punjabi2)
rm(Punjabi,Punjabi2)
amova_Y <- df2genind(Punjabi3[,-c(1,2,3,4,5,664420)], ploidy=1, pop = as.factor(Punjabi3$Pop), strata = Punjabi3[,c("Pop","Vs")] )
amova_Y <- as.genclone(amova_Y)
amova_Y
test_amova_y <- poppr.amova(amova_Y, ~Vs/Pop)
test_amova_y_sig <- randtest(test_amova_y, nrepet = 1000)
plot(test_amova_y_sig)
rm(amova_Y,Punjabi3,test_amova_y,test_amova_y_sig)
