library(dplyr)
library(knitr)
library(ggplot2)

######
#Corección para las regiones

Provincia<-c("Balochi","Gilgit","Khyber",
              "Punjab","Sindh")
p_val_prov<-c(0.2197,0.0039,0.0039,
         0.1388611,0.4895105)
B_P <- p.adjust(p = p_val_prov, method = "bonferroni") <= 0.05
BH_P <-p.adjust(p = p_val_prov, method = "BH") <= 0.05
kable(data.frame("Provincia" = Provincia, "p_value"= p_val_prov, "Bonferroni" = B_P,
                                     "Benjamini_Hochberg" = BH_P), align = "c")
##write.csv(data.frame("Provincia" = Provincia, "p_value"= p_val_prov, "Bonferroni" = B_P,
  #                   "Benjamini_Hochberg" = BH_P),
   #     file="Porvincias_AMOVA", row.names = F)

######
#Corección para las regiones

Ascendencia<-c("Medio Oriente: entre provincias","Medio Oriente: entre grupos","Medio Oriente: dentro de los grupos",
               "Asia Oriental: dentro de los grupos",
               "Europa: dentro de los grupos",
               "Asia Central y del Sur: entre provincias","Asia Central y del Sur: entre grupos","Asia Central y del Sur: dentro de los grupos")
p_val_ascen<-c(0.1638,0.3596,0.0679,
               0.3876,
               0.1478,
               0.0419, 0.2687,0.0419)
B_A <- p.adjust(p = p_val_ascen, method = "bonferroni") <= 0.05
BH_A <-p.adjust(p = p_val_ascen, method = "BH") <= 0.05
kable(data.frame("Ascendencia" = Ascendencia, "p_value"= p_val_ascen, "Bonferroni" = B_A,
                 "Benjamini_Hochberg" = BH_A), align = "c")

#write.csv(data.frame("Ascendencia" = Ascendencia, "p_value"= p_val_ascen, "Bonferroni" = B_A,
  #                   "Benjamini_Hochberg" = BH_A),
   #  file="Ascendencias_AMOVA", row.names = F)

######
#Corección para cada grupo

Grupo<-c("Brahui: contra los demás","Brahui: entre grupos","Brahui: dentro de los grupos",
         "Balochi: contra los demás","Balochi: entre grupos","Balochi: dentro de los grupos",
         "Hazara: contra los demás","Hazara: entre grupos","Hazara: dentro de los grupos",
         "Makrani: contra los demás","Makrani: entre grupos","Makrani: dentro de los grupos",
         "Sindhi: contra los demás","Sindhi: entre grupos","Sindhi: dentro de los grupos",
         "Pathan: contra los demás","Pathan: entre grupos","Pathan: dentro de los grupos",
         "Kalash: contra los demás","Kalash: entre grupos","Kalash: dentro de los grupos",
         "Burusho: contra los demás","Burusho: entre grupos","Burusho: dentro de los grupos",
         "Punjabi: contra los demás","Punjabi: entre grupos","Punjabi: dentro de los grupos")
p_val_grup<-c(0.4415,0.0009,0.0009,
               0.901,0.0009,0.0029,
               0.1328,0.0129,0.0029,
               1,0.0009,0.0019,
               0.6783,0.0009,0.0019,
               0.5404,0.0009,0.0009,
               0.3266,0.001,0.0009,
               0.2047,0.0049,0.0009,
              0.7942,0.001,0.0009)
B_G <- p.adjust(p = p_val_grup, method = "bonferroni") <= 0.05
BH_G <-p.adjust(p = p_val_grup, method = "BH") <= 0.05
kable(data.frame("Grupo" = Grupo, "p_value"= p_val_grup, "Bonferroni" = B_G,
                 "Benjamini_Hochberg" = BH_G), align = "c")
write.csv(data.frame("Grupo" = Grupo, "p_value"= p_val_grup, "Bonferroni" = B_G,
                     "Benjamini_Hochberg" = BH_G),
          file="Grupos_AMOVA", row.names = F)
