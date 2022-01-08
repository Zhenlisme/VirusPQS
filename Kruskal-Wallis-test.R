rm(list = ls())
gc()

install.packages("spdep")
install.packages("pgirmess")


library(spdep)
library(pgirmess)


setwd("D:/毕业文件/李振数据备份/课题数据备份/RevisedSeq/BasicStatistics/")


class_info=read.csv("../classinfo_viroid_satellite.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")
class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'

G2_density=read.csv("G2compare.csv",sep=",",header = T,stringsAsFactors = F)
G2_density=G2_density[G2_density$pattern=="G2",c(1,3)]
G2_density=merge(G2_density,unique(class_info[,c(3,5,6,9)]),by.x="taxid",by.y = "Taxid")
head(G2_density)

class(G2_density$GenomeType)

kruskal.test(real_density~GenomeType,data=G2_density)$p.value  ## 1.631712e-90

kruskal.test(real_density~Host,data=G2_density)$p.value  ## 1.279821e-20

G3_density=read.csv("G3compare.csv",sep=",",header = T,stringsAsFactors = F)
G3_density=G3_density[G3_density$pattern=="G3",c(1,3)]
G3_density=merge(G3_density,unique(class_info[,c(3,5,6,9)]),by.x="taxid",by.y = "Taxid")


kruskal.test(real_density~GenomeType,data=G3_density)$p.value  ## 1.519525e-202

kruskal.test(real_density~Host,data=G3_density)$p.value  ## 2.564951e-119





