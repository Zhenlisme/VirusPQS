rm(list=ls())
gc()
library(ggplot2)
library(ggsci)
library(ggpubr)

setwd("E:/PQS_in_Plantvirus/All_viruses/FullSegment/RevisedSeq/StructuralConservation/")

class_info=read.csv2("../BasicStatistics/class_add_taxid.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")
class_info[class_info$GenomeType=="circle-ssRNA","GenomeType"]="unknown"
class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'

HigherThan5=read.csv('alignedcount_higher_than_5.txt',stringsAsFactors = F,sep = ',',header = F)
colnames(HigherThan5)='AccessionNumber'
head(HigherThan5)

GGG_0=read.csv2("./SCI_0/GGG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GGG_0)=c("AccessionNumber","GGGCI_0")
GGG_5=read.csv2("./SCI_5/GGG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GGG_5)=c("AccessionNumber","GGGCI_5")
GGG_10=read.csv2("./SCI_10/GGG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GGG_10)=c("AccessionNumber","GGGCI_10")

GGGsci=cbind(GGG_0,GGG_5,GGG_10)[,c(1,2,4,6)]
GGGsci=tidyr::gather(GGGsci,key=pattern,value=con,-AccessionNumber)
head(GGGsci)

GGNG_0=read.csv2("./SCI_0/GGNG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GGNG_0)=c("AccessionNumber","GGNGCI_0")
GGNG_5=read.csv2("./SCI_5/GGNG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GGNG_5)=c("AccessionNumber","GGNGCI_5")
GGNG_10=read.csv2("./SCI_10/GGNG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GGNG_10)=c("AccessionNumber","GGNGCI_10")

GGNGsci=cbind(GGNG_0,GGNG_5,GGNG_10)[,c(1,2,4,6)]
GGNGsci=tidyr::gather(GGNGsci,key=pattern,value=con,-AccessionNumber)
head(GGNGsci)

GNGG_0=read.csv2("./SCI_0/GNGG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GNGG_0)=c("AccessionNumber","GNGGCI_0")
GNGG_5=read.csv2("./SCI_5/GNGG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GNGG_5)=c("AccessionNumber","GNGGCI_5")
GNGG_10=read.csv2("./SCI_10/GNGG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GNGG_10)=c("AccessionNumber","GNGGCI_10")

GNGGsci=cbind(GNGG_0,GNGG_5,GNGG_10)[,c(1,2,4,6)]
GNGGsci=tidyr::gather(GNGGsci,key=pattern,value=con,-AccessionNumber)
head(GNGGsci)

conservation_df=rbind(GGGsci,GGNGsci,GNGGsci)
conservation_df$con=as.numeric(conservation_df$con)
conservation_df=conservation_df[conservation_df$AccessionNumber %in% HigherThan5$AccessionNumber,]
head(conservation_df)
class(conservation_df$con)
mean_conservation=aggregate(conservation_df$con,list(conservation_df$AccessionNumber,conservation_df$pattern),mean)
colnames(mean_conservation)=c("AccessionNumber","pattern","Conservation")
head(mean_conservation)

spreaded_mean_conservation=tidyr::spread(mean_conservation,key=pattern,value=Conservation)
head(spreaded_mean_conservation)

GGNG=na.omit(spreaded_mean_conservation[,c(1,2,3,4,5,6,7)])
GNGG=na.omit(spreaded_mean_conservation[,c(1,2,3,4,8,9,10)])

GGNG=tidyr::gather(GGNG,key=pattern,value=Con,-AccessionNumber)
marker=as.data.frame(do.call(rbind,strsplit(GGNG$pattern,split = "_")))
GGNG$extension=marker$V2
GGNG$pattern=marker$V1

GNGG=tidyr::gather(GNGG,key=pattern,value=Con,-AccessionNumber)
marker=as.data.frame(do.call(rbind,strsplit(GNGG$pattern,split = "_")))
GNGG$extension=marker$V2
GNGG$pattern=marker$V1
head(GNGG)

GNGG=merge(GNGG,class_info[,c(1,5,6)],by.x='AccessionNumber',by.y = 'ncnumber')
GNGG$euorpro=ifelse(GNGG$Host %in% c("plant","animal","fungi","protist"),"eukaryote",
                     ifelse(GNGG$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))
levels(GNGG$extension)=c('extension=0','extension=10','extension=5')

GGNG=merge(GGNG,class_info[,c(1,5,6)],by.x='AccessionNumber',by.y = 'ncnumber')
GGNG$euorpro=ifelse(GGNG$Host %in% c("plant","animal","fungi","protist"),"eukaryote",
                    ifelse(GGNG$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))
levels(GGNG$extension)=c('extension=0','extension=10','extension=5')

head(GNGG)

ggplot(GNGG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim = c(0,1.2))+
  facet_grid(extension~GenomeType)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  ylab("Pattern Conservation Value")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("GGGCI","GNGGCI")),family="serif",size=1,map_signif_level = T,textsize = 10,
              y_position = 1.1,tip_length = 0.02,test = wilcox.test,test.args = list(paired=T,alternative="less"))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=24,face = "bold",colour = "black",family="serif"),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  guides(fill=F)+
  theme(strip.text = element_text(face = 'bold',size=rel(1.5),color = 'black',family = 'serif'))+
  scale_x_discrete(breaks=c('GGGCI','GNGGCI'),labels=c('G3','GNGG'))


ggplot(GNGG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim = c(0,1.2))+
  facet_grid(.~extension)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  ylab("Pattern Conservation Value")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("GGGCI","GNGGCI")),family="serif",size=1,map_signif_level = T,textsize = 10,
              y_position = 1.1,tip_length = 0.02,test = wilcox.test,test.args = list(paired=T,alternative="less"))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=24,face = "bold",colour = "black",family="serif"),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  guides(fill=F)+
  theme(strip.text = element_text(face = 'bold',size=rel(1.5),color = 'black',family = 'serif'))+
  scale_x_discrete(breaks=c('GGGCI','GNGGCI'),labels=c('G3','GNGG'))

ggsave("G3-GNGG-ConservationIndexCompare.pdf",width = 10,height =6)

ggplot(GGNG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim = c(0,1.2))+
  facet_grid(.~extension)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  ylab("Pattern Conservation Value")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("GGGCI","GGNGCI")),family="serif",size=1,map_signif_level = T,textsize = 10,
              y_position = 1.1,tip_length = 0.02,test = wilcox.test,test.args = list(paired=T,alternative="less"))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=24,face = "bold",colour = "black",family="serif"),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  guides(fill=F)+
  theme(strip.text = element_text(face = 'bold',size=rel(1.5),color = 'black',family = 'serif'))+
  scale_x_discrete(breaks=c('GGGCI','GGNGCI'),labels=c('G3','GGNG'))

ggsave("G3-GGNG-ConservationIndexCompare.pdf",width = 10,height =6)


ggplot(GNGG[GNGG$euorpro!="unclear",],aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim = c(0,1.2))+
  facet_grid(euorpro~extension)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  ylab("Pattern Conservation Value")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("GGGCI","GNGGCI")),family="serif",size=1,map_signif_level = T,textsize = 10,
              y_position = 1.05,tip_length = 0.02,test = wilcox.test,test.args = list(paired=T,alternative="less"))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=24,face = "bold",colour = "black",family="serif"),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  guides(fill=F)+
  theme(strip.text = element_text(face = 'bold',size=rel(1.5),color = 'black',family = 'serif'))+
  scale_x_discrete(breaks=c('GGGCI','GNGGCI'),labels=c('G3','GNGG'))

ggsave("euorpro_G3-GNGG_ConservationIndexCompare.pdf",width = 10,height =6)

ggplot(GGNG[GGNG$euorpro!="unclear",],aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim = c(0,1.2))+
  facet_grid(euorpro~extension)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  ylab("Pattern Conservation Value")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("GGGCI","GGNGCI")),family="serif",size=1,map_signif_level = T,textsize = 10,
              y_position = 1.05,tip_length = 0.02,test = wilcox.test,test.args = list(paired=T,alternative="less"))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=24,face = "bold",colour = "black",family="serif"),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  guides(fill=F)+
  theme(strip.text = element_text(face = 'bold',size=rel(1.5),color = 'black',family = 'serif'))+
  scale_x_discrete(breaks=c('GGGCI','GGNGCI'),labels=c('G3','GGNG'))

ggsave("euorpro_G3-GGNG_ConservationIndexCompare.pdf",width = 10,height =6)

