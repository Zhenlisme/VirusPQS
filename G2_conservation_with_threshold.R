rm(list=ls())
gc()
library(ggplot2)
library(ggsci)
library(ggpubr)

setwd("D:/毕业文件/李振数据备份/课题数据备份/RevisedSeq/StructuralConservation/")
class_info=read.csv("../classinfo_viroid_satellite.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'

HigherThan5=read.csv('alignedcount_higher_than_5.txt',stringsAsFactors = F,sep = ',',header = F)
colnames(HigherThan5)='AccessionNumber'
head(HigherThan5)

GG_0=read.csv2("./SCI_0/GG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GG_0)=c("AccessionNumber","GGCI_0")
GG_5=read.csv2("./SCI_5/GG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GG_5)=c("AccessionNumber","GGCI_5")
GG_10=read.csv2("./SCI_10/GG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GG_10)=c("AccessionNumber","GGCI_10")
GGsci=cbind(GG_0,GG_5,GG_10)[,c(1,2,4,6)]
GGsci=tidyr::gather(GGsci,key=pattern,value=con,-AccessionNumber)
head(GGsci)

GNG_0=read.csv2("./SCI_0/GNG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GNG_0)=c("AccessionNumber","GNGCI_0")
GNG_5=read.csv2("./SCI_5/GNG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GNG_5)=c("AccessionNumber","GNGCI_5")
GNG_10=read.csv2("./SCI_10/GNG.g4sci",sep = "\t",stringsAsFactors = F,header = F)[,c(1,6)]
colnames(GNG_10)=c("AccessionNumber","GNGCI_10")
GNGsci=cbind(GNG_0,GNG_5,GNG_10)[,c(1,2,4,6)]
GNGsci=tidyr::gather(GNGsci,key=pattern,value=con,-AccessionNumber)
head(GNGsci)

conservation_df=rbind(GGsci,GNGsci)
conservation_df$con=as.numeric(conservation_df$con)
head(conservation_df)
conservation_df=conservation_df[conservation_df$AccessionNumber %in% HigherThan5$AccessionNumber,]

mean_conservation=aggregate(conservation_df$con,list(conservation_df$AccessionNumber,conservation_df$pattern),mean)
colnames(mean_conservation)=c("AccessionNumber","pattern","Conservation")
head(mean_conservation)

spreaded_mean_conservation=tidyr::spread(mean_conservation,key=pattern,value=Conservation)
head(spreaded_mean_conservation)

GNG=na.omit(spreaded_mean_conservation)
GNG=tidyr::gather(GNG,key=pattern,value=Con,-AccessionNumber)

marker=as.data.frame(do.call(rbind,strsplit(GNG$pattern,split = "_")))
GNG$extension=marker$V2
GNG$pattern=marker$V1

head(GNG)
GNG=merge(GNG,class_info[,c(1,5,6)],by.x='AccessionNumber',by.y = 'ncnumber')
GNG$euorpro=ifelse(GNG$Host %in% c("plant","animal","alga","fungi","protist"),"eukaryote",
                   ifelse(GNG$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))
head(GNG)
levels(factor(GNG$extension))
GNG$extension=factor(GNG$extension,levels = c('0','5','10'),
                     labels = c('extension=0','extension=5','extension=10'))

ggplot(GNG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  #coord_cartesian(ylim = c(0,1.2))+
  facet_grid(.~extension)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  ylab("Pattern Conservation Value")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("GGCI","GNGCI")),family="serif",size=1,map_signif_level = T,textsize = 10,
              y_position = 0.92,tip_length = 0.02,test = wilcox.test,test.args = list(paired=T,alternative="greater"))+
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(size=26,face = "bold",colour = "black",family="serif"),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  guides(fill=F)+
  theme(strip.text = element_text(face = 'bold',size=rel(2),color = 'black',family = 'serif'))+
  scale_x_discrete(breaks=c('GGCI','GNGCI'),labels=c('G2','GNG'))

ggsave("G2-GNG-ConservationIndexCompare.pdf",width = 10,height =6)



head(GNG)
ggplot(GNG[GNG$euorpro!="unclear",],aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim = c(0,1.2))+
  facet_grid(euorpro~extension)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(strip.text.x = element_text(size = 26),
        strip.text.y = element_text(size = 26))+
  ylab("Pattern Conservation Value")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("GGCI","GNGCI")),family="serif",size=1,map_signif_level = T,textsize = 10,
              y_position = 1.05,tip_length = 0.02,test = wilcox.test,test.args = list(paired=T,alternative="greater"))+
  theme(axis.title.x=element_blank(),
        axis.title.y = element_text(size=26,face = "bold",colour = "black",family="serif"),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  guides(fill=F)+
  theme(strip.text = element_text(face = 'bold',size=rel(2),color = 'black',family = 'serif'))+
  scale_x_discrete(breaks=c('GGCI','GNGCI'),labels=c('G2','GNG'))+
  scale_y_continuous(breaks = c(0,0.5,1))

ggsave("euorpro_G2-GNG_ConservationIndexCompare.pdf",width = 8,height =8)

#############################################################################

ggplot(GNG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim = c(0,1.2))+
  facet_grid(extension~GenomeType)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  ylab("Pattern Conservation Value")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("GGCI","GNGCI")),family="serif",size=1,map_signif_level = T,textsize = 10,
              y_position = 1.1,tip_length = 0.02,test = wilcox.test,test.args = list(paired=T,alternative="greater",exact=F))+
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=26,face = "bold",colour = "black",family="serif"),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  guides(fill=F)+
  theme(strip.text = element_text(face = 'bold',size=rel(2),color = 'black',family = 'serif'))+
  scale_x_discrete(breaks=c('GGCI','GNGCI'),labels=c('G2','GNG'))

ggsave("G2-GNG-GTP.pdf",width = 15,height =15)

ggplot(GNG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim = c(0,1.2))+
  facet_grid(extension~Host)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  ylab("Pattern Conservation Value")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("GGCI","GNGCI")),family="serif",size=1,map_signif_level = T,textsize = 10,
              y_position = 1.1,tip_length = 0.02,test = wilcox.test,test.args = list(paired=T,alternative="greater",exact=F))+
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y = element_text(size=26,face = "bold",colour = "black",family="serif"),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  guides(fill=F)+
  theme(strip.text = element_text(face = 'bold',size=rel(2),color = 'black',family = 'serif'))+
  scale_x_discrete(breaks=c('GGCI','GNGCI'),labels=c('G2','GNG'))

ggsave("G2-GNG-Host.pdf",width = 15,height =15)



