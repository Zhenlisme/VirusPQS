rm(list=ls())
gc()

library(ggpubr)
library(ggplot2)
library(ggsci)

setwd("E:/PQS_in_Plantvirus/All_viruses/FullSegment/RevisedSeq/StructuralScore/")

class_info=read.csv2('../BasicStatistics/class_add_taxid.csv',sep = ',',
                     header = F,stringsAsFactors = F)
colnames(class_info)=c("AccessionNumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")
class_info[class_info$GenomeType=="circle-ssRNA","GenomeType"]="unknown"
class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','AccessionNumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'
head(class_info)

G4score=read.csv2("All_pattern_score.csv",sep = ",",header = T,stringsAsFactors = F)
head(G4score)
GNG=G4score[G4score$pattern=='GNG'|G4score$pattern=='CNC',c(1,6)]
GNG$pattern='GNG'
GG=G4score[G4score$pattern=='GG'|G4score$pattern=='CC',c(1,6)]
GG$pattern='GG'
GNG_simulate=G4score[G4score$pattern=='GNG_simulate'|G4score$pattern=='CNC_simulate',c(1,6)]
GNG_simulate$pattern='GNG_simulate'

G24score=rbind(GNG,GG,GNG_simulate)
G24score$G4HunterScore=as.numeric(G24score$G4HunterScore)
head(G24score)
class(G24score$G4HunterScore)

G24score=merge(G24score,class_info[,c('AccessionNumber','Taxid')])
head(G24score)


G24score_bytaxid=aggregate(G24score[,c(2)],list(G24score$pattern,G24score$Taxid),mean)
colnames(G24score_bytaxid)=c('type','Taxid','score')
head(G24score_bytaxid)
head(class_info)
G24score_bytaxid=tidyr::spread(G24score_bytaxid,key=type,value=score)

G24score_bytaxid=na.omit(G24score_bytaxid)


G24score_bytaxid=merge(G24score_bytaxid,unique(class_info[,c('GenomeType','Host','Taxid')]))
G24score_bytaxid=tidyr::gather(G24score_bytaxid,key=type,value=score,GG,GNG,GNG_simulate)

G24score_bytaxid$type=factor(G24score_bytaxid$type)
levels(G24score_bytaxid$type)=c('G2','GNG','S-G2')
G24score_bytaxid$type=factor(G24score_bytaxid$type,levels = c('GNG','G2','S-G2'),ordered = T)

head(G24score_bytaxid)

levels(G24score_bytaxid$type)

ggplot(G24score_bytaxid,aes(x=type,y=score,fill=type))+
  scale_fill_npg()+
  coord_cartesian(ylim = c(0,2.2))+
  ylab("G4Hunter Score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  stat_summary(fun = mean, geom = "point", aes(x=type,y=score),shape=16,size=2,color="blue",position = position_dodge(0.8))+
  guides(fill=F)+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("G2","GNG"),c('G2','S-G2'),c('GNG','S-G2')),
              family="serif",size=1,textsize = 10,
              map_signif_level = T,tip_length = 0.01,test = wilcox.test,y_position = c(1.7,1.8,2),
              test.args = list(paired=T,exact=F))


ggsave("G4HunterScore.pdf",width = 6,height = 6)

head(G24score_bytaxid)

ggplot(G24score_bytaxid,aes(x=type,y=score,fill=type))+
  scale_fill_npg()+
  coord_cartesian(ylim = c(0,3))+
  ylab("G4Hunter Score")+
  facet_wrap(.~GenomeType,ncol=4)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("G2","GNG"),c('G2','S-G2'),c('GNG','S-G2')),
              family="serif",size=0.8,textsize = 10,
              map_signif_level = T,tip_length = 0.02,test = wilcox.test,y_position = c(1.8,2.2,2.7),
              test.args = list(paired=T,exact=F))

ggsave("G4HunterScoreByGtp.pdf",width = 9,height = 9)

ggplot(G24score_bytaxid,aes(x=type,y=score,fill=type))+
  scale_fill_npg()+
  coord_cartesian(ylim = c(0,3))+
  ylab("G4Hunter Score")+
  facet_grid(.~Host)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c("G2","GNG"),c('G2','S-G2'),c('GNG','S-G2')),
              family="serif",size=0.8,textsize = 10,
              map_signif_level = T,tip_length = 0.02,test = wilcox.test,y_position = c(1.8,2.1,2.4),
              test.args = list(paired=T,exact=F))

ggsave("G4HunterScoreByHost.pdf",width = 9,height = 9)

