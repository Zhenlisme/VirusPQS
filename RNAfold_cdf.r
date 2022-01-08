rm(list=ls())
gc()

library(ggplot2)
library(ggpubr)
library(ggsci)
setwd('D:/毕业文件/李振数据备份/课题数据备份/RevisedSeq/RNAfold/')
class_info=read.csv("../classinfo_viroid_satellite.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'

##G2PQS

G2PQS_mfe=read.csv2('G2PQS.mfe.csv',sep = ',',stringsAsFactors = F,header = F)
colnames(G2PQS_mfe)=c('ncnumber','PQS','marker','mfe')
G2PQS_mfe$marker=ifelse(G2PQS_mfe$marker==1,'G4','notG4')

G2PQS_mfe=merge(G2PQS_mfe,class_info[,c(1,9)])


summary_count_G2=aggregate(G2PQS_mfe$Taxid,list(G2PQS_mfe$Taxid,G2PQS_mfe$marker),length)

colnames(summary_count_G2)=c('Taxid','Marker','Count')
head(summary_count_G2)
summary_count_G2=tidyr::spread(summary_count_G2,key=Marker,value=Count)
head(summary_count_G2)

summary_count_G2[is.na(summary_count_G2)]=0
summary_count_G2$G4pro=summary_count_G2$G4/(summary_count_G2$G4+summary_count_G2$notG4)

summary_count_G2=merge(summary_count_G2,unique(class_info[,c(5,6,9)]))
head(summary_count_G2)
summary_count_G2[summary_count_G2$Taxid=='1000373',]


summary_count_G2$europro=ifelse(summary_count_G2$Host %in% c('fungi',"plant","animal","protist"),"eukaryote",
                          ifelse(summary_count_G2$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))
head(summary_count_G2)
summary_count_G2$G4Type='G2'

##G3PQS

G3PQS_mfe=read.csv2('G3PQS.mfe.csv',sep = ',',stringsAsFactors = F,header = F)
colnames(G3PQS_mfe)=c('ncnumber','PQS','marker','mfe')
G3PQS_mfe$marker=ifelse(G3PQS_mfe$marker==1,'G4','notG4')
head(G3PQS_mfe)
G3PQS_mfe=merge(G3PQS_mfe,class_info[,c(1,9)])
summary_count_G3=aggregate(G3PQS_mfe$Taxid,list(G3PQS_mfe$Taxid,G3PQS_mfe$marker),length)

colnames(summary_count_G3)=c('Taxid','Marker','Count')
summary_count_G3=tidyr::spread(summary_count_G3,key=Marker,value=Count)
summary_count_G3[is.na(summary_count_G3)]=0
summary_count_G3$G4pro=summary_count_G3$G4/(summary_count_G3$G4+summary_count_G3$notG4)
summary_count_G3=merge(summary_count_G3,unique(class_info[,c(5,6,9)]))
head(summary_count_G3)

summary_count_G3$europro=ifelse(summary_count_G3$Host %in% c('fungi',"plant","animal","protist"),"eukaryote",
                             ifelse(summary_count_G3$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))
summary_count_G3$G4Type="G3"

head(summary_count_G3)
head(summary_count_G2)

Total_summary=rbind(summary_count_G2,summary_count_G3)

head(Total_summary)

wilcox.test(Total_summary[Total_summary$G4Type=='G2' & Total_summary$europro=='eukaryote',"G4pro"],
            Total_summary[Total_summary$G4Type=='G2' & Total_summary$europro=='prokaryote',"G4pro"],
            alternative = 'greater')$p.value  # 3.515955e-17

wilcox.test(Total_summary[Total_summary$G4Type=='G3' & Total_summary$europro=='eukaryote',"G4pro"],
            Total_summary[Total_summary$G4Type=='G3' & Total_summary$europro=='prokaryote',"G4pro"],
            alternative = 'greater')$p.value  # 1.119414e-25

ggplot(Total_summary[Total_summary$europro!='unclear',],
       aes(x=europro,y=G4pro,fill=europro))+
  facet_wrap(.~G4Type,scales = 'free')+
  coord_cartesian(ylim = c(0,1.3))+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('eukaryote','prokaryote')),test = 'wilcox.test',
              map_signif_level = T,tip_length = c(0.02),textsize = 10,y_position = 1.1,
              test.args=list(alternative='greater',exact=F),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  ylab('Proportion')+
  stat_boxplot(geom = "errorbar",width=0.2,size=1)+
  geom_boxplot(width=0.5,outlier.alpha = 0)+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))

ggsave('RNAfold_proportion.pdf',width = 6,height = 6)


####mfe compare

G2PQS_mfe$G4Type="G2"
G3PQS_mfe$G4Type="G3"
PQS_mfe=rbind(G2PQS_mfe[G2PQS_mfe$marker=='G4',c(4,5,6)],
              G3PQS_mfe[G3PQS_mfe$marker=='G4',c(4,5,6)])
PQS_mfe=merge(PQS_mfe,unique(class_info[,c('Host','Taxid')]))
PQS_mfe$mfe=as.numeric(PQS_mfe$mfe)
PQS_mfe$europro=ifelse(PQS_mfe$Host %in% c('fungi',"plant","animal","protist"),"eukaryote",
                                ifelse(PQS_mfe$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))
head(PQS_mfe)

options(scipen = 10,digits = 22)
wilcox.test(PQS_mfe[PQS_mfe$G4Type=='G2' & PQS_mfe$europro=='eukaryote','mfe'],
            PQS_mfe[PQS_mfe$G4Type=='G2' & PQS_mfe$europro=='prokaryote','mfe'],
            alternative = 'less')$p.value   ##0


wilcox.test(PQS_mfe[PQS_mfe$G4Type=='G3' & PQS_mfe$europro=='eukaryote','mfe'],
            PQS_mfe[PQS_mfe$G4Type=='G3' & PQS_mfe$europro=='prokaryote','mfe'],
            alternative = 'less')$p.value  ## 1.196696e-146

ggplot(PQS_mfe[PQS_mfe$europro!='unclear',],
       aes(mfe,color=europro))+
  facet_wrap(.~G4Type,scales='free')+
  coord_cartesian(xlim = c(-50,10))+
  #geom_boxplot(width=0.2,outlier.alpha = 0)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.position = "top",
        legend.title = element_blank())+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  xlab('Minimum Free Energy (kcal/mol)')+
  ylab("Cumulative Distribution")+
  stat_ecdf(size=0.5)+
  scale_color_manual(values = c("#E7B800","#00AFBB"))



ggsave('RNAfold_mfe.pdf',width = 7,height = 6)







