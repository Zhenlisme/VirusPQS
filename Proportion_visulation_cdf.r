rm(list=ls())
gc()

library(ggplot2)
library(ggpubr)



setwd("D:/毕业文件/李振数据备份/课题数据备份/RevisedSeq/random_g4hscore/")

class_info=read.csv("../classinfo_viroid_satellite.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'
head(class_info)

G2pvalue=read.csv("G24_pvalue_df.csv",sep = ',',header = T,stringsAsFactors = F)
G3pvalue=read.csv("G34_pvalue_df.csv",sep = ',',header = T,stringsAsFactors = F)


pvalue_df=rbind(G2pvalue,G3pvalue)
head(pvalue_df)

pvalue_df$explain=factor(pvalue_df$explain,levels = c('other','higher','lower'),ordered = T)
ggplot(pvalue_df[pvalue_df$euorpro!="unclear",],aes(x=euorpro,fill=explain))+
  facet_grid(.~G4type)+
  scale_fill_manual(values = c("#41B6E6","#FFB549",'#FF585D'))+
  geom_bar(width = .5,position = position_fill(0.5))+
  theme_set(theme_bw())+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_text(angle = -60,hjust = -.05),
        legend.title = element_blank(),
        legend.position = "top")+
  geom_text(aes(label=..count..),stat = "count",size=5,position = position_fill(0.5))+
  ylab("Proportion")

ggsave("shuffle_g4hscore_compare.pdf",width = 7,height = 7)



##For G2-PQS
table(pvalue_df[pvalue_df$G4type=='G2','explain'],
      pvalue_df[pvalue_df$G4type=='G2','euorpro'])

fisher.test(matrix(c(2144,2338,466,3360),ncol = 2,byrow = T),
            alternative = "greater")


##For G3-PQS
table(pvalue_df[pvalue_df$G4type=='G3','explain'],
      pvalue_df[pvalue_df$G4type=='G3','euorpro'])

fisher.test(matrix(c(1239,986,620,1068),ncol = 2,byrow = T),
            alternative = "greater")



##############compare G4Hunter score between pro and euro virus



G4score=read.csv2("All_pattern_score.csv",sep = ",",header = T,stringsAsFactors = F)
head(G4score)
real_g24=G4score[G4score$pattern=='GG'|G4score$pattern=='CC',c(1,6)]
real_g24$pattern='G2'
real_g34=G4score[G4score$pattern=='GGG'|G4score$pattern=='CCC',c(1,6)]
real_g34$pattern='G3'
real_g4=rbind(real_g24,real_g34)
head(real_g4)
real_g4=merge(real_g4,class_info[,c('ncnumber','Taxid')],by.x = "AccessionNumber",by.y = 'ncnumber')

head(class_info)

real_g4$G4HunterScore=as.numeric(real_g4$G4HunterScore)
head(real_g4)

real_g4=aggregate(real_g4[,c(2)],list(real_g4$Taxid,real_g4$pattern),mean)

colnames(real_g4)=c("Taxid","pattern","realscore")

head(class_info)
real_g4=merge(real_g4,unique(class_info[,c(6,9)]))

head(real_g4)
real_g4$europro=ifelse(real_g4$Host %in% c('fungi',"plant","animal","protist"),"eukaryote",
                       ifelse(real_g4$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))
head(real_g4)

unique(real_g4$pattern)

wilcox.test(real_g4[real_g4$pattern=='G2' & real_g4$europro=='eukaryote','realscore'],
            real_g4[real_g4$pattern=='G2' & real_g4$europro=='prokaryote','realscore'],
            alternative = 'greater')$p.value  ##2.542724e-240

wilcox.test(real_g4[real_g4$pattern=='G3' & real_g4$europro=='eukaryote','realscore'],
            real_g4[real_g4$pattern=='G3' & real_g4$europro=='prokaryote','realscore'],
            alternative = 'greater')$p.value  ##2.630734e-43


ggplot(real_g4[real_g4$europro!='unclear',],
       aes(realscore,color=europro))+
  facet_wrap(.~pattern)+
  coord_cartesian(xlim = c(0,4))+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.title = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  xlab('G4Hunter Score')+
  stat_ecdf(size=0.5)+
  ylab("Cumulative Distribution")+
  scale_color_manual(values = c("#E7B800","#00AFBB"))
  
ggsave('G4Hunter.pdf',width = 7,height = 6)

