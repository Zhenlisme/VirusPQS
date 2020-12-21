rm(list=ls())
gc()

library(vioplot)
library(ggplot2)
library(ggpubr)
library(ggsci)


setwd("E:/PQS_in_Plantvirus/All_viruses/FullSegment/RevisedSeq/BasicStatistics/")

class_info=read.csv2("class_add_taxid.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$GenomeType=="circle-ssRNA","GenomeType"]="unknown"
class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'

sequence_length=read.csv2("sequence_length.csv",sep=",",header = T,stringsAsFactors = F)
sequence_length=unique(sequence_length[,c(1,4)])
sequence_length$genomelength=2*sequence_length$genomelength
sequence_length=merge(sequence_length,class_info[,c('ncnumber','Taxid','Host','GenomeType')])
sequence_length=aggregate(sequence_length[,c('genomelength')],list(sequence_length$Taxid),sum)
colnames(sequence_length)=c('taxid','genomelength')
head(sequence_length)

Allpattern_df=read.csv2("long_pattern_count.csv",sep=",",header = T,stringsAsFactors = F)
Allpattern_df$G2=Allpattern_df$GG_count+Allpattern_df$CC_count
Allpattern_df$GNG=Allpattern_df$GNG_count+Allpattern_df$CNC_count
Allpattern_df$G3=Allpattern_df$GGG_count+Allpattern_df$CCC_count
Allpattern_df$GNGG=Allpattern_df$GNGG_count+Allpattern_df$CCNC_count
Allpattern_df$GGNG=Allpattern_df$GGNG_count+Allpattern_df$CNCC_count
Allpattern_df=Allpattern_df[,c("ncnumber","G2","GNG","G3","GNGG","GGNG")]
Allpattern_df=merge(Allpattern_df,class_info[,c('ncnumber','Taxid','Host','GenomeType')])
Allpattern_df=aggregate(Allpattern_df[,c(2:6)],list(Allpattern_df$Taxid),sum)
colnames(Allpattern_df)[1]='taxid'
Allpattern_df=merge(Allpattern_df,sequence_length)
Allpattern_df[,c(2:6)]=1000*Allpattern_df[,c(2:6)]/Allpattern_df$genomelength
head(Allpattern_df)

wilcox.test(Allpattern_df$G2,Allpattern_df$GNG,'greater')
wilcox.test(Allpattern_df$G3,Allpattern_df$GGNG,'less')
wilcox.test(Allpattern_df$G3,Allpattern_df$GNGG,'less')
#################################################GNG compare############################
random_shuffling=read.csv2("Allvirus_shuffled_count.csv",stringsAsFactors = F,header = T,sep = ",")

random_shuffling$G2=random_shuffling$GG_count+random_shuffling$CC_count
random_shuffling$GNG=random_shuffling$GNG_count+random_shuffling$CNC_count
random_shuffling$G3=random_shuffling$GGG_count+random_shuffling$CCC_count
random_shuffling$GNGG=random_shuffling$GNGG_count+random_shuffling$CCNC_count
random_shuffling$GGNG=random_shuffling$GGNG_count+random_shuffling$CNCC_count

random_shuffling=random_shuffling[,c("ncnumber","G2","GNG","G3","GNGG","GGNG")]
head(random_shuffling)
merged_taxid_shuffle_df=merge(random_shuffling,class_info[,c("ncnumber","Taxid")])
merged_taxid_shuffle_df$times=rep(1:1000,10394)
merged_taxid_shuffle_df=aggregate(merged_taxid_shuffle_df[,c(2:6)],
                                  list(merged_taxid_shuffle_df$Taxid,merged_taxid_shuffle_df$times),sum)

colnames(merged_taxid_shuffle_df)[1:2]=c('taxid','times')
merged_taxid_shuffle_df=merge(merged_taxid_shuffle_df,sequence_length)
merged_taxid_shuffle_df[,c(3,4,5,6,7)]=1000*merged_taxid_shuffle_df[,c(3,4,5,6,7)]/merged_taxid_shuffle_df$genomelength

merged_taxid_shuffle_df=merged_taxid_shuffle_df[,c(1,3:7)]
head(merged_taxid_shuffle_df)
gathered_shuffled_df=tidyr::gather(merged_taxid_shuffle_df,key=pattern,value=density,-taxid)
head(gathered_shuffled_df)
summary_shuffle=aggregate(gathered_shuffled_df$density,list(gathered_shuffled_df$taxid,gathered_shuffled_df$pattern),
                          function(x) c(mean(x),sd(x)))

summary_shuffle[,c(4,5)]=as.data.frame(summary_shuffle$x)
colnames(summary_shuffle)=c('taxid','pattern','x','density','sd')
summary_shuffle=summary_shuffle[,c(1,2,4,5)]
Allpattern_df=tidyr::gather(Allpattern_df[,1:6],key=pattern,value=density,-taxid)
names(summary_shuffle)[3]='shuffle_desnity'
names(Allpattern_df)[3]='real_density'
head(summary_shuffle)

head(Allpattern_df)

shuffle_and_real=merge(Allpattern_df,summary_shuffle,by=c('taxid','pattern'))
head(shuffle_and_real)
shuffle_and_real$zscore=(shuffle_and_real$real_density-shuffle_and_real$shuffle_desnity)/shuffle_and_real$sd
head(class_info)
head(shuffle_and_real)
write.table(shuffle_and_real,"zscore_total.csv",col.names = T,row.names = F,quote = F,sep = ',')
write.table(tidyr::spread(shuffle_and_real[,c(1,2,6)],key=pattern,value=zscore),"zscore.csv",col.names = T,row.names = F,quote = F,sep = ',')

shuffle_and_real=merge(shuffle_and_real,unique(class_info[,c('GenomeType',"Host",'Taxid')]),by.x = 'taxid',by.y = 'Taxid',all.x = T)
G2_compare=shuffle_and_real[shuffle_and_real$pattern=='G2'|shuffle_and_real$pattern=='GNG',]
G3_compare=shuffle_and_real[shuffle_and_real$pattern=='G3'|shuffle_and_real$pattern=='GGNG'|shuffle_and_real$pattern=='GNGG',]
G3_compare[is.na(G3_compare)]=0

droped_sets=unique(G3_compare[is.infinite(G3_compare$zscore),'taxid'])
G3_compare=G3_compare[which(!G3_compare$taxid%in%droped_sets),]
G3_compare[G3_compare$zscore==max(G3_compare$zscore),]

aggregate(G2_compare$real_density,list(G2_compare$pattern),mean)
aggregate(G3_compare$real_density,list(G3_compare$pattern),mean)
head(G2_compare)


G2_compare$europro=ifelse(G2_compare$Host %in% c('fungi',"plant","animal","protist"),"eukaryote",
                          ifelse(G2_compare$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

ggplot(G2_compare,aes(x=pattern,y=zscore,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim = c(-10,15))+
  ylab("z-score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  guides(fill=F)+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  geom_signif(comparisons = list(c("G2","GNG")),family="serif",size=1,textsize = 10,
              map_signif_level = T,y_position = 14,tip_length = 0.003,test = wilcox.test,
              test.args = list(alternative='greater',paired=T))+
  stat_boxplot(geom = "errorbar",width=0.2,size=1)+
  geom_boxplot(width=0.5,outlier.alpha = 0)
  #stat_summary(fun = mean, geom = "point", aes(x=pattern,y=zscore),shape=16,size=2,color="blue",position = position_dodge(0.8))

ggsave("G2PQS_zscore.pdf",width = 6,height =6)

head(G2_compare)

ggplot(G2_compare[G2_compare$europro!="unclear",],aes(x=pattern,y=zscore,fill=pattern))+
  scale_fill_npg()+
  facet_wrap(.~europro,ncol=4)+
  coord_cartesian(ylim = c(-12,12))+
  ylab("z-score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        legend.position = 'top',
        legend.title = element_blank())+
  geom_signif(comparisons = list(c("G2","GNG")),family="serif",size=1,textsize = 10,
              map_signif_level = T,y_position = 11,tip_length = 0.002,test = wilcox.test,
              test.args = list(alternative='greater',paired=T))+
  stat_boxplot(geom = "errorbar",width=0.2,size=1)+
  geom_boxplot(width=0.5,outlier.alpha = 0)

ggsave("europro_G2PQS_zscore.pdf",width = 6,height =6)


ggplot(G2_compare,aes(x=pattern,y=zscore,fill=pattern))+
  scale_fill_npg()+
  facet_wrap(.~GenomeType,ncol=4)+
  coord_cartesian(ylim = c(-10,20))+
  ylab("z-score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        legend.position = 'top',
        legend.title = element_blank())+
  geom_signif(comparisons = list(c("G2","GNG")),family="serif",size=1,textsize = 10,
              map_signif_level = T,y_position = 15,tip_length = 0.01,test = wilcox.test,
              test.args = list(alternative='greater',paired=T))+
  stat_boxplot(geom = "errorbar",width=0.2,size=1)+
  geom_boxplot(width=0.5,outlier.alpha = 0)

ggsave("G2PQS_zscore_Gtp.pdf",width = 9,height =9)

head(G2_compare)
ggplot(G2_compare,aes(x=pattern,y=zscore,fill=pattern))+
  scale_fill_npg()+
  facet_grid(.~Host)+
  coord_cartesian(ylim = c(-10,50))+
  ylab("z-score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        legend.position = 'top',
        legend.title = element_blank())+
  geom_signif(comparisons = list(c("G2","GNG")),family="serif",size=1,textsize = 10,
              map_signif_level = T,y_position = 45,tip_length = 0.005,test = wilcox.test,
              test.args = list(alternative='greater',paired=T))+
  stat_boxplot(geom = "errorbar",width=0.2,size=1)+
  geom_boxplot(width=0.5,outlier.alpha = 0)

ggsave("G2PQS_zscore_Host.pdf",width = 9,height =9)

G3violion=G3_compare[,c(1,2,6,7)]
G3violion=tidyr::spread(G3violion,key=pattern,value=zscore)
G3_compare$europro=ifelse(G3_compare$Host %in% c('fungi',"plant","animal","protist"),"eukaryote",
                          ifelse(G3_compare$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

head(G3_compare)

ggplot(G3_compare,aes(x=pattern,y=zscore,fill=pattern))+
  #geom_violin(trim=T,color="white")+
  #geom_boxplot(width=0.2,outlier.alpha = 0)+
  scale_fill_npg()+
  coord_cartesian(ylim = c(-5,10))+
  ylab("z-score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  guides(fill=F)+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  geom_signif(comparisons = list(c("G3","GGNG"),c("G3","GNGG")),family="serif",size=1,textsize = 10,
              map_signif_level = T,y_position = c(7,8.5),tip_length = 0.0005,test = wilcox.test,
              test.args = list(alternative='less',paired=T))+
  stat_boxplot(geom = "errorbar",width=0.2,size=1)+
  geom_boxplot(width=0.5,outlier.alpha = 0)
  
ggsave("G3PQS_zscore.pdf",width = 6,height =6)

ggplot(G3_compare[G3_compare$europro!="unclear",],aes(x=pattern,y=zscore,fill=pattern))+
  scale_fill_npg()+
  facet_grid(.~europro)+
  coord_cartesian(ylim = c(-5,7))+
  ylab("z-score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.title.x = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  geom_signif(comparisons = list(c("G3","GGNG"),c("G3","GNGG")),family="serif",size=1,textsize = 10,
              map_signif_level = T,y_position = c(5.5,6.5),tip_length = 0.0003,test = wilcox.test,
              test.args = list(paired=T,exact=F))+
  stat_boxplot(geom = "errorbar",width=0.2,size=1)+
  geom_boxplot(width=0.5,outlier.alpha = 0)

ggsave("europro_G3PQS_zscore.pdf",width = 6,height =6)


ggplot(G3_compare,aes(x=pattern,y=zscore,fill=pattern))+
  scale_fill_npg()+
  facet_grid(.~GenomeType)+
  coord_cartesian(ylim = c(-5,12))+
  ylab("z-score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.title.x = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  geom_signif(comparisons = list(c("G3","GGNG"),c("G3","GNGG")),family="serif",size=1,textsize = 10,
              map_signif_level = T,y_position = c(10,11.5),tip_length = 0.0003,test = wilcox.test,
              test.args = list(alternative='less',paired=T,exact=F))+
  stat_boxplot(geom = "errorbar",width=0.2,size=1)+
  geom_boxplot(width=0.5,outlier.alpha = 0)

ggsave("G3PQS_zscore_Gtp.pdf",width = 15,height =10)

ggplot(G3_compare,aes(x=pattern,y=zscore,fill=pattern))+
  scale_fill_npg()+
  facet_grid(.~Host)+
  coord_cartesian(ylim = c(-5,22))+
  ylab("z-score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.title.x = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  geom_signif(comparisons = list(c("G3","GGNG"),c("G3","GNGG")),family="serif",size=1,textsize = 10,
              map_signif_level = T,y_position = c(19,21),tip_length = 0.0003,test = wilcox.test,
              test.args = list(alternative='less',paired=T,exact=F))+
  stat_boxplot(geom = "errorbar",width=0.2,size=1)+
  geom_boxplot(width=0.5,outlier.alpha = 0)


ggsave("G3PQS_zscore_Host.pdf",width = 15,height =10)

