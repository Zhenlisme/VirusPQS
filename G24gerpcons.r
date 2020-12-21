library(ggplot2)
library(ggpubr)

rm(list=ls())
gc()

setwd('E:/PQS_in_Plantvirus/All_viruses/FullSegment/RevisedSeq/GerpScore/')

class_info=read.csv2("../BasicStatistics/class_add_taxid.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$GenomeType=="circle-ssRNA","GenomeType"]="unknown"
class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'


higher_than_5=read.csv('alignedcount_higher_than_5.txt',header = F,stringsAsFactors = F)
head(higher_than_5)


higher_than_5=read.csv('alignedcount_higher_than_5.txt',header = F,stringsAsFactors = F)
head(higher_than_5)

gerpscore=read.csv2('G4gerpscore.csv',header = F,sep = ',',stringsAsFactors = F)
colnames(gerpscore)=c('ncnumber','pattern','start','end','seq','GerpScore','grunscore','loopscore')
gerpscore$GerpScore=as.numeric(gerpscore$GerpScore)
head(gerpscore)
levels(factor(gerpscore$pattern))

gerpscore=gerpscore[gerpscore$ncnumber %in% higher_than_5$V1,]

GNG_gerp=gerpscore[gerpscore$pattern=='GNG'|gerpscore$pattern=='CNC',c(1,6)]
GNG_gerp=aggregate(GNG_gerp$GerpScore,list(GNG_gerp$ncnumber),mean)
colnames(GNG_gerp)=c('ncnumber','score')
head(GNG_gerp)

G2_gerp=gerpscore[gerpscore$pattern=='GG'|gerpscore$pattern=='CC',c(1,6)]
G2_gerp=aggregate(G2_gerp$GerpScore,list(G2_gerp$ncnumber),mean)
colnames(G2_gerp)=c('ncnumber','score')
head(G2_gerp)
gerpcompare=merge(G2_gerp,GNG_gerp,by='ncnumber')
colnames(gerpcompare)=c('ncnumber','G2','GNG')

head(gerpcompare)
gathered_gerpcompare=tidyr::gather(gerpcompare,key='pattern',value=score,G2,GNG)

head(gathered_gerpcompare)
gathered_gerpcompare=merge(gathered_gerpcompare,class_info[,c(1,5,6)])

gathered_gerpcompare$euorpro=ifelse(gathered_gerpcompare$Host %in% c("plant","animal","fungi","protist"),"eukaryote",
                     ifelse(gathered_gerpcompare$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

ggplot(gathered_gerpcompare,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim=c(-20,10))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G2','GNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.005,textsize = 10,y_position = 9,
              test.args=list(paired=T,alternative='greater'),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  #stat_boxplot(geom = "errorbar",width=0.4,size=1)+
  #geom_boxplot(outlier.alpha = 1)+
  ylab('RS Score')+
  guides(fill=F)

ggsave("gerpcompare.pdf",width = 4,height =6)

ggplot(gathered_gerpcompare[gathered_gerpcompare$euorpro!="unclear",],aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim=c(-20,10))+
  facet_grid(.~euorpro,scales = "free_y")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G2','GNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.005,textsize = 10,y_position = 9,
              test.args=list(paired=T,alternative='greater'),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  ylab('RS Score')

ggsave("euorpro_gerpcompare.pdf",width = 6,height =6)

ggplot(gathered_gerpcompare,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  facet_wrap(.~GenomeType,ncol = 4)+
  coord_cartesian(ylim=c(-20,12))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G2','GNG')),test = 'wilcox.test',map_signif_level = T,
              tip_length = 0.005,textsize = 10,y_position = 8,
              test.args=list(paired=T,alternative='greater',exact=F),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  ylab('RS Score')

ggsave("gerpcompare_byGtp.pdf",width = 8,height =8)


ggplot(gathered_gerpcompare,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  facet_wrap(.~Host,ncol = 4)+
  coord_cartesian(ylim=c(-20,12))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G2','GNG')),test = 'wilcox.test',map_signif_level = T,
              tip_length = 0.005,textsize = 10,y_position = 8,
              test.args=list(paired=T,alternative='greater',exact=F),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  ylab('RS Score')

ggsave("gerpcompare_byHost.pdf",width = 8,height =8)


#############################################compare loops and gruns conservation##################################################

head(gerpscore)
gruns_and_loops=gerpscore[gerpscore$pattern=='GG'|gerpscore$pattern=='CC',c(1,7,8)]
gruns_and_loops$grunscore=as.numeric(gruns_and_loops$grunscore)
gruns_and_loops$loopscore=as.numeric(gruns_and_loops$loopscore)
head(gruns_and_loops)

gruns_and_loops=aggregate(gruns_and_loops[,c(2,3)],list(gruns_and_loops$ncnumber),mean)
colnames(gruns_and_loops)[1]=c('ncnumber')
head(gruns_and_loops)

head(class_info)

gruns_and_loops=tidyr::gather(gruns_and_loops,key='pattern',value='score',grunscore,loopscore)
gruns_and_loops=merge(gruns_and_loops,class_info[,c(1,5,6)])

head(gruns_and_loops)

gruns_and_loops$euorpro=ifelse(gruns_and_loops$Host %in% c("plant","animal","fungi","protist"),"eukaryote",
                               ifelse(gruns_and_loops$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim=c(-20,10))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('grunscore','loopscore')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.005,textsize = 10,
              y_position = 9,test.args=list(paired=T,alternative='greater'),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  ylab('RS Score for G2-PQS')+
  guides(fill=F)+
  scale_x_discrete(breaks=c('grunscore','loopscore'),labels=c('G-track','Loop'))

ggsave("gtrack_loop_conservation_gerp.pdf",width = 4,height =6)

gruns_and_loops$pattern=ifelse(gruns_and_loops$pattern=='loopscore','Loop','G-track')

head(gruns_and_loops)

ggplot(gruns_and_loops[gruns_and_loops$euorpro!="unclear",],aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim=c(-20,10))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  facet_grid(.~euorpro)+
  geom_signif(comparisons = list(c('G-track','Loop')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.005,textsize = 10,
              y_position = 9,test.args=list(paired=T,alternative='greater'),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
  ylab('RS Score for G2-PQS')

ggsave("euorpro_gtrack_loop_conservation_gerp.pdf",width = 6,height =6)


ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  facet_wrap(.~Host,ncol = 4)+
  coord_cartesian(ylim=c(-20,12))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G-track','Loop')),test = 'wilcox.test',y_position = 9,
              map_signif_level = T,tip_length = 0.01,textsize = 10,
              test.args=list(paired=T,alternative='greater',exact=F),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())+
  ylab('RS score for G2-PQS')


ggsave("gtrack_loop_conservation_byHost.pdf",width = 8,height=8)


ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  facet_wrap(.~GenomeType,ncol = 4)+
  coord_cartesian(ylim=c(-20,12))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G-track','Loop')),test = 'wilcox.test',y_position = 9,
              map_signif_level = T,tip_length = 0.01,textsize = 10,
              test.args=list(paired=T,alternative='greater',exact=F),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())+
  ylab('RS score for G2-PQS')


ggsave("gtrack_loop_conservation_byGtp.pdf",width = 8,height=8)



