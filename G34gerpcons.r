library(ggplot2)
library(ggpubr)
library(ggsci)
rm(list=ls())
gc()

setwd('E:/课题数据备份/RevisedSeq/GerpScore/')

class_info=read.csv("../classinfo_viroid_satellite.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'


higher_than_5=read.csv('alignedcount_higher_than_5.txt',header = F,stringsAsFactors = F)
head(higher_than_5)

gerpscore=read.csv2('G4gerpscore.csv',header = F,sep = ',',stringsAsFactors = F)
colnames(gerpscore)=c('ncnumber','pattern','start','end','seq','GerpScore','grunscore','loopscore')
gerpscore$GerpScore=as.numeric(gerpscore$GerpScore)
head(gerpscore)
levels(factor(gerpscore$pattern))

gerpscore=gerpscore[gerpscore$ncnumber %in% higher_than_5$V1,]

GGNG_gerp=gerpscore[gerpscore$pattern=='GGNG'|gerpscore$pattern=='CNCC',c(1,6)]
GGNG_gerp=aggregate(GGNG_gerp$GerpScore,list(GGNG_gerp$ncnumber),mean)
colnames(GGNG_gerp)=c('ncnumber','score')
GGNG_gerp$type="GGNG"

GNGG_gerp=gerpscore[gerpscore$pattern=='GNGG'|gerpscore$pattern=='CCNC',c(1,6)]
GNGG_gerp=aggregate(GNGG_gerp$GerpScore,list(GNGG_gerp$ncnumber),mean)
colnames(GNGG_gerp)=c('ncnumber','score')
GNGG_gerp$type="GNGG"
head(GNGG_gerp)

G3_gerp=gerpscore[gerpscore$pattern=='GGG'|gerpscore$pattern=='CCC',c(1,6)]
G3_gerp=aggregate(G3_gerp$GerpScore,list(G3_gerp$ncnumber),mean)
colnames(G3_gerp)=c('ncnumber','score')
G3_gerp$type="G3"
head(G3_gerp)

gerpcompare=rbind(G3_gerp,GGNG_gerp,GNGG_gerp)
head(gerpcompare)

spread_gerpcompare=tidyr::spread(gerpcompare,key=type,value=score)

spread_gerpcompare=merge(spread_gerpcompare,class_info[,c(1,5,6)])
spread_gerpcompare$euorpro=ifelse(spread_gerpcompare$Host %in% c("plant","alga","animal","fungi","protist"),"eukaryote",
                    ifelse(spread_gerpcompare$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

head(spread_gerpcompare)
GGNG=na.omit(spread_gerpcompare[,c(1,2,3,5,6,7)])
GGNG=tidyr::gather(GGNG,key=pattern,value=Con,G3,GGNG)
head(GGNG)


GNGG=na.omit(spread_gerpcompare[,c(1,2,4,5,6,7)])
GNGG=tidyr::gather(GNGG,key=pattern,value=Con,G3,GNGG)
head(GNGG)


ggplot(GGNG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim=c(-20,10))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GGNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.005,textsize = 10,y_position = 9,
              test.args=list(paired=T),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  ylab('RS Score')+
  guides(fill=F)

ggsave("GGNG_gerpcompare.pdf",width = 4,height =6)

ggplot(GNGG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim=c(-20,10))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GNGG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.005,textsize = 10,y_position = 9,
              test.args=list(paired=T),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  ylab('RS Score')+
  guides(fill=F)

ggsave("GNGG_gerpcompare.pdf",width = 4,height = 6)

head(GGNG)
ggplot(GGNG[GGNG$euorpro!="unclear",],aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim=c(-20,10))+
  facet_grid(.~euorpro,scales = "free_y")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GGNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.005,textsize = 10,y_position = 9,
              test.args=list(paired=T),family = 'serif',size = 1)+
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
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  ylab('RS Score')

ggsave("GGNG_euorpro_gerpcompare.pdf",width = 6,height =6)

ggplot(GNGG[GNGG$euorpro!="unclear",],aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim=c(-20,10))+
  facet_grid(.~euorpro,scales = "free_y")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GNGG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.005,textsize = 10,y_position = 9,
              test.args=list(paired=T),family = 'serif',size = 1)+
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
        axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  ylab('RS Score')

ggsave("GNGG_euorpro_gerpcompare.pdf",width = 6,height =6)

#######################################################

ggplot(GGNG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~Host,ncol = 4)+
  coord_cartesian(ylim=c(-20,10))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GGNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T),family = 'serif',size = 1)+
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
  ylab('RS Score')

ggsave("G3-GGNG-gerp_byHost.pdf",width = 8,height = 8)

ggplot(GGNG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~GenomeType,ncol = 4)+
  coord_cartesian(ylim=c(-20,10))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GGNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T),family = 'serif',size = 1)+
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
  ylab('RS Score')

ggsave("G3-GGNG-gerp_byGTP.pdf",width = 8,height = 8)

ggplot(GNGG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~Host,ncol = 4)+
  coord_cartesian(ylim=c(-20,10))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GNGG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T),family = 'serif',size = 1)+
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
  ylab('RS Score')

ggsave("G3-GNGG-gerp_byHost.pdf",width = 8,height = 8)

ggplot(GNGG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~GenomeType,ncol = 4)+
  coord_cartesian(ylim=c(-20,10))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GNGG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T),family = 'serif',size = 1)+
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
  ylab('RS Score')

ggsave("G3-GNGG-gerp_byGTP.pdf",width = 8,height = 8)

#############################################compare loops and gruns conservation##################################################

head(gerpscore)
gruns_and_loops=gerpscore[gerpscore$pattern=='GGG'|gerpscore$pattern=='CCC',c(1,7,8)]
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

gruns_and_loops$euorpro=ifelse(gruns_and_loops$Host %in% c("plant","alga","animal","fungi","protist"),"eukaryote",
                               ifelse(gruns_and_loops$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim=c(-20,10))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('grunscore','loopscore')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.005,textsize = 10,
              y_position = 9,test.args=list(paired=T,alternative='greater'),family = 'serif',size = 1)+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size=1.2))+
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  ylab('RS Score for G3-PQS')+
  guides(fill=F)+
  scale_x_discrete(breaks=c('grunscore','loopscore'),labels=c('G-track','Loop'))

ggsave("G3_gtrack_loop_conservation_gerp.pdf",width = 4,height =6)

gruns_and_loops$pattern=ifelse(gruns_and_loops$pattern=='loopscore','Loop','G-track')

head(gruns_and_loops)

ggplot(gruns_and_loops[gruns_and_loops$euorpro!="unclear",],aes(x=pattern,y=score,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
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
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank())+
  ylab('RS Score for G3-PQS')

ggsave("G3_euorpro_gtrack_loop_conservation_gerp.pdf",width = 6,height =6)


ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
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
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())+
  ylab('RS score for G3-PQS')


ggsave("G3_gtrack_loop_conservation_byHost.pdf",width = 8,height=8)


ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
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
  theme(legend.text = element_text(face = "bold",size=26,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())+
  ylab('RS score for G3-PQS')


ggsave("G3_gtrack_loop_conservation_byGtp.pdf",width = 8,height=8)



