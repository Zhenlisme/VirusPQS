rm(list = ls())
gc()

library(ggplot2)
library(ggpubr)


setwd('E:/课题数据备份/RevisedSeq/G4Phastcons/')

class_info=read.csv("../classinfo_viroid_satellite.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'

higher_than_5=read.csv('alignedcount_higher_than_5.txt',header = F,stringsAsFactors = F)
head(higher_than_5)

seqphast=read.csv2('All_pattern_phast.csv',header = F,sep = ',',stringsAsFactors = F)
colnames(seqphast)=c('ncnumber','pattern','start','end','seq','phastscore','grunscore','loopscore')
seqphast$phastscore=as.numeric(seqphast$phastscore)
head(seqphast)

seqphast=seqphast[seqphast$ncnumber %in% higher_than_5$V1,]
length(unique(seqphast$ncnumber))

GGNG_phast=seqphast[seqphast$pattern=='GGNG'|seqphast$pattern=='CNCC',c(1,6)]
GGNG_phast=aggregate(GGNG_phast$phastscore,list(GGNG_phast$ncnumber),mean)
colnames(GGNG_phast)=c('ncnumber','score')
GGNG_phast$type="GGNG"

GNGG_phast=seqphast[seqphast$pattern=='GNGG'|seqphast$pattern=='CCNC',c(1,6)]
GNGG_phast=aggregate(GNGG_phast$phastscore,list(GNGG_phast$ncnumber),mean)
colnames(GNGG_phast)=c('ncnumber','score')
GNGG_phast$type="GNGG"
head(GNGG_phast)

G3_phast=seqphast[seqphast$pattern=='GGG'|seqphast$pattern=='CCC',c(1,6)]
G3_phast=aggregate(G3_phast$phastscore,list(G3_phast$ncnumber),mean)
colnames(G3_phast)=c('ncnumber','score')
G3_phast$type="G3"
head(G3_phast)

phastcompare=rbind(G3_phast,GGNG_phast,GNGG_phast)
head(phastcompare)

spread_phastcompare=tidyr::spread(phastcompare,key=type,value=score)

spread_phastcompare=merge(spread_phastcompare,class_info[,c(1,5,6)])
spread_phastcompare$euorpro=ifelse(spread_phastcompare$Host %in% c("plant","alga","animal","fungi","protist"),"eukaryote",
                                  ifelse(spread_phastcompare$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

head(spread_phastcompare)
GGNG=na.omit(spread_phastcompare[,c(1,2,3,5,6,7)])
GGNG=tidyr::gather(GGNG,key=pattern,value=Con,G3,GGNG)
head(GGNG)


GNGG=na.omit(spread_phastcompare[,c(1,2,4,5,6,7)])
GNGG=tidyr::gather(GNGG,key=pattern,value=Con,G3,GNGG)
head(GNGG)

ggplot(GGNG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  geom_violin(trim=T,color="white")+
  coord_cartesian(ylim = c(0,1.1))+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GGNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,y_position = 1.05,
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
  ylab('Phastcons Score')+
  guides(fill=F)

ggsave("GGNG_phastcompare.pdf",width = 4,height =6)

ggplot(GNGG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  geom_violin(trim=T,color="white")+
  coord_cartesian(ylim = c(0,1.1))+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GNGG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,y_position = 1.05,
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
  ylab('Phastcons Score')+
  guides(fill=F)

ggsave("GNGG_phastcompare.pdf",width = 4,height =6)


ggplot(GGNG[GGNG$euorpro!="unclear",],aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim=c(0,1.1))+
  facet_grid(.~euorpro,scales = "free_y")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GGNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,y_position = 1.02,
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
  ylab('Phastcons Score')

ggsave("GGNG_euorpro_phastcompare.pdf",width = 6,height =6)

ggplot(GNGG[GNGG$euorpro!="unclear",],aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim=c(0,1.1))+
  facet_grid(.~euorpro,scales = "free_y")+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GNGG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,y_position = 1.02,
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
  ylab('Phastcons Score')

ggsave("GNGG_euorpro_phastcompare.pdf",width = 6,height =6)

###########################################

head(GGNG)
levels(factor(GGNG$pattern))

ggplot(GGNG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~Host,ncol = 4)+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GGNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T,alternative='less'),family = 'serif',size = 1)+
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
  ylab('Phastcons Score')+
  scale_y_continuous(breaks = c(0,0.5,1))

ggsave("G3-GGNG-phastconscompare_byHost.pdf",width = 8,height = 8)

ggplot(GGNG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~GenomeType,ncol = 3)+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GGNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T,alternative='less'),family = 'serif',size = 1)+
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
  ylab('Phastcons Score')+
  scale_y_continuous(breaks = c(0,0.5,1))

ggsave("G3-GGNG-phastconscompare_byGTP.pdf",width = 8,height = 8)

ggplot(GNGG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~Host,ncol = 4)+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GNGG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T,alternative='less'),family = 'serif',size = 1)+
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
  ylab('Phastcons Score')+
  scale_y_continuous(breaks = c(0,0.5,1))

ggsave("G3-GNGG-phastconscompare_byHost.pdf",width = 8,height = 8)

ggplot(GNGG,aes(x=pattern,y=Con,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~GenomeType,ncol = 4)+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G3','GNGG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T,alternative='less'),family = 'serif',size = 1)+
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
  ylab('Phastcons Score')+
  scale_y_continuous(breaks = c(0,0.5,1))

ggsave("G3-GNGG-phastconscompare_byGTP.pdf",width = 8,height = 8)

#############################################compare loops and gruns conservation##################################################
head(seqphast)
gruns_and_loops=seqphast[seqphast$pattern=='GGG'|seqphast$pattern=='CCC',c(1,7,8)]
gruns_and_loops$grunscore=as.numeric(gruns_and_loops$grunscore)
gruns_and_loops$loopscore=as.numeric(gruns_and_loops$loopscore)
head(gruns_and_loops)

gruns_and_loops=aggregate(gruns_and_loops[,c(2,3)],list(gruns_and_loops$ncnumber),mean)
colnames(gruns_and_loops)[1]=c('ncnumber')
head(gruns_and_loops)

gruns_and_loops=tidyr::gather(gruns_and_loops,key='pattern',value='score',grunscore,loopscore)

gruns_and_loops=merge(gruns_and_loops,class_info[,c(1,5,6)])

head(gruns_and_loops)
gruns_and_loops$euorpro=ifelse(gruns_and_loops$Host %in% c("plant","alga","animal","fungi","protist"),"eukaryote",
       ifelse(gruns_and_loops$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('grunscore','loopscore')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T,alternative='greater'),family = 'serif',size = 1)+
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
  ylab('Phastcons Score for G3-PQS')+
  guides(fill=F)+
  scale_x_discrete(breaks=c('grunscore','loopscore'),labels=c('G-track','Loop'))

ggsave("G3_gtrack_loop_conservation.pdf",width = 4,height =6)

head(gruns_and_loops)
gruns_and_loops$pattern=ifelse(gruns_and_loops$pattern=='loopscore','Loop','G-track')

ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~GenomeType,ncol = 4)+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G-track','Loop')),test = 'wilcox.test',
              map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T,alternative='greater'),family = 'serif',size = 1)+
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
  ylab('Phastcons Score for G3-PQS')+
  scale_y_continuous(breaks = c(0,0.5,1))
ggsave("G3_gtrack_loop_conservation_byGTP.pdf",width = 8,height=8)

ggplot(gruns_and_loops[gruns_and_loops$euorpro!="unclear",],aes(x=pattern,y=score,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~euorpro,ncol = 4)+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G-track','Loop')),test = 'wilcox.test',
              map_signif_level = T,tip_length = 0.02,textsize = 10,
              test.args=list(paired=T,exact=F,alternative='greater'),family = 'serif',size = 1)+
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
  ylab('Phastcons Score for G3-PQS')

ggsave("G3_euorpro_gtrack_loop_conservation.pdf",width = 6,height=6)

ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  facet_wrap(.~Host,ncol = 4)+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G-track','Loop')),test = 'wilcox.test',
              map_signif_level = T,tip_length = 0.02,textsize = 10,
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
  ylab('Phastcons Score for G3-PQS')+
  scale_y_continuous(breaks = c(0,0.5,1))


ggsave("G3_gtrack_loop_conservation_byHost.pdf",width = 8,height=8)

