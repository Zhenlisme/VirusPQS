rm(list = ls())
gc()

library(ggplot2)
library(ggpubr)


setwd('E:/PQS_in_Plantvirus/All_viruses/FullSegment/RevisedSeq/G4Phastcons/')

class_info=read.csv2("../class_add_taxid.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$GenomeType=="circle-ssRNA","GenomeType"]="unknown"
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
length(higher_than_5$V1)

GNG_phast=seqphast[seqphast$pattern=='GNG'|seqphast$pattern=='CNC',c(1,6)]
GNG_phast=aggregate(GNG_phast$phastscore,list(GNG_phast$ncnumber),mean)
colnames(GNG_phast)=c('ncnumber','score')
head(GNG_phast)

G2_phast=seqphast[seqphast$pattern=='GG'|seqphast$pattern=='CC',c(1,6)]
G2_phast=aggregate(G2_phast$phastscore,list(G2_phast$ncnumber),mean)
colnames(G2_phast)=c('ncnumber','score')
head(G2_phast)
phastcompare=merge(G2_phast,GNG_phast,by='ncnumber')
colnames(phastcompare)=c('ncnumber','G2','GNG')

head(phastcompare)
gathered_phastcompare=tidyr::gather(phastcompare,key='pattern',value=score,G2,GNG)

head(gathered_phastcompare)
gathered_phastcompare=merge(gathered_phastcompare,class_info[,c(1,5,6)])

gathered_phastcompare$euorpro=ifelse(gathered_phastcompare$Host %in% c("plant","animal","fungi","protist"),"eukaryote",
                                     ifelse(gathered_phastcompare$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

ggplot(gathered_phastcompare,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G2','GNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
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
  ylab('Phastcons Score')+
  guides(fill=F)

ggsave("phastconscompare.pdf",width = 4,height =6)

ggplot(gathered_phastcompare[gathered_phastcompare$euorpro!="unclear",],aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  facet_wrap(.~euorpro,ncol = 4)+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G2','GNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
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
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  ylab('Phastcons Score')

ggsave("euorpro_phastconscompare.pdf",width = 6,height =6)


ggplot(gathered_phastcompare,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  facet_wrap(.~GenomeType,ncol = 4)+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G2','GNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
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
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  ylab('Phastcons Score')

ggsave("phastconscompare_byGtp.pdf",width = 8,height = 8)

ggplot(gathered_phastcompare,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  facet_wrap(.~Host,ncol = 4)+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('G2','GNG')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
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
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank())+
  ylab('Phastcons Score')

ggsave("phastconscompare_byHost.pdf",width = 8,height = 8)

#############################################compare loops and gruns conservation##################################################
head(seqphast)
gruns_and_loops=seqphast[seqphast$pattern=='GG'|seqphast$pattern=='CC',c(1,7,8)]
gruns_and_loops$grunscore=as.numeric(gruns_and_loops$grunscore)
gruns_and_loops$loopscore=as.numeric(gruns_and_loops$loopscore)
head(gruns_and_loops)

gruns_and_loops=aggregate(gruns_and_loops[,c(2,3)],list(gruns_and_loops$ncnumber),mean)
colnames(gruns_and_loops)[1]=c('ncnumber')
head(gruns_and_loops)

gruns_and_loops=tidyr::gather(gruns_and_loops,key='pattern',value='score',grunscore,loopscore)

gruns_and_loops=merge(gruns_and_loops,class_info[,c(1,5,6)])

head(gruns_and_loops)
gruns_and_loops$euorpro=ifelse(gruns_and_loops$Host %in% c("plant","animal","fungi","protist"),"eukaryote",
       ifelse(gruns_and_loops$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
  coord_cartesian(ylim=c(0,1.19))+
  geom_violin(trim=T,color="white")+
  geom_boxplot(width=0.2,outlier.alpha = 0)+
  geom_signif(comparisons = list(c('grunscore','loopscore')),test = 'wilcox.test',map_signif_level = T,tip_length = 0.02,textsize = 10,
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
  ylab('Phastcons Score for G2-PQS')+
  guides(fill=F)+
  scale_x_discrete(breaks=c('grunscore','loopscore'),labels=c('G-track','Loop'))

ggsave("gtrack_loop_conservation.pdf",width = 4,height =6)

head(gruns_and_loops)
gruns_and_loops$pattern=ifelse(gruns_and_loops$pattern=='loopscore','Loop','G-track')

ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
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
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())+
  ylab('Phastcons Score for G2-PQS')
ggsave("gtrack_loop_conservation_byGTP.pdf",width = 8,height=8)

ggplot(gruns_and_loops[gruns_and_loops$euorpro!="unclear",],aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
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
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())+
  ylab('Phastcons Score for G2-PQS')

ggsave("euorpro_gtrack_loop_conservation.pdf",width = 6,height=6)

ggplot(gruns_and_loops,aes(x=pattern,y=score,fill=pattern))+
  scale_fill_npg()+
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
  theme(legend.text = element_text(face = "bold",size=24,family="serif"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "top")+
  theme(axis.text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank())+
  ylab('Phastcons Score for G2-PQS')


ggsave("gtrack_loop_conservation_byHost.pdf",width = 8,height=8)

