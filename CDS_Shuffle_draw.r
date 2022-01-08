rm(list=ls())
gc()


library(ggplot2)
#library(ggpubr)
#library(ggsci)
setwd('D:/毕业文件/李振数据备份/课题数据备份/RevisedSeq/CDS shuffling/')

class_info=read.csv("classinfo_viroid_satellite.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'
class_info$GenomeType=factor(class_info$GenomeType,levels = c('dsDNA','dsDNA-RT','dsRNA','ssDNA','ssRNA-RT','ssRNA(-)','ssRNA(+)','satellite','viroid','unclear'),ordered = T)
head(class_info)

G2_compare=read.csv2("G24_pvalues.txt",sep=',',header = T,stringsAsFactors = F)
G2_compare$zscore=as.numeric(G2_compare$zscore)
comparision_G2=tidyr::spread(G2_compare[,c(1,2,6)],key=pattern,value="zscore")
wilcox.test(comparision_G2$G2,comparision_G2$GNG,alternative = 'greater',paired = T)

ggplot(G2_compare,aes(zscore,color=pattern))+
  geom_density()+
  #stat_ecdf(size=0.5)+
  coord_cartesian(xlim = c(-4,8))+
  ylab("Cumulative Distribution")+
  xlab("Z-score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  scale_color_manual(values = c("#E7B800","#00AFBB"))+
  theme(legend.position = "top",
        legend.title = element_blank())+
  theme(panel.background = element_blank())+
  scale_x_continuous(limits = c(-4,8),breaks=seq(-4,8,2))

ggsave("G2PQS_zscore.pdf",width = 6,height =6)
head(class_info)


colnames(G2_compare)[1]='proid'
G2_compare$ncnumber=as.vector(sapply(G2_compare$proid,function(x){strsplit(x,split = "_cds_")[[1]][1]}))

G2_compare=merge(G2_compare,class_info[c('ncnumber','Host')])
G2_compare$europro=ifelse(G2_compare$Host %in% c('fungi',"plant","alga","animal","protist"),"eukaryote",
                          ifelse(G2_compare$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

head(G2_compare)

eukaryote_G2=G2_compare[G2_compare$europro=='eukaryote',]
prokaryote_G2=G2_compare[G2_compare$europro=='prokaryote',]
head(G2_compare)

spread_eukaryote=tidyr::spread(eukaryote_G2[c(2,3,7)],key=pattern,value=zscore)
spread_prokaryote=tidyr::spread(prokaryote_G2[c(2,3,7)],key=pattern,value=zscore)
spread_eukaryote=na.omit(spread_eukaryote)
spread_prokaryote=na.omit(spread_prokaryote)
spread_eukaryote=spread_eukaryote[is.finite(spread_eukaryote$G2) & is.finite(spread_eukaryote$GNG),]
spread_prokaryote=spread_prokaryote[is.finite(spread_prokaryote$G2) & is.finite(spread_prokaryote$GNG),]

spread_prokaryote[is.infinite(spread_prokaryote),]

options(digits = 10,scipen = 5)
wilcox.test(spread_eukaryote$G2,spread_eukaryote$GNG,alternative = 'greater',paired = T)$p.value  ##3.62535714e-13
wilcox.test(spread_prokaryote$G2,spread_prokaryote$GNG,alternative = 'greater',paired = T)$p.value ##0

spread_eukaryote=tidyr::gather(spread_eukaryote,key=pattern,value=zscore,c(2,3))

head(spread_eukaryote)
G2_spread_eukaryotic_plot=ggplot(spread_eukaryote,aes(zscore,color=pattern))+
  geom_density()+
  ylab("Density")+
  xlab("Z-score (eukaryotic virus)")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  scale_color_manual(values = c("#E7B800","#00AFBB"))+
  theme(legend.position = "top",
        legend.title = element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.position = c(0.85, 0.7),
        legend.title = element_blank())+
  coord_cartesian(xlim = c(-4,8))
  #scale_x_continuous(limits = c(-4,8),breaks=seq(-4,8,2))
print(G2_spread_eukaryotic_plot)

spread_prokaryote=tidyr::gather(spread_prokaryote,key=pattern,value=zscore,c(2,3))

head(spread_prokaryote)

G2_spread_prokaryote_plot=ggplot(spread_prokaryote,aes(zscore,color=pattern))+
  geom_density()+
  ylab("Density")+
  xlab("Z-score (prokaryotic virus)")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  scale_color_manual(values = c("#E7B800","#00AFBB"))+
  theme(legend.position = "top",
        legend.title = element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.position = c(0.85, 0.7),
        legend.title = element_blank())+
  coord_cartesian(xlim = c(-4,8))

print(G2_spread_prokaryote_plot)


cowplot::plot_grid(G2_spread_eukaryotic_plot,G2_spread_prokaryote_plot,ncol = 1)
ggsave("G2PQS_zscore_Host.pdf",width = 6,height =6)



#####################G34##############################

G3_compare=read.csv2("G34_pvalues.txt",sep=',',header = T,stringsAsFactors = F)
G3_compare$zscore=as.numeric(G3_compare$zscore)
head(G3_compare)
comparision_G3=tidyr::spread(G3_compare[,c(1,2,6)],key=pattern,value="zscore")
wilcox.test(comparision_G3$G3,comparision_G3$GGNG,alternative = 'greater',paired = T)
wilcox.test(comparision_G3$G3,comparision_G3$GNGG,alternative = 'greater',paired = T)


ggplot(G3_compare,aes(zscore,color=pattern))+
  stat_ecdf(size=0.5)+
  coord_cartesian(xlim = c(-4,6))+
  ylab("Cumulative Distribution")+
  xlab("Z-score")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  scale_color_manual(values = c("#E7B800","#D54C4C","#50CB93"))+
  theme(legend.position = "top",
        legend.title = element_blank())+
  theme(panel.background = element_blank())+
  scale_x_continuous(limits = c(-4,6),breaks=seq(-4,6,2))

ggsave("G3PQS_zscore.pdf",width = 6,height =6)
head(class_info)

colnames(G3_compare)[1]='proid'
G3_compare$ncnumber=as.vector(sapply(G3_compare$proid,function(x){strsplit(x,split = "_cds_")[[1]][1]}))
G3_compare=merge(G3_compare,class_info[c('ncnumber','Host')])
G3_compare$europro=ifelse(G3_compare$Host %in% c('fungi',"plant","alga","animal","protist"),"eukaryote",
                          ifelse(G3_compare$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))
head(G3_compare)
G3_compare_m=tidyr::spread(G3_compare[,c(1,2,3,7)],key=pattern,value=zscore)
G2_compare_m=tidyr::spread(G2_compare[,c(1,2,3,7)],key=pattern,value=zscore)
head(G2_compare_m)

PQS=merge(G3_compare_m,G2_compare_m[,c(2,3,4)],by='proid')

PQS=merge(PQS,unique(class_info[,c(1,3,4,5,6,8,9)]))
head(PQS)
write.table(PQS,file = "PQS_zscore.csv",sep = ',',quote = F,col.names = T,row.names = F)

head(class_info)

G3_euka_compare=G3_compare[G3_compare$europro=='eukaryote',c(2,3,7)]
G3_pro_compare=G3_compare[G3_compare$europro=='prokaryote',c(2,3,7)]


spread_eukaryote=tidyr::spread(G3_euka_compare,key=pattern,value=zscore)
spread_prokaryote=tidyr::spread(G3_pro_compare,key=pattern,value=zscore)
spread_eukaryote=na.omit(spread_eukaryote)
spread_prokaryote=na.omit(spread_prokaryote)
spread_eukaryote=spread_eukaryote[is.finite(spread_eukaryote$G3) & is.finite(spread_eukaryote$GGNG) & is.finite(spread_eukaryote$GNGG),]
spread_prokaryote=spread_prokaryote[is.finite(spread_prokaryote$G3) & is.finite(spread_prokaryote$GGNG) & is.finite(spread_prokaryote$GNGG),]

wilcox.test(spread_eukaryote$G3,spread_eukaryote$GNGG,alternative = 'greater',paired=T)$p.value ##0
wilcox.test(spread_eukaryote$G3,spread_eukaryote$GGNG,alternative = 'greater',paired = T)$p.value ##0

wilcox.test(spread_prokaryote$G3,spread_prokaryote$GGNG,alternative = 'greater',paired = T)$p.value  ##0
wilcox.test(spread_prokaryote$G3,spread_prokaryote$GNGG,alternative = 'greater',paired = T)$p.value  ##0

spread_eukaryote=tidyr::gather(spread_eukaryote,key=pattern,value=zscore,c(2,3,4))

spread_eukaryote_G3=ggplot(spread_eukaryote,aes(zscore,color=pattern))+
  geom_density()+
  #facet_grid(.~europro)+
  #stat_ecdf(size=0.5)+
  coord_cartesian(xlim = c(-2,2))+
  ylab("Density")+
  xlab("Z-score (eukaryotic virus)")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  scale_color_manual(values = c("#E7B800","#D54C4C","#50CB93"))+
  theme(legend.position = "top",
        legend.title = element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.position = c(0.85, 0.7),
        legend.title = element_blank())


spread_prokaryote=tidyr::gather(spread_prokaryote,key=pattern,value=zscore,c(2,3,4))
spread_prokaryote_G3=ggplot(spread_prokaryote,aes(zscore,color=pattern))+
  geom_density()+
  #facet_grid(.~europro)+
  #stat_ecdf(size=0.5)+
  coord_cartesian(xlim = c(-2,2))+
  ylab("Density")+
  xlab("Z-score (prokaryotic virus)")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  scale_color_manual(values = c("#E7B800","#D54C4C","#50CB93"))+
  theme(legend.position = "top",
        legend.title = element_blank())+
  theme(panel.background = element_blank())+
  theme(legend.position = c(0.85, 0.7),
        legend.title = element_blank())

cowplot::plot_grid(spread_eukaryote_G3,spread_prokaryote_G3,ncol = 1)

ggsave("G3PQS_zscore_Host.pdf",width = 6,height =6)
