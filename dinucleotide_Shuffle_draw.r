rm(list=ls())
gc()


library(ggplot2)
library(cowplot)


setwd('D:/毕业文件/李振数据备份/课题数据备份/RevisedSeq/dinucleotide_shuffle/')

class_info=read.csv("../classinfo_viroid_satellite.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'
class_info$GenomeType=factor(class_info$GenomeType,levels = c('dsDNA','dsDNA-RT','dsRNA','ssDNA','ssRNA-RT','ssRNA(-)','ssRNA(+)','satellite','viroid','unclear'),ordered = T)
head(class_info)

G2_compare=read.csv2("G2_compare.csv",sep=',',header = T,stringsAsFactors = F)
G2_compare$zscore=as.numeric(G2_compare$zscore)

G3_compare=read.csv2("G3_compare.csv",sep=',',header = T,stringsAsFactors = F)
G3_compare$zscore=as.numeric(G3_compare$zscore)


head(G2_compare)

comparision_G2=tidyr::spread(G2_compare[,c(1,2,6,9)],key=pattern,value="zscore")

euro_comparision_G2=comparision_G2[comparision_G2$europro=='eukaryote',]
pro_comparision_G2=comparision_G2[comparision_G2$europro=='prokaryote',]
wilcox.test(euro_comparision_G2$G2,euro_comparision_G2$GNG,paired = T,alternative = 'greater')$p.value  ##1.6486246170996337e-56
wilcox.test(pro_comparision_G2$G2,pro_comparision_G2$GNG,paired = T,alternative = 'greater')$p.value  ## 0

head(G3_compare)
G3_compare$europro=ifelse(G3_compare$Host %in% c('fungi',"plant","alga","animal","protist"),"eukaryote",
                          ifelse(G3_compare$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))
comparision_G3=tidyr::spread(G3_compare[,c(1,2,6,9)],key=pattern,value="zscore")
head(comparision_G3)

wilcox.test(comparision_G3$G3,comparision_G3$GGNG,alternative = 'greater',paired = T)
wilcox.test(comparision_G3$G3,comparision_G3$GNGG,alternative = 'greater',paired = T)

euro_comparision_G3=comparision_G3[comparision_G3$europro=='eukaryote',]
pro_comparision_G3=comparision_G3[comparision_G3$europro=='prokaryote',]

wilcox.test(euro_comparision_G3$G3,euro_comparision_G3$GGNG,paired = T,alternative = 'greater')$p.value #0.0000058437914770202239

wilcox.test(euro_comparision_G3$G3,euro_comparision_G3$GNGG,paired = T,alternative = 'greater')$p.value #0.0000000025000861496898385

wilcox.test(pro_comparision_G3$G3,pro_comparision_G3$GGNG,paired = T,alternative = 'less')$p.value #1.437650551270429e-115
wilcox.test(pro_comparision_G3$G3,pro_comparision_G3$GNGG,paired = T,alternative = 'less')$p.value #3.8661781073817959e-30


ggplot(test,aes(zscore,color=pattern))+
  geom_density()+
  coord_cartesian(xlim = c(-10,20))
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
  theme(panel.background = element_blank())


ggplot(G2_compare,aes(zscore,color=pattern))+
  #stat_ecdf(size=0.5)+
  geom_density()+
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
  theme(panel.background = element_blank())

ggsave("G2PQS_zscore.pdf",width = 6,height =6)

head(pro_comparision_G2)
pro_comparision_G2=tidyr::gather(pro_comparision_G2,key=pattern,value=zscore,c(3,4))
euro_comparision_G2=tidyr::gather(euro_comparision_G2,key=pattern,value=zscore,c(3,4))



euro_G2_plot=ggplot(euro_comparision_G2,aes(zscore,color=pattern))+
  #facet_grid(.~europro,scales = 'free_y')+
  coord_cartesian(xlim = c(-4,15))+
  #stat_ecdf(size=0.5)+
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
  theme(legend.position = c(0.85, 0.7),
        legend.title = element_blank())+
  theme(panel.background = element_blank())

pro_G2_plot=ggplot(pro_comparision_G2,aes(zscore,color=pattern))+
  #facet_grid(.~europro,scales = 'free_y')+
  #coord_cartesian(xlim = c(-4,10))+
  #stat_ecdf(size=0.5)+
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
  theme(legend.position = c(0.85, 0.7),
        legend.title = element_blank())+
  theme(panel.background = element_blank())

plot_grid(euro_G2_plot,pro_G2_plot,ncol = 1)
ggsave("G2PQS_zscore_Host.pdf",width = 6,height =6)

################################################################################################
head(G3_compare)
comparision_G3=tidyr::spread(G3_compare[,c(1,2,6)],key=pattern,value="zscore")
wilcox.test(comparision_G3$G3,comparision_G3$GGNG,alternative = 'greater',paired = T)
wilcox.test(comparision_G3$G3,comparision_G3$GNGG,alternative = 'greater',paired = T)

ggplot(G3_compare,aes(zscore,color=pattern))+
  #stat_ecdf(size=0.5)+
  geom_density()+
  coord_cartesian(xlim = c(-10,40))+
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
  theme(panel.background = element_blank())

ggsave("G3PQS_zscore.pdf",width = 6,height =6)



head(pro_comparision_G3)
pro_comparision_G3=tidyr::gather(pro_comparision_G3,key=pattern,value=zscore,c(3,4,5))
euro_comparision_G3=tidyr::gather(euro_comparision_G3,key=pattern,value=zscore,c(3,4,5))

euro_G3_plot=ggplot(euro_comparision_G3,aes(zscore,color=pattern))+
  #facet_grid(.~europro,scales = 'free_y')+
  coord_cartesian(xlim = c(-4,10))+
  #stat_ecdf(size=0.5)+
  geom_density()+
  ylab("Density")+
  xlab("Z-score (eukaryotic virus)")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  scale_color_manual(values = c("#E7B800","#D54C4C","#50CB93"))+
  theme(legend.position = c(0.8, 0.7),
        legend.title = element_blank())+
  theme(panel.background = element_blank())

pro_G3_plot=ggplot(pro_comparision_G3,aes(zscore,color=pattern))+
  #facet_grid(.~europro,scales = 'free_y')+
  coord_cartesian(xlim = c(-4,10))+
  #stat_ecdf(size=0.5)+
  geom_density()+
  ylab("Density")+
  xlab("Z-score (prokaryotic virus)")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.text = element_text(size=26,face = "bold",colour = "black",family="serif"),
        text = element_text(size=26,face = "bold",colour = "black",family="serif"))+
  scale_color_manual(values = c("#E7B800","#D54C4C","#50CB93"))+
  theme(legend.position = c(0.8, 0.7),
        legend.title = element_blank())+
  theme(panel.background = element_blank())

plot_grid(euro_G3_plot,pro_G3_plot,ncol = 1)
ggsave("G3PQS_zscore_Host.pdf",width = 6,height =6)


plot_grid(euro_G2_plot,euro_G3_plot,pro_G2_plot,pro_G3_plot,ncol = 2)

ggsave("PQS_zscore_europro_disshuffle.pdf",width = 8,height =8)

