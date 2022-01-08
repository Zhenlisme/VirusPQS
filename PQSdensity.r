rm(list = ls())
gc()

library(ggpubr)
library(RColorBrewer)
library(ggplot2)
library(ggsci)

setwd("D:/毕业文件/李振数据备份/课题数据备份/RevisedSeq/BasicStatistics/")

class_info=read.csv("../classinfo_viroid_satellite.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'

head(class_info)
class_info$GenomeType=factor(class_info$GenomeType,
      levels = c('dsDNA','dsDNA-RT','dsRNA','ssDNA','ssRNA-RT','ssRNA(-)','ssRNA(+)','satellite','viroid','unclear'),ordered = T)

unique_taxid=unique(class_info[,c('GenomeType','Host','Taxid')])

aggregate(unique_taxid$Taxid,list(unique_taxid$GenomeType,unique_taxid$Host),length)
       

ggplot(unique(class_info[,c('GenomeType','Host','Taxid')]),aes(x=GenomeType,fill=Host))+
  geom_bar()+
  scale_fill_d3()+
  ylab("Species Count")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  theme(axis.title.x = element_blank(),legend.position = "top",
        axis.text.y = element_text(size=24,colour = "black",family="serif"))+
  theme(axis.text.x = element_text(size=24,colour = "black",family="serif",angle = -60,hjust = -.05),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))

ggsave("virus_species.pdf",width = 8,height =8)


sequence_length=read.csv2("sequence_length.csv",sep=",",header = T,stringsAsFactors = F)
sequence_length=unique(sequence_length[,c(1,4)])
sequence_length$genomelength=2*sequence_length$genomelength
head(sequence_length)


Allpattern_df=read.csv2("long_pattern_count.csv",sep=",",header = T,stringsAsFactors = F)
head(Allpattern_df)

Allpattern_df$G2=Allpattern_df$GG_count+Allpattern_df$CC_count
#Allpattern_df$GNG=Allpattern_df$GNG_count+Allpattern_df$CNC_count
Allpattern_df$G3=Allpattern_df$GGG_count+Allpattern_df$CCC_count
#Allpattern_df$GNGG=Allpattern_df$GNGG_count+Allpattern_df$CCNC_count
#Allpattern_df$GGNG=Allpattern_df$GGNG_count+Allpattern_df$CNCC_count

Allpattern_df=Allpattern_df[,c("ncnumber","G2","G3")]

Allpattern_df=merge(Allpattern_df,class_info[,c("ncnumber","Taxid")],all.y = T)
head(Allpattern_df)
Allpattern_df=merge(Allpattern_df,sequence_length,all.x = T)

head(Allpattern_df)

merged_taxid_df=aggregate(Allpattern_df[,c(2,3,5)],list(Allpattern_df$Taxid),sum)
colnames(merged_taxid_df)[1]='Taxid'
head(merged_taxid_df)

merged_taxid_df$G2=1000*merged_taxid_df$G2/merged_taxid_df$genomelength
merged_taxid_df$G3=1000*merged_taxid_df$G3/merged_taxid_df$genomelength

merged_taxid_df[merged_taxid_df$Taxid=='1000373',]

head(class_info)

head(unique(class_info[,c('GenomeType','Host','Taxid')]))

merged_taxid_df=merge(unique(class_info[,c('GenomeType','Host','Taxid')]),merged_taxid_df,by='Taxid',all.y = T)
head(merged_taxid_df)
PQSdensity=tidyr::gather(merged_taxid_df,key=pattern,value=density,G2,G3)
head(PQSdensity)

tapply(PQSdensity$density, list(PQSdensity$pattern,PQSdensity$GenomeType), mean)
tapply(PQSdensity$density, list(PQSdensity$pattern,PQSdensity$Host), mean)
tapply(PQSdensity$density, list(PQSdensity$pattern), mean)

head(PQSdensity)
unique(PQSdensity$Host)

PQSdensity$europro=ifelse(PQSdensity$Host %in% c('fungi',"plant","animal","protist","alga"),"eukaryote",
                          ifelse(PQSdensity$Host %in% c("bacteria","archaea"),"prokaryote","unclear"))

tapply(PQSdensity$density, list(PQSdensity$pattern,PQSdensity$europro), mean)
tapply(PQSdensity$density, list(PQSdensity$pattern,PQSdensity$Host), mean)

write.table(t(tapply(PQSdensity$density, list(PQSdensity$pattern,PQSdensity$GenomeType), mean)),
            "genometype_density.csv",col.names = T,row.names = T,quote = F,sep = ",")

ggplot(PQSdensity[PQSdensity$europro!="unclear",],aes(x=europro,y=density,fill=europro))+
  scale_fill_manual(values = c('#FF585D',"#FFB549"))+
  geom_violin()+
  ylab("Frequency (PQSs/1000 nt)")+
  facet_grid(pattern~.,scales = 'free_y')+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  guides(fill=F)+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=24,colour = "black",family="serif"))+
  theme(axis.text.x = element_text(size=24,colour = "black",family="serif",angle = -60,hjust = -.05),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  stat_summary(fun = mean, geom = "point", 
               aes(x=europro,y=density),shape=16,size=2,color="blue",position = position_dodge(0.8))
  
ggsave("Eurorpro_density.pdf",width = 6,height =6)


ggplot(PQSdensity,aes(x=GenomeType,y=density,fill=pattern))+
  scale_fill_npg()+
  geom_violin()+
  facet_grid(pattern~.,scales = 'free')+
  ylab("Frequency (PQSs/1000 nt)")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  guides(fill=F)+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=24,colour = "black",family="serif"))+
  theme(axis.text.x = element_text(size=24,colour = "black",family="serif",angle = -60,hjust = -.05),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  stat_summary(fun = mean, geom = "point", 
               aes(x=GenomeType,y=density),shape=16,size=2,color="blue",position = position_dodge(0.8))

ggsave("Genometype_density.pdf",width = 6,height =6)

head(PQSdensity)

PQSdensity$Host=factor(PQSdensity$Host,levels = c("animal","plant","alga","fungi","protist","archaea","bacteria","unclear"),ordered = T)
levels(PQSdensity$Host)


ggplot(PQSdensity,aes(x=Host,y=density,fill=pattern))+
  scale_fill_npg()+
  geom_violin()+
  facet_grid(pattern~.,scales = 'free')+
  ylab("Frequency (PQSs/1000 nt)")+
  theme_set(theme_bw())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black",size = 1.2))+
  guides(fill=F)+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=24,colour = "black",family="serif"))+
  theme(axis.text.x = element_text(size=24,colour = "black",family="serif",angle = -60,hjust = -.05),
        text = element_text(size=24,face = "bold",colour = "black",family="serif"))+
  stat_summary(fun = mean, geom = "point", 
               aes(x=Host,y=density),shape=16,size=2,color="blue",position = position_dodge(0.8))

ggsave("HostType_density.pdf",width = 6,height =6)



