
rm(list = ls())
gc()

library(ggpubr)
library(ggplot2)

setwd("E:/课题数据备份/RevisedSeq/GenomeFeature//")

class_info=read.csv("../classinfo_viroid_satellite.csv",sep = ",",header = F,stringsAsFactors = F)
colnames(class_info)=c("ncnumber","TrueSpeciesName","Family","Genus","GenomeType","Host","Segment","Species","Taxid")

class_info[class_info$Segment=='non','Taxid']=paste(class_info[class_info$Segment=='non','Taxid'],
                                                    class_info[class_info$Segment=='non','ncnumber'],sep = '-')
class_info[class_info$GenomeType=='unknown','GenomeType']='unclear'

head(class_info)
class_info$GenomeType=factor(class_info$GenomeType,
                             levels = c('dsDNA','dsDNA-RT','dsRNA','ssDNA','ssRNA-RT','ssRNA(-)','ssRNA(+)','satellite','viroid','unclear'),ordered = T)

genomfeature_PQS=read.csv2("genomefeature_PQS.txt",sep = "\t",stringsAsFactors = F,header = F)
colnames(genomfeature_PQS)=c("ncnumber","genestart","genend","genomefeature","gene_strand","PQSstart","PQSend","PQStype","PQSstrand")
genomfeature_PQS$genelength=genomfeature_PQS$genend-genomfeature_PQS$genestart+1

G2_non_intersect=read.csv2("non_intersect.bed6",sep = "\t",stringsAsFactors = F,header = F)
colnames(G2_non_intersect)=c("ncnumber","genestart","genend","genomefeature","value","gene_strand")

G2_non_intersect$PQSstart=0
G2_non_intersect$PQSend=0
G2_non_intersect$PQStype="G2"
G2_non_intersect$PQSstrand='.'
G2_non_intersect$genelength=G2_non_intersect$genend-G2_non_intersect$genestart+1

G3_non_intersect=G2_non_intersect
G3_non_intersect$PQStype="G3"

non_intersect=rbind(G2_non_intersect,G3_non_intersect)[,c(1,2,3,4,6,7,8,9,10,11)]

genomfeature_PQS=rbind(genomfeature_PQS,non_intersect)
unclear_features=c("gap","sequence_alteration","modified_DNA_base","operon","polyA_signal_sequence","sequence_feature","transcript","polyA_site","sequence_conflict",
                  "sequence_difference","sequence_uncertainty","DIRECT","INVERTED")
genomfeature_PQS=genomfeature_PQS[! genomfeature_PQS$genomefeature %in% unclear_features,]

genomfeature_PQS[genomfeature_PQS$genomefeature=='LONG_TERMINAL_REPEAT','genomefeature']='long_terminal_repeat'
head(genomfeature_PQS)

drow_function=function(genomfeature_PQS,PQS,GTP){
  if(GTP!=""){
  genomfeature_PQS=genomfeature_PQS[genomfeature_PQS$ncnumber %in% class_info[class_info$GenomeType==GTP,"ncnumber"],]
  }
  
  genomfeature_PQS$PQScount=ifelse(genomfeature_PQS$gene_strand==genomfeature_PQS$PQSstrand,1,0)
  #head(genomfeature_PQS)
  genomfeature_PQS=genomfeature_PQS[,c(1,2,3,4,6,7,8,10,11)]
  genomfeature_PQS$genemarker=paste(genomfeature_PQS$ncnumber,genomfeature_PQS$genestart,
                            genomfeature_PQS$genend,genomfeature_PQS$genomefeature,genomfeature_PQS$genelength,sep = "-")
  #head(genomfeature_PQS)
  genomfeature_PQS_summary=aggregate(genomfeature_PQS$PQScount,list(genomfeature_PQS$genemarker,genomfeature_PQS$PQStype),sum)
  #head(genomfeature_PQS_summary)
  colnames(genomfeature_PQS_summary)=c("genemarker","PQStype","Freq")
  gene_and_length=sapply(as.character(genomfeature_PQS_summary$genemarker), function(x){strsplit(x,split = "-")[[1]][c(4:5)]})
  
  gene_and_length=as.data.frame(matrix(gene_and_length,ncol = 2,byrow = T))
  genomfeature_PQS_summary$Locus=gene_and_length$V1
  genomfeature_PQS_summary$length=gene_and_length$V2
  genomfeature_PQS_summary$density=1000*genomfeature_PQS_summary$Freq/as.numeric(genomfeature_PQS_summary$length)
  #head(genomfeature_PQS_summary)
  
  genomfeature_PQS_summary=aggregate(genomfeature_PQS_summary$density,
                                     list(genomfeature_PQS_summary$PQStype,genomfeature_PQS_summary$Locus),mean)
  colnames(genomfeature_PQS_summary)=c("PQStype","Locus","Density")
  genomfeature_PQS_summary$Density=round(genomfeature_PQS_summary$Density,3)
  head(genomfeature_PQS_summary)
  
  aggregate(genomfeature_PQS_summary$Density,list(genomfeature_PQS_summary$PQStype),mean)
  
  genomfeature_PQS_summary=genomfeature_PQS_summary[genomfeature_PQS_summary$Density>0,]

  if(PQS=="G2"){
    genomfeature_PQS_summary=genomfeature_PQS_summary[genomfeature_PQS_summary$PQStype=="G2",]
    genomfeature_PQS_summary$Locus=factor(genomfeature_PQS_summary$Locus,
                    levels = genomfeature_PQS_summary[order(genomfeature_PQS_summary$Density),"Locus"],ordered=T)
  pG2=ggplot(genomfeature_PQS_summary,aes(x=Locus,y=Density))+
    ylim(0,22)+
    geom_bar(stat = "identity")+
    coord_flip()+
    ylab("Frequency (G2-PQSs/1,000 nt)")+
    theme_set(theme_bw())+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour = "black",size = 1.2))+
    theme(axis.title.y = element_blank())+
    theme(axis.text = element_text(face = "bold",colour = "black",family="serif"),
          text = element_text(face = "bold",colour = "black",family="serif",size = 18),
          legend.position = 'top',
          legend.title = element_blank())+
    labs(title=GTP)+
    theme(plot.title = element_text(hjust = 0.5,size=18))+
    geom_text(aes(label=Density,vjust=0.5,hjust=0))
  return(pG2)}
  if(PQS=="G3"){
    genomfeature_PQS_summary=genomfeature_PQS_summary[genomfeature_PQS_summary$PQStype=="G3",]
    genomfeature_PQS_summary$Locus=factor(genomfeature_PQS_summary$Locus,
                     levels = genomfeature_PQS_summary[order(genomfeature_PQS_summary$Density),"Locus"],ordered=T)
  pG3=ggplot(genomfeature_PQS_summary,aes(x=Locus,y=Density))+
    ylim(0,1.2)+
    geom_bar(stat = "identity")+
    coord_flip()+
    ylab("Frequency (G3-PQSs/1,000 nt)")+
    theme_set(theme_bw())+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour = "black",size = 1.2))+
    theme(axis.title.y = element_blank())+
    theme(axis.text = element_text(face = "bold",colour = "black",family="serif"),
          text = element_text(face = "bold",colour = "black",family="serif",size=18),
          legend.position = 'top',
          legend.title = element_blank())+
    labs(title=GTP)+
    theme(plot.title = element_text(hjust = 0.5,size = 18))+
    geom_text(aes(label=Density,vjust=0.5,hjust=0))
  return(pG3)}
}


drow_function(genomfeature_PQS,"G2","")
ggsave("all_genomefeature_G2.pdf",width = 6,height = 10)

drow_function(genomfeature_PQS,"G3","")
ggsave("all_genomefeature_G3.pdf",width = 6,height = 10)


drow_function(genomfeature_PQS,"G2","dsDNA")
ggsave("G2-dsDNA.pdf",width = 10,height = 9)
drow_function(genomfeature_PQS,"G2","ssRNA(+)")
ggsave("G2-ssRNA_p.pdf",width = 10,height = 6)
drow_function(genomfeature_PQS,"G2","ssRNA(-)")
ggsave("G2-ssRNA_n.pdf",width = 10,height = 6)
drow_function(genomfeature_PQS,"G2","dsDNA-RT")
ggsave("G2-dsDNART.pdf",width = 10,height = 5)
drow_function(genomfeature_PQS,"G2","dsRNA")
ggsave("G2-dsRNA.pdf",width = 10,height = 6)
drow_function(genomfeature_PQS,"G2","ssDNA")
ggsave("G2-ssDNA.pdf",width = 10,height = 9)
drow_function(genomfeature_PQS,"G2","ssRNA-RT")
ggsave("G2-ssRNART.pdf",width = 10,height = 8)
drow_function(genomfeature_PQS,"G2","satellite")
ggsave("G2-satellite.pdf",width = 10,height = 5)
drow_function(genomfeature_PQS,"G2","viroid")
ggsave("G2-viroid.pdf",width = 10,height = 3)

drow_function(genomfeature_PQS,"G3","dsDNA")
ggsave("G3-dsDNA.pdf",width = 10,height = 8)
drow_function(genomfeature_PQS,"G3","ssRNA(+)")
ggsave("G3-ssRNA_p.pdf",width = 10,height = 5)
drow_function(genomfeature_PQS,"G3","ssRNA(-)")
ggsave("G3-ssRNA_n.pdf",width = 10,height = 5)
drow_function(genomfeature_PQS,"G3","dsDNA-RT")
ggsave("G3-dsDNART.pdf",width = 10,height = 2)
drow_function(genomfeature_PQS,"G3","dsRNA")
ggsave("G3-dsRNA.pdf",width = 10,height = 2)
drow_function(genomfeature_PQS,"G3","ssDNA")
ggsave("G3-ssDNA.pdf",width = 10,height = 6)
drow_function(genomfeature_PQS,"G3","ssRNA-RT")
ggsave("G3-ssRNART.pdf",width = 10,height = 6)

#drow_function(genomfeature_PQS,"G3","satellite")
#drow_function(genomfeature_PQS,"G3","viroid")

levels(factor(genomfeature_PQS$genomefeature))





