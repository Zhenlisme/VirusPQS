#This script was used for calculating the G4score for each g4seq.

###extract from the program G4Hunter
###################################################################
library(stringi)
require(S4Vectors)

args=commandArgs(T)

G4translate <- function(x)		
  # x is a Rle of a sequence
{
  xres=x
  runValue(xres)[runValue(x)=='C' & runLength(x)>3] <- -4
  runValue(xres)[runValue(x)=='C' & runLength(x)==3] <- -3
  runValue(xres)[runValue(x)=='C' & runLength(x)==2] <- -2
  runValue(xres)[runValue(x)=='C' & runLength(x)==1] <- -1
  runValue(xres)[runValue(x)=='G' & runLength(x)>3] <- 4
  runValue(xres)[runValue(x)=='G' & runLength(x)==3] <- 3
  runValue(xres)[runValue(x)=='G' & runLength(x)==2] <- 2
  runValue(xres)[runValue(x)=='G' & runLength(x)==1] <- 1
  runValue(xres)[runValue(x)!='C' & runValue(x)!='G'] <- 0
  Rle(as.numeric(xres))
}

G4Hscore <- function(y)		
  # y can be DNAString or a DNAStringSet or a simple char string.
{
    y2 <- Rle(strsplit(as.character(y),NULL)[[1]])
    y3 <- G4translate(y2)
    return(abs(mean(y3)))
  }

inpdf=read.csv2(args[1],sep = '\t',stringsAsFactors = F,header = F)

scoredf=tapply(inpdf$V3,list(inpdf$V1,inpdf$V2),function(x){
  if(stri_length(x)==0){score=0}
  else{score=mean(sapply(strsplit(x,split = ',')[[1]],G4Hscore))}
  return(round(score,3))}
  )

write.table(as.data.frame(scoredf),file = args[2],row.names = T,col.names = F,sep = ',',quote=F)
