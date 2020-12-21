import os,re
import random,argparse
from collections import defaultdict

parser=argparse.ArgumentParser(description="To shuffle genomes and detect the GCPQS.")
parser.add_argument("-g","--fasta",required=True,type=str,help="The genome sequence,fasta format.")
parser.add_argument("-m","--method",required=True,type=str,choices=["shuffle","detect"],help="To shuffle the genome or not before you detect the sequence")
parser.add_argument("-o","--opf",required=True,type=str,help="The output file name.")
parser.add_argument("-l","--length",action='store_true',help="To calculate the g4 length.")
Args=parser.parse_args()

class GCdetection():
    def __init__(self,genomefile):
        self.genome_dict={}
        with open(genomefile,'r') as F:
            for line in F:
                if line.startswith(">"):
                    key=line.split(" ")[0].strip(">\n")
                    self.genome_dict[key]=""
                else:
                    self.genome_dict[key]+=line.rstrip().upper()
    def shuffling(self,sequence):
        seqlist=list(sequence)
        random.shuffle(seqlist)
        return "".join(seqlist)
    def detection(self,sequence):
        GG_count = len(re.findall('(?:G{2}.{1,7}){3,}G{2}.',sequence, re.I))
        CC_count = len(re.findall('(?:C{2}.{1,7}){3,}C{2}.',sequence, re.I))
        GNG_count = len(re.findall('(?:G.G.{0,6}){3,}G.G', sequence, re.I))
        CNC_count = len(re.findall('(?:C.C.{0,6}){3,}C.C', sequence, re.I))
        GGG_count = len(re.findall('(?:G{3}.{1,7}){3,}G{3}.',sequence, re.I))
        CCC_count = len(re.findall('(?:C{3}.{1,7}){3,}C{3}.',sequence, re.I))
        GNGG_count = len(re.findall('(?:G.GG.{0,6}){3,}G.GG', sequence, re.I))
        CCNC_count = len(re.findall('(?:CC.C.{0,6}){3,}CC.C', sequence, re.I))
        GGNG_count = len(re.findall('(?:GG.G.{0,6}){3,}GG.G', sequence, re.I))
        CNCC_count = len(re.findall('(?:C.CC.{0,6}){3,}C.CC', sequence, re.I))
        return  ",".join([str(GG_count),str(CC_count),str(GNG_count), str(CNC_count),str(GGG_count),str(CCC_count),str(GNGG_count), str(CCNC_count),
                  str(GGNG_count), str(CNCC_count)])
                                    
    def G4length(self,sequence):
        GG_length = [i.end()-i.start() for i in re.finditer('(?:G{2}.{1,7}){3,}G{2}.',sequence, re.I)]
        CC_length = [i.end()-i.start() for i in re.finditer('(?:C{2}.{1,7}){3,}C{2}.',sequence, re.I)]
        GNG_length = [i.end()-i.start() for i in re.finditer('(?:G.G.{0,6}){3,}G.G', sequence, re.I)]
        CNC_length = [i.end()-i.start() for i in re.finditer('(?:C.C.{0,6}){3,}C.C', sequence, re.I)]        
        GGG_length = [i.end()-i.start() for i in re.finditer('(?:G{3}.{1,7}){3,}G{3}.',sequence, re.I)]
        CCC_length = [i.end()-i.start() for i in re.finditer('(?:C{3}.{1,7}){3,}C{3}.',sequence, re.I)]
        GNGG_length = [i.end()-i.start() for i in re.finditer('(?:G.GG.{0,6}){3,}G.GG', sequence, re.I)]
        CCNC_length = [i.end()-i.start() for i in re.finditer('(?:CC.C.{0,6}){3,}CC.C', sequence, re.I)] 
        GGNG_length = [i.end()-i.start() for i in re.finditer('(?:GG.G.{0,6}){3,}GG.G', sequence, re.I)]
        CNCC_length = [i.end()-i.start() for i in re.finditer('(?:C.CC.{0,6}){3,}C.CC', sequence, re.I)]
        return ",".join(str(round(sum(i)/len(i),3)) if len(i)!=0 else '0' for i in [GG_length,CC_length,GNG_length,CNC_length,GGG_length,CCC_length,GNGG_length,CCNC_length,GGNG_length,CNCC_length])
        
    def random_shuffle(self,opfile):
        random_dict=defaultdict(list)
        for i in range(1000):
            for key in self.genome_dict:
                random_seq=self.shuffling(self.genome_dict[key])
                if Args.length:
                    random_dict[key].append(self.G4length(random_seq))
                else:
                    random_dict[key].append(self.detection(random_seq))
        with open(opfile,'w') as F:
            F.write("ncnumber,GG_count,CC_count,GNG_count,CNC_count,GGG_count,CCC_count,GNGG_count,CCNC_count,GGNG_count,CNCC_count\n")
            for key in random_dict:
                [F.write(",".join([key,t])+"\n") for t in random_dict[key]]
    def RealDetection(self,opfile):
        oplines=[]
        oplines.append("ncnumber,GG_count,CC_count,GNG_count,CNC_count,GGG_count,CCC_count,GNGG_count,CCNC_count,GGNG_count,CNCC_count\n")
        for key in self.genome_dict:
            oplines.append(",".join([key,self.detection(self.genome_dict[key])])+"\n")
        with open(opfile,'w') as F:
            F.writelines(oplines)
if __name__ == "__main__":
    main_process=GCdetection(genomefile=os.path.abspath(Args.fasta))
    if Args.method=="shuffle":
        main_process.random_shuffle(opfile=os.path.abspath(Args.opf))
    else:
        main_process.RealDetection(opfile=os.path.abspath(Args.opf))
    
