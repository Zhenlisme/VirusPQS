import os,re
import random,argparse
from collections import defaultdict

parser=argparse.ArgumentParser(description="To shuffle genomes and detect the GCPQS.")
parser.add_argument("-g","--fasta",required=True,type=str,help="The genome sequence,fasta format.")
parser.add_argument("-grun","--grun",required=True,type=int,help="The length of grun.")
parser.add_argument("-o","--opf",required=True,type=str,help="The output file name.")
Args=parser.parse_args()

class GCdetection():
    def __init__(self,genomefile):
        self.genome_dict={}
        self.seqdf=os.path.abspath(Args.opf)+".seq"
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
        GG = ",".join(re.findall('(?:G{%s}.{1,7}){3,}G{%s}'%(Args.grun,Args.grun),sequence, re.I))
        CC = ",".join(re.findall('(?:C{%s}.{1,7}){3,}C{%s}'%(Args.grun,Args.grun),sequence, re.I))
        return ','.join([GG,CC]).strip(',')
                                    
    def random_shuffle(self):
        oplines=[]
        for key in self.genome_dict:
            for i in range(1,1001):
                random_seq=self.shuffling(self.genome_dict[key])
                oplines.append("\t".join([key,str(i),self.detection(random_seq)])+'\n')
        return oplines
    def G4Hscore_cal(self,opdf):
        with open(self.seqdf,'w') as F:
            F.writelines(self.random_shuffle())
        os.system("Rscript /home/zhluo/Project/lizhen_lncRNA/Allvirus/G4ReFinder/random_g4hunter/g4hunterscore_random.r %s %s"%(self.seqdf,opdf))
        os.system("rm -f %s"%self.seqdf)
        

if __name__ == "__main__":
    main_process=GCdetection(genomefile=os.path.abspath(Args.fasta))
    main_process.G4Hscore_cal(opdf=os.path.abspath(Args.opf))
    
