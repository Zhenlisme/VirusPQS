import os,re
import argparse
from collections import defaultdict

class CompareDetection():
    def __init__(self,genomefile):
        self.genome_dict={}
        with open(genomefile,'r') as F:
            for line in F:
                if line.startswith(">"):
                    key=line.split(" ")[0].strip(">\n")
                    self.genome_dict[key]=""
                else:
                    self.genome_dict[key]+=line.rstrip()
    def detection(self,key,sequence):
        opline=[]
        pattern = '(?:%s.{0,6}){3,}%s'
        GNG_count = re.finditer(pattern % ("G.G","G.G"), sequence, re.I)
        CNC_count = re.finditer(pattern % ("C.C","C.C"), sequence, re.I)
        GNGG_count = re.finditer(pattern % ("G.GG","G.GG"), sequence, re.I)
        CCNC_count = re.finditer(pattern % ("CC.C","CC.C"), sequence, re.I)
        GGNG_count = re.finditer(pattern % ("GG.G","GG.G"), sequence, re.I)
        CNCC_count = re.finditer(pattern % ("C.CC","C.CC"), sequence, re.I)

        GNG_opt=[",".join([key,"GNG",str(i.span()[0]+1),str(i.span()[1]),i.group()]) for i in GNG_count]
        CNC_opt=[",".join([key,"CNC",str(i.span()[0]+1),str(i.span()[1]),i.group()]) for i in CNC_count]
        GNGG_opt=[",".join([key,"GNGG",str(i.span()[0]+1),str(i.span()[1]),i.group()]) for i in GNGG_count]
        CCNC_opt=[",".join([key,"CCNC",str(i.span()[0]+1),str(i.span()[1]),i.group()]) for i in CCNC_count]
        GGNG_opt=[",".join([key,"GGNG",str(i.span()[0]+1),str(i.span()[1]),i.group()]) for i in GGNG_count]
        CNCC_opt=[",".join([key,"CNCC",str(i.span()[0]+1),str(i.span()[1]),i.group()]) for i in CNCC_count]

        opline.extend(GNG_opt)
        opline.extend(CNC_opt)        
        opline.extend(GNGG_opt)
        opline.extend(CCNC_opt)
        opline.extend(GGNG_opt)
        opline.extend(CNCC_opt)

        return "\n".join(opline)
    def RealDetection(self,opfile):
        oplines=[]
        for key in self.genome_dict:
            oplines.append(self.detection(key,self.genome_dict[key])+"\n")
        with open(opfile,'w') as F:
            F.writelines(oplines)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To get the location of comparision. First-writtern@20190729 by ZhenLi.")
    parser.add_argument("-f", "--genomefile", type=str,required=True,help="Please type in the file of genome file.")
    parser.add_argument("-o", "--opfile", required=True, type=str, help="The output file name.")
    Args = parser.parse_args()
    D=CompareDetection(os.path.abspath(Args.genomefile))
    D.RealDetection(os.path.abspath(Args.opfile))


