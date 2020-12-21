from collections import defaultdict
import os,re

def Coording(wigdir,anotationfile,opfile):
    wigdict={}
    for f in os.listdir(wigdir):
        nccore=f.split(".")[0]
        wigfile="/".join([wigdir,f])
        with open(wigfile,'r') as F:
            con_score=[float(line.rstrip()) for line in F if re.match('\d',line)]
            wigdict[nccore]=con_score
    keylist=tuple(map(lambda x:x.split(".")[0],wigdict.keys()))
    opf=open(opfile,'w')
    with open(anotationfile,'r') as F:
        for line in F:
            ncnumber,start,end,anotation=line.rstrip().split("\t")
            if ncnumber.split(".")[0] not in keylist:
                continue
            convalue=wigdict[ncnumber.split(".")[0]][int(start)-1:int(end)]
            value_list,coord_list=[],[]
            for idx in range(0,len(convalue),20):
                value=convalue[idx:idx+20]
                average=round(sum(value)/len(value),4)
                value_list.append(average)
                coord_list.append([idx+int(start),idx+len(value)+int(start)-1])
            try:
                if coord_list[-1][-1]-coord_list[-1][0]==0:
                    coord_list[-2][-1]=coord_list[-1][-1]
                    value_list[-2]=round((value_list[-2]*20+value_list[-1])/21,4)
                    del(coord_list[-1],value_list[-1])
                opf.writelines(list(map(lambda a,b:"\t".join([ncnumber,str(a[0]),str(a[1]),str(b),anotation])+"\n",coord_list,value_list)))
            except:
                print("\n",line,convalue,coord_list,sep="\n")
    opf.close()

if __name__ == '__main__':
    os.chdir("/home/zhluo/Project/lizhen_lncRNA/Allvirus/G4ReFinder/G4phastcons/")

    Coording(wigdir="./WigFile/",
             anotationfile="/home/zhluo/Project/lizhen_lncRNA/Allvirus/G4ReFinder/sorted_gene_numbered.bed",
             opfile="Allvirus_LSC_coord.xlsx")

