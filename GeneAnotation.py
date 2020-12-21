from collections import defaultdict
import re,os

class PQSinGene_coord():
    def __init__(self,gff,PQSfile=""):
        self.gff=gff
        self.PQSfile=PQSfile

    def anotate_nongene(self):
        anotation_dir = defaultdict(lambda: defaultdict(list))
        with open(self.gff, 'r') as F:
            for line in F:
                if line.startswith("##sequence-region"):
                    MatchKey, region_start, region_end = line.rstrip().split(" ")[1:]
                    anotation_dir[MatchKey]['region'] = [int(region_start), int(region_end)]
                    continue
                MatchKey = re.match(r'[A-Za-z]+?_?\d+\.?\d?', line, re.I)
                if MatchKey:
                    if line.split("\t")[2] == "CDS":
                        genename=re.findall('[\t;]gene=(.+?)[;\n]',line)
                        if genename:
                            genename=genename[0]
                        else:
                            genename = re.findall(r'Parent=(.+?);', line, re.I)[0] if "Parent=" in line else \
                                re.findall(r'Name=(.+?);', line, re.I)[0]
                        genename = genename.replace("gene-","").replace(' ','-')
                        strandtype = line.split("\t")[6]
                        start, end = [int(i) for i in line.split("\t")[3:5]]
                        anotation_dir[MatchKey[0]][genename].append((start, end, 'CDS', strandtype))
        for ncnumber in anotation_dir:
            if not [i for i in anotation_dir[ncnumber] if i!="region"]:
                continue
            left_limit, right_limit = anotation_dir[ncnumber]['region']
            codings=[[location[0:2] for location in anotation_dir[ncnumber][gene]]
                     for gene in anotation_dir[ncnumber] if gene != "region"]
            codings=sorted([b for a in codings for b in a], key=lambda X: (X[0], X[1]))
            # noncodings = [(codings[i][1] + 1, codings[i + 1][0] - 1) for i in range(len(codings))
            #              if len(codings) > i + 1 and codings[i + 1][0] > codings[i][1]]
            noncodings=[]
            for i in range(len(codings)):
                if len(codings) > i + 1 and codings[i + 1][0] > codings[i][1]:
                    maxvalue=max([a[1] for a in codings[:i+1]])
                    if maxvalue<codings[i + 1][0]:
                        noncodings.append([maxvalue+1, codings[i + 1][0] - 1])
            if codings:
                pemax=max([i[1] for i in codings])
                if codings[0][0]>left_limit:
                    noncodings.insert(0,(left_limit,codings[0][0]-1))
                if pemax<right_limit:
                    noncodings.append((pemax+1,right_limit))
            else:
                noncodings.insert(0,(left_limit,right_limit))
            anotation_dir[ncnumber]['noncoding']=[i for i in noncodings if i[0]<i[1]]
        return {nc:[(b[0],b[1],a,b[-1]) if a!="noncoding" else (b[0],b[1],a) for a in anotation_dir[nc]
                    for b in anotation_dir[nc][a] if a!="region"] for nc in anotation_dir }

os.chdir("/home/zhluo/Project/lizhen_lncRNA/Allvirus/G4ReFinder/")
nongene=PQSinGene_coord(gff="../Allvirus.gff3").anotate_nongene()
oplines=[]
for nc in nongene:
    genelist=sorted(nongene[nc],key=lambda x:x[0])
    i = 1
    for gene in genelist:
        if len(gene)==3:
            ano="".join(['noncoding','(',str(i),')'])
            i+=1
        else:
            ano="".join([gene[2],'(',gene[-1],')'])
        oplines.append("\t".join([nc, str(gene[0]), str(gene[1]), ano])+"\n")
with open("Allvirus_ano.bed",'w') as F:
    F.writelines(oplines)

