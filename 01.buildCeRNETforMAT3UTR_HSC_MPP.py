__version__="0.9.3"

import numpy as np
import os,glob
import os.path
import sys
from scipy.stats.stats import pearsonr
from optparse import OptionParser

import math
import scipy.misc as sc

#miRNA binding sites overlap pvalue
def calOverlapPval(sizeN,sizeA,sizeB,observationl):
    mostPossiblel=min(sizeA,sizeB)
    allProbsInLog=[math.log(sc.comb(sizeA, l, exact=True))+math.log(sc.comb((sizeN-sizeA), (sizeB-l), exact=True))-math.log(sc.comb(sizeN,sizeB)) for l in range(observationl, mostPossiblel+1)]
    return sum([math.exp(i) for i in allProbsInLog])

#return overlap(+int,overlapStartlocation,overlapEndlocation) or dist(no overlap,distance:-int)
def getOverlapOrDist(a, b):
    overlapOrDist=min(a[1], b[1]) - max(a[0], b[0])
    if overlapOrDist<0:
        return (overlapOrDist, None, None)
    else:
        return (overlapOrDist, max(a[0], b[0]), min(a[1], b[1]))

CeRNACorrFileName="ceRNACorr_HSC_MPP.txt"
ceRNACutoff=0.6 #pearson correlation

#source code
usage = "usage: %prog -a arg1 -b arg2 -c arg3 -d arg4 -e arg5 -f arg6 -g arg7 -i arg8 -j arg9 -o arg10"
parser = OptionParser()
#output
parser.add_option("-a", "--threeUTRInfo", dest="threeUTRInfo", help="[REQUIRED] file with gene models", metavar="FILE")
parser.add_option("-b", "--geneExprNormal", dest="geneExprNormal", help="[REQUIRED] file with gene expression", metavar="FILE")
parser.add_option("-c", "--miRSitesbyTx", dest="miRSitesbyTx", help="[REQUIRED] file with miRNA binding site information", metavar="FILE")
parser.add_option("-d", "--miRToInclude", dest="miRToInclude", help="[OPTIONAL] file with miRNA to be used", metavar="FILE", default="")
parser.add_option("-e", "--threeUS", dest="threeUS", help="[REQUIRED] file with 3UTR shortening information", metavar="FILE")
#input
parser.add_option("-o", "--APA_ceRNA", dest="APA_ceRNA", help="[REQUIRED] outfile with 3UTR shortening genes and their ceRNA partners", metavar="FILE")

(options, args) = parser.parse_args()
if not options.threeUTRInfo:   # if filename is not given
    parser.error('options.threeUTRInfo not given')
if not options.geneExprNormal:   # if filename is not given
    parser.error('options.geneExprNormal not given')
if not options.miRSitesbyTx:   # if filename is not given
    parser.error('options.miRSitesbyTx not given')
if not options.threeUS:   # if filename is not given
    parser.error('options.threeUS not given')
#miRToInclude??
if not options.APA_ceRNA:   # if filename is not given
    parser.error('options.ceRNA not given')

#input    
threeUTRInfoFile=options.threeUTRInfo
miRSitesbyTxFile=options.miRSitesbyTx
miRToIncludeFile=options.miRToInclude
geneExprNormal=options.geneExprNormal
threeUS=options.threeUS
#output
outAPA_ceRNAFileName=options.APA_ceRNA

#gene model information
refSeqInfos={}
refSeqInfoFile=open(threeUTRInfoFile, "r")
header=refSeqInfoFile.readline()
for aRefSeqInfo in refSeqInfoFile:
    txName,NCBIID,chrID,threeUTR_Start,threeUTR_End=aRefSeqInfo.rstrip("\n").split("\t")

    if txName not in refSeqInfos:
        refSeqInfos[txName]=[(NCBIID, chrID,threeUTR_Start,threeUTR_End)]
    else:
        refSeqInfos[txName].append((NCBIID, chrID,threeUTR_Start,threeUTR_End))
refSeqInfoFile.close()

#miRNA used
miRsNotHighLow=[]
if miRToIncludeFile!="":
    inFile=open(miRToIncludeFile, "r")
    for aLine in inFile:
        fields=aLine.rstrip("\n").split("\t")
        miRsNotHighLow.append(fields[0])
    inFile.close()

#miRNA binding sites for each transcript
miRSitesbyTx,allMiRs={},[]#holding miR binding sites in each gene's 3UTR. 
inFile=open(miRSitesbyTxFile, "r")
for aLine in inFile:
    fields=aLine.rstrip("\n").split("\t")
    if fields[0] not in miRSitesbyTx:
        miRSitesbyTx[fields[0]]=[]
    miRSitesbyTx[fields[0]].append(fields[2])
    allMiRs.append(fields[2])
inFile.close()

#ExprsNormal = CDSExprsNormal
CDSExprsNormal={}
exonExprFile=open(geneExprNormal, "r")
for aLine in exonExprFile:
    fields=aLine.rstrip("\n").split("\t")
    gene=fields[0]
    if gene.startswith("?"):
        continue
    CDSExprsNormal[gene]=[float(i) for i in fields[1:]]
exonExprFile.close()

print("[1. MAKING OVERLAP BETWEEN TRANSCRIPT PAIRS]")
import itertools
numAllPairs=0
list2Sort=[]
outFile=open("txPairWMirSiteShare0_HSC_MPP.txt", "w")
for txPair in itertools.combinations(miRSitesbyTx, r=2):
    if txPair[0] not in refSeqInfos:
        continue
    if txPair[1] not in refSeqInfos:
        continue
        
    miRs4Tx1=set(miRSitesbyTx[txPair[0]])
    miRs4Tx2=set(miRSitesbyTx[txPair[1]])
    if ((len(miRs4Tx1)<=5) or (len(miRs4Tx2)<=5)):
        continue
    commonMiRs=miRs4Tx1&miRs4Tx2
    numAllPairs+=1
    if len(commonMiRs)>0:
        sigOverlap=calOverlapPval(len(allMiRs), int(len(miRs4Tx1)), int(len(miRs4Tx2)), int(len(commonMiRs)))
        list2Sort.append((txPair[0], txPair[1], sigOverlap))
        outFile.write(refSeqInfos[txPair[0]][0][0]+"\t"+refSeqInfos[txPair[1]][0][0]+"\t"+txPair[0]+"\t"+txPair[1]+"\t"+str(int(len(miRs4Tx1)))+"\t"+str(len(miRs4Tx2))+"\t"+str(sigOverlap)+"\n")
outFile.close()

listSorted=sorted(list2Sort, key=lambda elmnt:elmnt[2])
alpha=0.05
numSigPairs=0
k=1
for tx1,tx2,sigOverlap in listSorted:
    if sigOverlap>(k*alpha)/float(numAllPairs):
        numSigPairs=k-1
        break
    k+=1
if numSigPairs==0:
    numSigPairs=k-1
    
print("[2. BUILDING CERNA NETWORK]")
from scipy.stats.stats import pearsonr

lineIdx=0
corrOut=open(CeRNACorrFileName, "w")
corrBWGenes=[]
for tx1,tx2,sigOverlap in listSorted:
    if lineIdx>numSigPairs:
        corrOut.close()
        break
    #tx1,tx2,lenMiRs4Gene1,lenMiRs4Gene2,lenCommonMiRs,sigOverlap=aLine.rstrip("\n").split(",")
    if ((tx1 not in refSeqInfos) or (tx2 not in refSeqInfos)):
        lineIdx+=1
        continue
    gene1,gene2=refSeqInfos[tx1][0][0],refSeqInfos[tx2][0][0]
    if gene1==gene2:
        lineIdx+=1
        continue
        
    if ((gene1 not in CDSExprsNormal) or (gene2 not in CDSExprsNormal)):
        lineIdx+=1
        continue
    corrmRNANormal=pearsonr(CDSExprsNormal[gene1], CDSExprsNormal[gene2])
    corrOut.write(",".join([gene1,gene2,tx1,tx2,str(corrmRNANormal[0])])+"\n")
    corrBWGenes.append((gene1,gene2,tx1,tx2,corrmRNANormal[0]))
    lineIdx+=1
corrOut.close()

print("[3. COLLECTING CERNA PARTNERS OF 3'UTR SHORTENING GENES]")#the output file will be used for parameter values of -e in MAT3UTR
sigAPAGenes=[]
sigAPAGenesFile=open(threeUS, "r")
header=sigAPAGenesFile.readline()
for aLine in sigAPAGenesFile:
    fields=aLine.rstrip("\n").split("\t")
    sigAPAGenes.append(refSeqInfos[fields[0]][0][0])
sigAPAGenesFile.close()
sigAPAGenes=set(sigAPAGenes)

sigAPACeRNAPartners={}
for gene1,gene2,tx1,tx2,corrmRNANormal in corrBWGenes:
    if corrmRNANormal<ceRNACutoff:
        continue
    isGene1APA=(gene1 in sigAPAGenes)
    isGene2APA=(gene2 in sigAPAGenes)
        
    if isGene1APA==True and isGene2APA==True:
        continue
    
    sigAPA,ceRNAPartner="",""
    if isGene1APA==True:
        sigAPA=gene1
        ceRNAPartner=gene2
    elif isGene2APA==True:
        sigAPA=gene2
        ceRNAPartner=gene1
    else:
        continue

    if ceRNAPartner not in sigAPACeRNAPartners:
        sigAPACeRNAPartners[ceRNAPartner]=[]
    sigAPACeRNAPartners[ceRNAPartner].append(sigAPA)

outFile=open(outAPA_ceRNAFileName, "w")
for ceRNAPartner in sigAPACeRNAPartners:
    outFile.write(ceRNAPartner+"\t"+",".join(set(sigAPACeRNAPartners[ceRNAPartner]))+"\n")
outFile.close()
