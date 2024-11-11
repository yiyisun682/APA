__version__="0.9.4"

import os
from optparse import OptionParser
usage = "usage: %prog -o arg1"
parser = OptionParser()
parser.add_option("-o", "--miRSitesbyTx", dest="miRSitesbyTx", help="[REQUIRED] name of output file having miRNA binding information for each transcript", metavar="FILE", default="miRSitesbyTx.txt")
(options, args) = parser.parse_args()

if not options.miRSitesbyTx:   # if filename is not given
    parser.error('options.miRSitesbyTx not given')
    
miRSitesbyTxFileName=options.miRSitesbyTx

#it returns overlap (+int, overlapStart, overlapEnd) or dist (-int, None, None)
def getOverlapOrDist(a, b):
    overlapOrDist=min(a[1], b[1]) - max(a[0], b[0])
    if overlapOrDist<0:
        return (overlapOrDist, None, None)
    else:
        return (overlapOrDist, max(a[0], b[0]), min(a[1], b[1]))
        
refSeqInfos,txsByGene,txsByGene={},{},{}
refSeqInfoFile=open("refGene_hg38.txt", "r")
header=refSeqInfoFile.readline()
for aRefSeqInfo in refSeqInfoFile:
    fields=aRefSeqInfo.rstrip("\n").split("\t")
    txName,chrID,strand,txStart,txEnd,cdsStart,cdsEnd,NCBIID=[fields[i] for i in [1,2,3,4,5,6,7,12]]
    if txName not in refSeqInfos:
        refSeqInfos[txName]=[(NCBIID, chrID,strand,int(txStart),int(txEnd),int(cdsStart),int(cdsEnd))]
    else:
        refSeqInfos[txName].append((NCBIID, chrID,strand,int(txStart),int(txEnd),int(cdsStart),int(cdsEnd)))
    if NCBIID not in txsByGene:
        txsByGene[NCBIID]=[]
    if txName not in txsByGene[NCBIID]:
        txsByGene[NCBIID].append(txName)
refSeqInfoFile.close()

TomiRFam={}
miRInfoFile=open("targetScan_miR_Family_Info.txt", "r")
header=miRInfoFile.readline()
for anMiR in miRInfoFile:
    fields=anMiR.rstrip("\n").split("\t")
    miRFam,speciesID,miRBaseID=[fields[i] for i in [0,2,3]]
    #assert speciesID=="9606", "this should be about HG19 (9606)"
    TomiRFam[miRBaseID]=miRFam
miRInfoFile.close()

miRsNotHighLow=[]
miRsNotHighLowFile=open("miRNA_detected_in_HSPC.txt", "r")
for aLine in miRsNotHighLowFile:
    miRsNotHighLow.append(aLine.rstrip("\n").split("\t")[0].replace('r', 'R'))
miRsNotHighLowFile.close()
miRsNotHighLow=set(miRsNotHighLow)
print(str(len(miRsNotHighLow)))

print("[targetScan]")
miRSitesbyGene,miRSitesbyTx,allMiRs,allInteractions={},{},[],[]#holding miR binding sites in each gene's 3UTR. 
miRbindingSitesFile=open("Predicted_Targets_Info.txt", "r")#in miR format
header=miRbindingSitesFile.readline()
for aBinding in miRbindingSitesFile:
    fields=aBinding.rstrip("\n").split("\t")
    miRFamID,GeneSymbol,transID,speciesID,UTRStart,UTREnd,PCT=[fields[i] for i in [0,2,3,4,5,6,10]]
    if speciesID!="9606":
        continue
    if transID not in refSeqInfos:
        continue
    if PCT=="NULL":
        continue
    if float(PCT)<=0:#
        continue
    if GeneSymbol not in txsByGene:
        continue
    if transID not in txsByGene[GeneSymbol]:
        continue
    
    NCBIID, chrID,strand,txStart,txEnd,cdsStart,cdsEnd=refSeqInfos[transID][0]
    if miRFamID not in miRsNotHighLow:
        continue
    if strand=="+":
        coordStart,coorEnd=cdsEnd+int(UTRStart),cdsEnd+int(UTREnd)
    else:
        coordStart,coorEnd=cdsStart-int(UTREnd),cdsStart-int(UTRStart)
    if transID not in miRSitesbyTx:
        miRSitesbyTx[transID]=[]
    miRSitesbyTx[transID].append((GeneSymbol,miRFamID,chrID,int(coordStart),int(coorEnd), strand))
    if GeneSymbol not in miRSitesbyGene:
        miRSitesbyGene[GeneSymbol]=[]
    miRSitesbyGene[GeneSymbol].append((miRFamID,chrID,int(coordStart),int(coorEnd), strand))
    allMiRs.append(miRFamID)
miRbindingSitesFile.close()

print("[tarBase]")
miRbindingSitesFile=open("TarBase_V5.0.txt", "r")#in miR format
header=miRbindingSitesFile.readline()
for aBinding in miRbindingSitesFile:
    fields=aBinding.rstrip("\n").split("\t")
    Support_Type,Organism,miRNA,Gene,ChrLoc=[fields[i] for i in [3,4,5,7,10]]
    if Organism!="Human":
        continue
    if Support_Type!="TRUE":
        continue
    if Gene not in txsByGene:
        continue
    if 'hsa-'+miRNA not in TomiRFam:
            continue
    miRFamID=TomiRFam['hsa-'+miRNA]
    if miRFamID not in miRsNotHighLow:
        continue
    chr,loc=ChrLoc.split(":")
    seedStart,seedEnd=[int(i.rstrip('\"').replace(",", "")) for i in loc.split("-")]
    if Gene not in miRSitesbyGene:
        miRSitesbyGene[Gene]=[]
    miRSitesbyGene[Gene].append((miRFamID,chr,int(seedStart),int(seedEnd), "?"))
    for tx4Gene in txsByGene[Gene]:
        NCBIID,chrID,strand,txStart,txEnd,cdsStart,cdsEnd=refSeqInfos[tx4Gene][0]
        if tx4Gene not in miRSitesbyTx:
            miRSitesbyTx[tx4Gene]=[(Gene,miRFamID,chr,int(seedStart),int(seedEnd), "?")]
            continue
        for geneID2Cmp, miRFamID2Cmp, chrID2Cmp, coordStart2Cmp, coorEnd2Cmp, strand2Cmp in miRSitesbyTx[tx4Gene]:
            if not ((getOverlapOrDist((coordStart2Cmp, coorEnd2Cmp), (int(seedStart),int(seedEnd)))[0]>0) and (miRFamID2Cmp==miRFamID)):#if it's not overlap
                allMiRs.append(miRFamID)
                if tx4Gene not in miRSitesbyTx:
                    miRSitesbyTx[tx4Gene]=[]
                miRSitesbyTx[tx4Gene].append((Gene,miRFamID,chr,int(seedStart),int(seedEnd), "?"))
                break
miRbindingSitesFile.close()

print("[miRecords]")
miRbindingSitesFile=open("miRecords_version4.txt", "r")#in miR format
header=miRbindingSitesFile.readline()
for aBinding in miRbindingSitesFile:
    fields=aBinding.rstrip("\n").split("\t")
    Organism,txID,miRName,reporter=[fields[i] for i in [5,3,6,8]]
    if Organism!="Homo sapiens":
        continue
    if ((reporter!="luciferase") and (reporter!="GFP")):#only taking luciferase and GFP
        continue
    if txID not in refSeqInfos:
        continue
    if miRName not in TomiRFam:
        continue
    miRFamID=TomiRFam[miRName]
    if miRFamID not in miRsNotHighLow:
        continue
    if txID not in miRSitesbyTx:
        miRSitesbyTx[txID]=[]
    if ((txID in miRSitesbyTx) and (miRFamID in [i[1] for i in miRSitesbyTx[txID]])):
        continue
    NCBIID,chrID,strand,txStart,txEnd,cdsStart,cdsEnd=refSeqInfos[txID][0]
    if NCBIID not in miRSitesbyGene:
        miRSitesbyGene[NCBIID]=[]
    miRSitesbyGene[NCBIID].append((miRFamID,chrID, -2,-2, "?"))
    miRSitesbyTx[txID].append((refSeqInfos[txID][0][0],miRFamID,chrID, -2,-2, "?"))
    allMiRs.append(miRFamID)
miRbindingSitesFile.close()

print("[miRTarBase]")
miRbindingSitesFile=open("hsa_MTI.txt", "r")#in miR format
header=miRbindingSitesFile.readline()
for aBinding in miRbindingSitesFile:
    fields=aBinding.rstrip("\n").split("\t")
    miRNA,species,targetGene,experiments=[fields[i] for i in [1,2,3,6]]
    if species!="Homo sapiens":
        continue
    if ((experiments=="Microarray") or (experiments=="Proteomics")):
        continue
    if miRNA not in TomiRFam:
        continue
    miRFamID=TomiRFam[miRNA]
    if miRFamID not in miRsNotHighLow:
        continue
    if targetGene not in txsByGene:
        continue
    if targetGene not in miRSitesbyGene:
        miRSitesbyGene[targetGene]=[]
    miRSitesbyGene[targetGene].append((miRFamID,chrID, -3,-3, "?"))
    for tx4Gene in txsByGene[targetGene]:
        if tx4Gene not in refSeqInfos:
            continue
        NCBIID,chrID,strand,txStart,txEnd,cdsStart,cdsEnd=refSeqInfos[tx4Gene][0]
        if tx4Gene not in miRSitesbyTx:
            miRSitesbyTx[tx4Gene]=[]
        if ((tx4Gene in miRSitesbyTx) and (miRFamID in [i[1] for i in miRSitesbyTx[tx4Gene]])):
            continue
        miRSitesbyTx[tx4Gene].append((targetGene,miRFamID,chrID,-3,-3,"?"))
        allMiRs.append(miRFamID)
miRbindingSitesFile.close()

print("[AGO-CLIP]")
AGO_CLIPFile=open("ncomms3730-s4.txt", "r")
header=AGO_CLIPFile.readline()
for aLine in AGO_CLIPFile:
    fields=aLine.rstrip("\n").split("\t")
    chrID,seedStart,seedEnd,strand,miRFamID,geneID,occurrence=[fields[i] for i in [0, 3, 4, 5, 6, 10, 11]]
    if int(occurrence)<3:#corresponding to q-val<0.05
        continue
    if miRFamID not in miRsNotHighLow:
        continue
    if geneID not in txsByGene:
        continue
    if geneID not in miRSitesbyGene:
        miRSitesbyGene[geneID]=[]
    miRSitesbyGene[geneID].append((miRFamID,chrID, int(seedStart),int(seedEnd), strand))
    for tx4Gene in txsByGene[geneID]:
        if tx4Gene not in miRSitesbyTx:
            miRSitesbyTx[tx4Gene]=[(geneID,miRFamID,chrID,int(seedStart),int(seedEnd), strand)]
            continue

        miRSet=set([(i[1],i[3],i[4]) for i in miRSitesbyTx[tx4Gene] if ((i[3]>0) and (i[4]>0))])
        for miRFamID2Cmp, coordStart2Cmp, coorEnd2Cmp in miRSet:
            if not ((getOverlapOrDist((coordStart2Cmp, coorEnd2Cmp), (int(seedStart),int(seedEnd)))[0]>0) and (miRFamID2Cmp==miRFamID)):#if it's not overlap
                allMiRs.append(miRFamID)
                if tx4Gene not in miRSitesbyTx:
                    miRSitesbyTx[tx4Gene]=[]
                miRSitesbyTx[tx4Gene].append((geneID,miRFamID,chrID,int(seedStart),int(seedEnd), strand))
                break
AGO_CLIPFile.close()

miRSitesbyGeneFile=open("miRSitesbyTx.txt", "w")
for txName in miRSitesbyTx:
    miRFamIDs=set([i[1] for i in miRSitesbyTx[txName]])
    for geneID,miRFamID,chrID,seedStart,seedEnd,strand in miRSitesbyTx[txName]:
        miRSitesbyGeneFile.write(txName+"\t"+geneID+"\t"+miRFamID+"\t"+chrID+"\t"+str(seedStart)+"\t"+str(seedEnd)+"\t"+strand+"\n")
miRSitesbyGeneFile.close()
