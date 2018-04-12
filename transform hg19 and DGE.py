# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 14:31:42 2017

@author: bmaier
"""

from pylab import array
import os
import pandas as pd
import datetime

pd.options.mode.chained_assignment = None  # default='warn'

#get file names and paths and set output path
inputFile_hg19="hg19_geneMapping_forESAT.txt"
inputFile_DGE="norm_reads_window_DGE_christoph.csv"
outputFile_hg19="hg19andreadpos_transformed.txt"
outputFile_DGE="norm_reads_window_DGE_christoph_indexed.txt"

currentDir=os.path.dirname(os.path.abspath(__file__))
inputFileandPath_hg19=os.path.join(currentDir,'Input',inputFile_hg19)
inputFileandPath_DGE=os.path.join(currentDir,'Input',inputFile_DGE)

ouputFileandPath_hg19=os.path.join(currentDir,'Output',outputFile_hg19)
outputFileandPath_DGE=os.path.join(currentDir,'Output',outputFile_DGE)

#open files and declare and preallocate some variables
TexonStartsList=pd.Series()
TexonEndsList=pd.Series()
    #load hg19 data
hg19table =  pd.read_csv(inputFileandPath_hg19, sep="\t", header=0)
nums=list(range(1,len(hg19table)+1))
nonums=[]
Index=pd.Series(nums,name='Index')

exonLength=pd.Series(nonums,name='exonLength')
intronLength=pd.Series(nonums,name='intronLength')
readStart=pd.Series(nonums,name='readStart')
readEnd=pd.Series(nonums,name='readEnd')
readItem=pd.Series(nonums,name='readItem')
hg19tableI=pd.concat([Index,hg19table.loc[:,:"exonEnds"], exonLength,intronLength, readStart, readEnd, readItem, hg19table.loc[:,"score":]],axis=1)

DGEtable =  pd.read_csv(inputFileandPath_DGE, sep=",", header=0)  
IndexDGE=pd.Series(list(range(0,len(DGEtable))),name='Index')
DGEtableI = pd.concat([IndexDGE, DGEtable],axis=1)

procisof=0

for element in hg19tableI.Index:
    procisof+=1
    
    TexonStartsList=""
    TexonEndsList=""
    TexonLengthList=""
    TintronLengthList=""
    readStartList=""
    readEndList=""
    readItemList=""

    strand=array(hg19tableI.strand[hg19tableI.Index==element], dtype=str)[0]
    exonCount=hg19tableI.exonCount[hg19tableI.Index==element]
    symbol=array(hg19tableI.name2[hg19tableI.Index==element], dtype=str)[0]
    
    if strand == "+":
        txStart=int(hg19tableI.txStart[hg19tableI.Index==element])-1
        hg19tableI.cdsStart[hg19tableI.Index==element]=int(hg19tableI.cdsStart[hg19tableI.Index==element])-txStart
        hg19tableI.cdsEnd[hg19tableI.Index==element]=int(hg19tableI.cdsEnd[hg19tableI.Index==element])-txStart
        exonStarts=hg19tableI.exonStarts[hg19tableI.Index==element]
        exonStartsList=exonStarts.str.split(',', expand=True)
        exonEnds=hg19tableI.exonEnds[hg19tableI.Index==element]
        exonEndsList=exonEnds.str.split(',', expand=True)
        hg19tableI.txEnd[hg19tableI.Index==element]=int(hg19tableI.txEnd[hg19tableI.Index==element])-txStart
        
        for i in range(exonCount):
            TexonStart=int(exonStartsList[i])-txStart
            TexonEnd=int(exonEndsList[i])-txStart
            TexonLength=TexonEnd-TexonStart
            if i< max(range(exonCount))-1:
                TintronLength=int(exonStartsList[i+1])-int(exonEndsList[i])
                TintronLengthList +=str(TintronLength)+","
                TexonLengthList +=str(TexonLength)+","
                TexonStartsList +=str(TexonStart)+","
                TexonEndsList +=str(TexonEnd)+","
            if i== max(range(exonCount)):

                TexonLengthList +=str(TexonLength)+","
                TexonStartsList +=str(TexonStart)+","
                TexonEndsList +=str(TexonEnd)+","
            
        readStartCol=DGEtable.start[DGEtable.Symbol==symbol]-txStart
        readEndCol=DGEtable.end[DGEtable.Symbol==symbol]-txStart
        for k in readStartCol:
            if k > 0:
                readStartList +=str(k)+","
                readItemList +=str(readStartCol[readStartCol==k].index[0])+"," 
        for k in readEndCol:
            if k > 0:
                readEndList +=str(k)+","
             
    else:
        txStart=int(hg19tableI.txEnd[hg19tableI.Index==element])+1
        cdsStart=txStart-int(hg19tableI.cdsEnd[hg19tableI.Index==element])
        hg19tableI.cdsEnd[hg19tableI.Index==element]=txStart-int(hg19tableI.cdsStart[hg19tableI.Index==element])
        hg19tableI.cdsStart[hg19tableI.Index==element]=cdsStart
        exonStarts=hg19tableI.exonEnds[hg19tableI.Index==element]
        exonStartsList=exonStarts.str.split(',', expand=True)
        exonEnds=hg19tableI.exonStarts[hg19tableI.Index==element]
        exonEndsList=exonEnds.str.split(',', expand=True)
        hg19tableI.txEnd[hg19tableI.Index==element]=txStart-int(hg19tableI.txStart[hg19tableI.Index==element])
        for i in range(exonCount-1,-1,-1):
            TexonStart=txStart-int(exonStartsList[i])
            TexonEnd=txStart-int(exonEndsList[i])
            TexonLength=TexonEnd-TexonStart
            if i>= 1:
                TintronLength=int(exonEndsList[i])-int(exonStartsList[i-1])
                TintronLengthList +=str(TintronLength)+","
                TexonStartsList +=str(TexonStart)+","
                TexonEndsList +=str(TexonEnd)+","
                TexonLengthList +=str(TexonLength)+","
            if i == 0:
                TexonStartsList +=str(TexonStart)+","
                TexonEndsList +=str(TexonEnd)+","
                TexonLengthList +=str(TexonLength)+","
                
        readStartCol=txStart-DGEtable.end[DGEtable.Symbol==symbol]
        readEndCol=txStart-DGEtable.start[DGEtable.Symbol==symbol]
        for k in readStartCol:
            if k > 0:
                readStartList +=str(k)+","
                readItemList +=str(readStartCol[readStartCol==k].index[0])+"," 
        for k in readEndCol:
            if k> 0:
                readEndList +=str(k)+","
     
    hg19tableI.exonStarts[hg19tableI.Index==element]=TexonStartsList.strip(',')
    hg19tableI.exonEnds[hg19tableI.Index==element]=TexonEndsList.strip(',')
    hg19tableI.exonLength[hg19tableI.Index==element]=TexonLengthList.strip(',')
    hg19tableI.intronLength[hg19tableI.Index==element]=TintronLengthList.strip(',')
    hg19tableI.readStart[hg19tableI.Index==element]=readStartList.strip(',')
    hg19tableI.readEnd[hg19tableI.Index==element]=readEndList.strip(',')
    hg19tableI.readItem[hg19tableI.Index==element]=readItemList.strip(',')
    hg19tableI.txStart[hg19tableI.Index==element]=1
    
    

    if procisof==100:
            print('At: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "  " +str(element) +" elements processed so far")
            #print str(element)+" elements processed so far" 
            procisof = 0


###############################################################
hg19tableI.to_csv(ouputFileandPath_hg19,sep='\t',index=False)
DGEtableI.to_csv(ouputFileandPath_DGE, sep='\t', index=False) #outcomment this or test reasons#
###############################################################
###############################################################
            #break           #uncomment this for test reasons#
###############################################################    
#next step would be to splice the DNA - this should then be stored in a new table