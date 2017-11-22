'''
Created on Jun 13, 2017

@author: Guyling1
'''
from __future__ import division
import re 
import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt


def parseRealValuesInFile(filePath):
    dict={}
    effects=re.findall('P\d_*P*_*\d*__[ACTG]+_[ACTG]_-*\d+.\d*',filePath)
    motifs=[re.findall('P\d_*P*_*\d*__[ACTG]+_[ACTG]',eff) for eff in effects]
    motifs=["".join(m) for m in motifs]
    muts=[m[-1] for m in motifs]#indicates the mutation type
    motifs=[m[:-2] for m in motifs]#cleaving the mutation type from the effect
    fileMutationType=re.findall('_[ACTG][ACTG]_[ab]',filePath)[0][1:3]
    coeffs=[f.split('_')[-1] for f in effects]
    for i in range(len(motifs)):
        dict[(motifs[i],muts[i])]=coeffs[i]
    for key in dict.keys():#remove the effects not for this mutation type analysis
        if key[1]!=fileMutationType[1]:
            del dict[key]
    return dict

def readAnalysisFile(filePath,threshold=0.6):
    truePositive=0
    falsePositive=0
    falseNegative=0
    dataAnalysis=pd.read_csv(filePath)
    effects=parseRealValuesInFile(filePath)
    motifs=[ef[0] for ef in effects.keys()]
    #searching to see if the motifs were found. TP
    for motif in motifs:
        data=dataAnalysis.loc[dataAnalysis['name'] == motif]
        if data['gammaMean'].values[0]>=threshold:
            truePositive+=1
        else:
            falseNegative+=1
    allPos=dataAnalysis.loc[dataAnalysis['gammaMean'] >=threshold]['name'].values
    allNeg=dataAnalysis.loc[dataAnalysis['gammaMean'] <threshold]['name'].values
    for pos in allPos:#searching for positive-TP=FP
        if pos not in motifs:
            falsePositive+=1
    trueNegative=len(allNeg)-falseNegative
    return truePositive,falsePositive,falseNegative,trueNegative

def folderAnalyze(folder,): 
    os.chdir(folder)
    x=[]
    y=[]
    files=glob.glob("*bestChainsSummary.csv")
    for i in np.linspace(start=1,stop=0):
        print i
        results=[]
        for f in files:
            results.append(readAnalysisFile(os.path.join(folder,f),threshold=i))
        tp=[r[0] for r in results]
        fp=[r[1] for r in results]
        fn=[r[2] for r in results]
        tn=[r[3] for r in results]
        print "tp is :{}".format(sum(tp))
        print "fp is :{}".format(sum(fp))
        print "fn is :{}".format(sum(fn))
        print "tn is :{}".format(sum(tn))
        accuracy=(sum(tp)+sum(tn))/(sum(tp)+sum(tn)+sum(fp)+sum(fn))
        print "accuracy is {}".format(accuracy)
        print"sensitivity is {}".format(sum(tp)/(sum(tp)+sum(fn)))
        x.append(sum(tn)/(sum(tn)+sum(fp)))
        y.append(sum(tp)/(sum(tp)+sum(fn)))
        print "specificity is {}".format(sum(tn)/(sum(tn)+sum(fp)))
    plt.scatter(x, y)
    plt.title("ROC")
    plt.show()
        
    
    
print folderAnalyze(r'C:\Users\Guyling1\Documents\guyReserch\MCMCsimulationsResults\normdDet')