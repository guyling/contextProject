'''
Created on Jan 16, 2017

@author: Guyling1
'''
import evoModel
import simulations
from simulations import NUCS
import sys
print "ok"
import time
import pandas as pd
import numpy as np
import numpy.random as nr
from featureExtraction import featureExtraction
global NUCS
NUCS=['A','C','G','T']
global NUCSNOC
NUCSNOC=['A','G','T']
import glob
import os
from GammaPriorTest import fullAnalysis
"""
mod=evoModel.evoModel('JC',rateParam={'mu':5*10**-5})
from SimulationGenerator import createSimulationRandomly,analyseSample
i=int(sys.argv[1])
createSimulationRandomly(oldNuc=i%4)
"""

"""
#analysis of multiple files in the foler location
location='/sternadi/home/volume1/guyling/MCMC/dataSimulations/dataSets/'
os.chdir(location)
i=int(sys.argv[1])
print "index is : "+str(i)
files=glob.glob("*.csv")
filesMut=[(f,NUCSNOC[j]) for f in files for j in range(3)]
fileMut=filesMut[i]
curFile,curMutType=fileMut[0],'C'+fileMut[1]
print fileMut
fullAnalysis(location+curFile,curMutType)
"""
"""
sim=simulations.EvolustionSimulations(mod,time=15,sampleNum=1000,chroLen=800)
sim.parseChromosome(r'C:\Users\Guyling1\ContextProject\AggarwalaPaper\sabin2.full_genome.U882C.A2973G.C4905U.C5526U.fasta')
sim.addMotif(([-1,-1,1,2,-1],3),5)
#sim.addMotif(([-1,1,1,1,-1],3),-2)   
#sim.addMotif(([3,-1,1,-1,-1],3),1) 
sim.initializeProb()
#sim.setName('testRun')
tic=time.time()
#sim.evolve()
moranMatrix=sim.moranModelByPosition()
sim.toFreqsFile(moranMatrix, r'C:\Users\Guyling1\Documents\guyReserch\moranFreqs.freqs')
fc=featureExtraction(r'C:\Users\Guyling1\Documents\guyReserch\moranFreqs.freqs',[(748,7371)],5)
fc.createRegressionTable()  
fc.regTable.to_csv(r'C:\Users\Guyling1\Documents\guyReserch\moranFreqsTry.csv')
#sim.ouputTofasta(r'C:\Users\Guyling1\Documents\guyReserch\outputSimulator.txt')
"""
 

"""
#createSimulationRandomly(index,numOfMotifs=int(i/20.)+1)



"""
"""
if newNuc==1:
    newNuc=3
location='/sternadi/home/volume1/guyling/MCMC/dataSimulations'
file="{}motif_{}.csv".format(motifNumber,simulationNumber)
analyseSample(file,newNuc=newNuc)

createSimulationRandomly("1motif_{}".format(i),numOfMotifs=1)
createSimulationRandomly("2motif_{}".format(i),numOfMotifs=2)
createSimulationRandomly("3motif_{}".format(i),numOfMotifs=3) 

     
sim=simulations.EvolustionSimulations(mod,time=15,sampleNum=10000,chroLen=800)
sim.parseChromosome(r'C:\Users\Guyling1\ContextProject\AggarwalaPaper\sabin2.full_genome.U882C.A2973G.C4905U.C5526U.fasta')
sim.addMotif(([-1,-1,1,2,-1],3),2)
sim.addMotif(([-1,1,1,1,-1],3),-2)   
sim.addMotif(([3,-1,1,-1,-1],3),1)  
#sim.parseChromosome(r'C:\Users\Guyling1\Documents\guyReserch\sabin.txt')
sim.initializeProb()
#sim.setName('testRun')
tic=time.time()
#sim.evolve()
moranMatrix=sim.moranModelByPosition()
sim.toFreqsFile(moranMatrix, r'C:\Users\Guyling1\Documents\guyReserch\moranFreqs.freqs')
fc=featureExtraction(r'C:\Users\Guyling1\Documents\guyReserch\moranFreqs.freqs',[(748,7371)],5)
fc.createRegressionTable()  
fc.regTable.to_csv(r'C:\Users\Guyling1\Documents\guyReserch\moranFreqsTry.csv')
#sim.ouputTofasta(r'C:\Users\Guyling1\Documents\guyReserch\outputSimulator.txt')
tac=time.time()
print tac-tic
fc=featureExtraction(r'/media/jenia/ESD-USB/AggarwalaPaper/SABIN/P7.freqs',[(748,7371)],5)
fc.createRegressionTable()  
fc.regTable.to_csv(r'/media/jenia/ESD-USB/AggarwalaPaper/SABIN/P7_regressionFile.csv')
"""
global MT
MT=['AG','CT','GA''TC']
i=int(sys.argv[1])
fullAnalysis(r'/sternadi/home/volume1/guyling/MCMC/OPV/dataset/P2_regressionFile.csv', MT[i])
