'''
Created on 23.08.2017

@author: marisakoe
'''

import pandas as pd
from numpy import *

from tcoffee import tCoffee

def compute_without_synonyms(data_list,tree, sounds, pmiScores, gp1, gp2 ):      
    ##create a pandas data frame out of the data
    wl = pd.DataFrame(data_list, columns=['doculect','ASJP'])
    ##drop the double entries in the doculect
    wl = wl.drop_duplicates('doculect')
      
    ##create a pmiDict
    pmiDict = dict()
    for s1 in sounds:
        for s2 in sounds:
            pmiDict[s1,s2] = pmiScores[s1].ix[s2]
    
    
    ##align the data with tCoffee
    alg = tCoffee(wl,tree,pmiDict,gp1,gp2, sounds)
    
    ## Data Frame for the global alignments
    algDF = pd.DataFrame([list(x) for x in alg.alignment.values],index=alg.doculect.values,columns = range(len(alg.alignment.values[0])))
    
    ##get the taxa
    taxa = alg.doculect.values
     
     
    ##create the binary matrix
    binMtx = pd.DataFrame(index=taxa)
    ##for each column in the global alignment in the dataframe
    for i in algDF.columns:
        ##get the alignments in column i for all doculects
        cl = algDF[i]
        ##get the individual states (sounds) excluding gaps
        states = unique([s for s in cl.values if s!='-'])
        ##create a new pandas frame for each colum, if doculect has the state assign 1 otherwise 0
        clBin = pd.DataFrame(array([(cl==s).values for s in states],int).T, index=cl.index)
        ##concatenate the matrices to a bigger one
        binMtx = pd.concat([binMtx,clBin],axis=1)
         
    
    pad = max(map(len,taxa))+5
    
    return pad, taxa,binMtx
        
def compute_with_synonyms(data_list,tree, sounds, pmiScores, gp1, gp2):
    ##create a pandas data frame out of the data
    wl = pd.DataFrame(data_list, columns=['doculect','ASJP'])
      
    ##create a pmiDict
    pmiDict = dict()
    for s1 in sounds:
        for s2 in sounds:
            pmiDict[s1,s2] = pmiScores[s1].ix[s2]
    
    
    ##align the data with tCoffee
    alg = tCoffee(wl,tree,pmiDict,gp1,gp2, sounds)

    ## Data Frame for the global alignments
    algDF = pd.DataFrame([list(x) for x in alg.alignment.values],index=alg.doculect.values,columns = range(len(alg.alignment.values[0])))
    ##get the taxa
    taxa = alg.doculect.values

    ##create the binary matrix
    binMtx = pd.DataFrame(index=taxa)
    
    ##for each column in the global alignment in the dataframe
    for i in algDF.columns:
        
        ##get the alignments in column i for all doculects
        cl = algDF[i]
        
        ##get the individual states (sounds) excluding gaps
        states = unique([s for s in cl.values if s!='-'])
        #print states
        ##create a new pandas frame for each column, if doculect has the state assign 1 otherwise 0
        clBin = pd.DataFrame(array([(cl==s).values for s in states],int).T, index=cl.index)

        ##concatenate the matrices to a bigger one
        binMtx = pd.concat([binMtx,clBin],axis=1)
         
    pad = max(map(len,taxa))+5
    
    return pad, taxa,binMtx

if __name__ == '__main__':
    pass