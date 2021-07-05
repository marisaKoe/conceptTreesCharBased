'''
Created on 17.08.2017

@author: marisakoe
'''

import subprocess, os, codecs, glob, shutil

import pandas as pd
from collections import defaultdict
from ete3 import Tree
from numpy import *

from nelex import compute_without_synonyms, compute_with_synonyms



def read_isoLangs(list_leaves):
    '''
    read the matching from the iso code to the language names
    '''
    f = open("input/isoLangNames.csv")
    raw_data = f.readlines()
    f.close()
    
    ##dictionary with key = iso value = lang
    isoLang_dict = defaultdict()
    for line in raw_data:
        line = line.strip().split(",")
        if line[1] in list_leaves:
            isoLang_dict[line[0]] = line[1]
        
    return isoLang_dict


def read_tree():
    '''
    read the tree in ete3 format and extracts the number of leaves for further analysis
    '''
    ##guide tree = language tree automatically computed with pmi
    tree = Tree('input/pmiTree2.tre')
    
    list_leaves = []
    for leaf in tree.iter_leaves():
        list_leaves.append(leaf.name)
        
    return list_leaves, tree

def read_nelex(list_leaves):
    '''
    read the nelex database
    '''
    
    f = open("input/northeuralex-cldf.tsv")
    raw_data = f.readlines()
    f.close()
    
    ##get the dictionary with the mappings
    isoLangdict = read_isoLangs(list_leaves)
    print len(list_leaves)
    ##dictionary with key = concept value=tuple with (lang,asjp_word)
    nelex_dict = defaultdict(list)
    unique_langs_concept = defaultdict(list)
    for l in raw_data[1:]:
        line = l.split("\t")
        ##glottocode
        #glot = line[0]
        ##iso code
        iso_code = line[1]
        
        ##concept
        concept = line[2]
        if " " in concept:
            concept = concept.replace(" ","")
        ##ipa
        #ipa = line[3]
        ##asjp
        asjp_word = line[4]
        ##sca
        #sca = line[5]
        ##dolgo
        #dolgo = line[6]
        
        ##ignoring the languages, which do not have an isocode in the mapping, they are therefore also not in the tree
        if iso_code in isoLangdict:
            nelex_dict[concept].append((isoLangdict[iso_code],asjp_word))
            if not isoLangdict[iso_code] in unique_langs_concept[concept]:
                unique_langs_concept[concept].append(isoLangdict[iso_code])
        

    return nelex_dict, unique_langs_concept
        

def write_dataMatrix_withoutSynonyms(fout1, dataMatrix,numLangs, numSoundPairs):
    '''
    write the data matrix to a file ready for Mr Bayes, which is similar to a nexus file
    :param fout1: the name of the output file
    :param dataMatrix: the datamatrix as a dict of dicts
    :param numLangs: the number of languages in this concept
    :param numSoundPairs: the number of sound pairs
    '''
    #open the file
    fout = codecs.open(fout1,"wb","utf-8")
    #write the first lines of the nexus file
    fout.write("#NEXUS"+"\n"+"\n")
    fout.write("BEGIN DATA;"+"\n"+"DIMENSIONS ntax="+str(numLangs)+" NCHAR="+str(numSoundPairs)+";\n"+"FORMAT DATATYPE=Restriction GAP=- MISSING=? interleave=yes;\n"+"MATRIX\n\n")
    for l in dataMatrix.index:
        fout.write(l.ljust(40)+"\t")
        fout.write(''.join(map(str,dataMatrix.ix[l].values))+'\n')

    fout.write("\n"+";\n"+"END;")
    fout.close()

def write_dataMatrix_withSynonyms(fout1, dataMatrix,numLangs, numSoundPairs):
    '''
    write the data matrix to a file ready for Mr Bayes, which is similar to a nexus file
    :param fout1: the name of the output file
    :param dataMatrix: the datamatrix as a dict of dicts
    :param numLangs: the number of languages in this concept
    :param numSoundPairs: the number of sound pairs
    '''
    #open the file
    fout = codecs.open(fout1,"wb","utf-8")
    #write the first lines of the nexus file
    fout.write("#NEXUS"+"\n"+"\n")
    fout.write("BEGIN DATA;"+"\n"+"DIMENSIONS ntax="+str(numLangs)+" NCHAR="+str(numSoundPairs)+";\n"+"FORMAT DATATYPE=Restriction GAP=- MISSING=? interleave=yes;\n"+"MATRIX\n\n")
    for l, vals in dataMatrix.items():
        fout.write(l.ljust(40)+"\t")
        fout.write(''.join(map(str,vals))+'\n')

    fout.write("\n"+";\n"+"END;")
    fout.close()


if __name__ == '__main__':
    

    ##read the tree to get the list of leaves for reading nelex data and sort it out correctly   
    list_leaves, dummy_tree = read_tree()
    print list_leaves
    ##get the data as dictionary with key = concept value= list of tuples (lang, asjp_word), get a dictionary with unique languages per concept key=concept val=list of langs
    nelex_dict, unique_langs_concept = read_nelex(list_leaves)

    ##get the pmi scores
    pmiScores = pd.read_csv('input/pmiScores.csv',index_col=0)
    ## get the gap penalties
    gp1,gp2 = pd.read_csv('input/gapPenalties.csv',header=None,squeeze=True,index_col=0).values[:2]
    ##get the sounds
    sounds = array(pmiScores.index)
    count = 0
    ##for each concept, align the data via tcoffee, save the matrix and reconstruct a concept tree
    for concept, data_list in nelex_dict.items():
 
         
        count += 1
        #if not concept in cur_trees :
        print concept

        ##read the tree and list of leaves again in case we need to prune the tree
        list_leaves, tree = read_tree()
           
        ##prune the tree if the number of leaves is not the same
        if len(list_leaves) != len(unique_langs_concept[concept]):
            tree.prune(unique_langs_concept[concept])
   
           
           
           
           
           
        ##########with synonyms###############
        pad, taxa,binMtx = compute_with_synonyms(data_list, tree, sounds, pmiScores, gp1, gp2)
        taxa = list(set(taxa))
        ##try to work with synonyms!!!!!!!!!!not perfect comment out if not wanted
        dataMtx = dict()
        for l in binMtx.index:
            ##if the language is not in the dictionary
            if not l in dataMtx:
                ##check the dimension of the array (for synonyms its 2)
                dimension = binMtx.ix[l].values.ndim
                ##if 1, take the array into the dict
                if dimension == 1:
                    dataMtx[l] = binMtx.ix[l].values
                      
                ##if not merge the arrays into one keeping the 1s
                elif dimension == 2:
                    array1 = binMtx.ix[l].values[0]
                    array2 = binMtx.ix[l].values[1]
                    for idx1, val in ndenumerate(array2):
                        if val == 1:
                            other_val = array1[idx1]
                            if other_val == 0:
                                array1[idx1] = val
                       
                    dataMtx[l] = array1
            ##if the language is already in the dict
            else:
                ##merge the 2 dimensional array
                synonym_array = binMtx.ix[l].values
                sym1 = synonym_array[0]
                sym2 = synonym_array[1]
                for idx2, val2 in ndenumerate(sym2):
                    if val2 == 1:
                        diff_val = sym1[idx2]
                        if diff_val == 0:
                            sym1[idx2] = val2
                  
                ##merge the first and the second entry in the dict to get one arry, keeping 1s
                np_array = dataMtx[l]
                for idx, val in ndenumerate(sym1):
                    if val == 1:
                        old_val = np_array[idx]
                        if old_val == 0:
                            np_array[idx] = val
                          
                dataMtx[l] = np_array

        ##create nexus files for further analysis with mrBayes
        ##output filename
        fout = "output/tcoffee_data_matrices_withSynonyms/"+concept+'_binMtx.phy'
        write_dataMatrix_withSynonyms(fout, dataMtx, str(len(taxa)), str(len(binMtx.columns)))
           
        ##work with tempfile to create iqTrees
        tempdir = "output/temp"
        with open(tempdir+"/binMtx.phy","w") as f:
            f.write(str(len(taxa))+' '+str(len(binMtx.columns))+'\n')
            for l, vals in dataMtx.items():
                f.write(l.ljust(pad))
                f.write(''.join(map(str,vals))+'\n')
        ##create a tree with iqtree
        '''
        Wenn es das naechste Mal wieder nicht tut, Matrizen speichern und in einem extra Programm in iqtree einlesen und ausfuehren
        -s indicates the input alignment file (binary matrix)
        -t specifies the starting tree (BioNJ tree)
        -bb ultrafast bootstrap with 1000 replicates -wbt to write all replicates and -wbtl replicates with branch lengths
        -nt defines the cores used for the multi cores variant
        '''
        p = subprocess.Popen('iqtree -s binMtx.phy -t BIONJ -bb 1000 -wbtl -safe -nt AUTO > /dev/null',shell=True,cwd=tempdir)
        os.waitpid(p.pid,0)
        #p = sdtout.close()
        
        ##files can be moved, so the tree does not be loaded and written again
        ##concept tree .treefile
        shutil.move(tempdir+'/binMtx.phy.treefile','output/iqTrees/'+concept+'_iqTree.tre')
        ##consensus tree .contree
        shutil.move(tempdir+'/binMtx.phy.contree','consensusTrees/'+concept+'_consensus.tre')
        ##bootstrap replicates .ufboot
        shutil.move(tempdir+'/binMtx.phy.ufboot','bootstrapReplicates/'+concept+'_bootstrapReplicates.ufboot')
        ##the consensus tree is in file .phy.consensus!!!!!!!!!!!!!!!!!!!!
        ##
        #conceptTree = Tree(tempdir+'/binMtx.phy.treefile')
        for filename in os.listdir(tempdir):
            os.remove(os.path.join(tempdir,filename))
           
        #conceptTree.write(format=1,outfile='output/iqTrees/'+concept+'_iqTree.tre')
      
        print count
     
    
    