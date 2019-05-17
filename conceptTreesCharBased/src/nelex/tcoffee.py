'''
Created on 17.08.2017

@author: marisakoe
'''

import pandas as pd
from ete3 import Tree
from numpy import *
from helper_methods import *

def tCoffee(wl, guideTree, lodict, gp1, gp2, sounds):
    """- wl is a pandas data frame containing the columns 'doculect' and 'ASJP'.
    - The ASJP-column must contain non-emtpy ASJP strings.
    - guideTree is a binary-branching rooted ete2-tree. It's leaf names must coincide
      exactly with the entries occurring in the 'doculect' column of wl.
    - Synonyms, i.e., multiple rows for the same doculect, are allowed.
    - lodict is a dictionary which gives a log-odds score for each pair of the 41 ASJP
      sound classes.
    - gp1 and gp2 are gap penalties, i.e., non-positive real numbers.
    The function returns a pandas data frame with the columns 'doculect' and 'alignment'.
    """
    tree = guideTree.copy()
    wlI = wl.copy()
    wlI['id'] = range(len(wlI))
    wlI.index = wlI.id.values

    for l in tree.get_leaf_names():
        ids = wlI[wlI.doculect==l].id.values
        for i in ids:
            (tree&l).add_child(name=i)
    ##wlI.ASJP.values = a list with all words (including the synonyms)
    lib = createExtendedLibrary(wlI.ASJP.values,lodict,gp1,gp2,sounds)
    wTree = Tree(tree.write(format=9))
    wTree.standardize()
    for nd in wTree.traverse('postorder'):
        if nd.is_leaf():
            nd.add_feature('algn',array([wlI.ASJP[int(nd.name)]]))
            nd.add_feature('nTaxa',array([nd.name]))
            #print nd.name
        else:
            dl,dr = nd.get_children()
            b1,b2 = dl.algn,dr.algn
            a = nwBlock(b1,b2,lib)[0]
            nd.add_feature('algn',a)
            nd.add_feature('nTaxa',concatenate([dl.nTaxa,dr.nTaxa]))
    algn = pd.DataFrame(wTree.algn,index=wTree.nTaxa,columns=['alignment'])
    algn.index = array(algn.index,int)
    algn['doculect'] = wlI.doculect[algn.index]
    algn = algn.ix[wlI.index][['doculect','alignment']]
    
    return algn




def createExtendedLibrary(words,lodict,gp1,gp2,sounds):
    """
    Takes a list of sequences and returns an extended library in the
    sense of the T-Coffee algorithm. An extended library is a dictionary with
    sequence pairs as keys and a score matrix as values.
    For a pair of sequences x,y and a corresponding score matrix M,
    M[i,j] is the score for aligning x[i] with y[j].
    """
    library = createLibrary(words,lodict,gp1,gp2,sounds)
    extLibrary = dict()
    for w1 in words:
        for w2 in words:
            dm = zeros((len(w1),len(w2)))
            for w3 in words:
                a1,s1 = library[w1,w3]
                a2,s2 = library[w3,w2]
                dm += (s1+s2)*dot(a1,a2)
            extLibrary[w1,w2] = dm
    #print extLibrary
    return extLibrary

def createLibrary(words,lodict,gp1,gp2,sounds):
    """
    Takes a list of sequences and returns a library in the sense of the
    T-Coffee algorithm. A library is a dictionary with sequence pairs
    as keys and pairwise alignments in matrix format as columns.
    """
    ##library key = tuple of a wordpair, value=tuple with an array with 1 for presence and 0 for absence and the hamming distance
    library=dict()
    for w1 in words:
        for w2 in words:
            if (w2,w1) in library:
                x = library[w2,w1]
                library[w1,w2] = x[0].T,x[1]
            else:
                a1,a2 = nw(w1,w2,lodict,gp1,gp2)[1]
                library[w1,w2] = algnMtx([a1,a2],sounds),(1-sHamming(a1,a2))
    #print library
    return library


if __name__ == '__main__':
    pass