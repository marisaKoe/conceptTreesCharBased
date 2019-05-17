'''
Created on 17.08.2017

@author: marisakoe


'''
import re
from numpy import *

def nw(x,y,lodict,gp1,gp2):
    """
    Needleman-Wunsch algorithm for pairwise string alignment
    with affine gap penalties.
    'lodict' must be a dictionary with all symbol pairs as keys
    and match scores as values.
    gp1 and gp2 are gap penalties for opening/extending a gap.
    Returns the alignment score and one optimal alignment.
    """
    n,m = len(x),len(y)
    dp = zeros((n+1,m+1))
    pointers = zeros((n+1,m+1),int)
    for i in xrange(1,n+1):
        dp[i,0] = dp[i-1,0]+(gp2 if i>1 else gp1)
        pointers[i,0]=1
    for j in xrange(1,m+1):
        dp[0,j] = dp[0,j-1]+(gp2 if j>1 else gp1)
        pointers[0,j]=2
    for i in xrange(1,n+1):
        for j in xrange(1,m+1):
            match = dp[i-1,j-1]+lodict[x[i-1],y[j-1]]
            insert = dp[i-1,j]+(gp2 if pointers[i-1,j]==1 else gp1)
            delet = dp[i,j-1]+(gp2 if pointers[i,j-1]==2 else gp1)
            dp[i,j] = max([match,insert,delet])
            pointers[i,j] = argmax([match,insert,delet])
    alg = []
    i,j = n,m
    while(i>0 or j>0):
        pt = pointers[i,j]
        if pt==0:
            i-=1
            j-=1
            alg = [[x[i],y[j]]]+alg
        if pt==1:
            i-=1
            alg = [[x[i],'-']]+alg
        if pt==2:
            j-=1
            alg = [['-',y[j]]]+alg
    return dp[-1,-1],array([''.join(x) for x in array(alg).T])




def algnMtx(al, sounds):
    """
    Takes a pairwise alignment (i.e. a pair of gapped strings with identical length)
    as input and returns a matrix representation M as output.
    The matrix M is defined as M[i,j] = 1 if x[i] is matched with y[j]
    in the alignment, 0 else (where x,y are the two ungapped strings to be aligned).
    """
    w1 = ''.join(array([s for s in al[0] if s!='-']))
    w2 = ''.join(array([s for s in al[1] if s!='-']))
    dm = zeros((len(w1),len(w2)),int)
    i,j = 0,0
    for s1,s2 in array([list(w) for w in al]).T:
        if s1 in sounds:
            if s2 in sounds:
                dm[i,j] += 1
                i+=1
                j+=1
            else:
                i+=1
        else:
            j+=1
    return dm

def sHamming(x,y):
    """
    Takes two gapped strings and returns the hamming distance between
    them. Positions containing a gap in at least one string are ignored.
    """
    w1,w2 = array(list(x)),array(list(y))
    r= mean(w1[(w1!='-')*(w2!='-')]!=w2[(w1!='-')*(w2!='-')])
    if isnan(r): r=1.
    return r

def nwBlock(b1,b2,lib):
    """
    Needleman-Wunsch alignment of two aligned blocks b1 and b2,
    using the scores in the extended library lib.
    """
    def pos(gappedString,i):
        """
        Returns the index of gappedString[i] in the
        ungapped version thereof.
        if gappedString[i] is a gap, returns -1
        """
        if gappedString[i]!='-':
            return 0 if i==0 else sum(array(list(gappedString[:i]))!='-')
        else:
            return -1
    words1 = array([re.sub("-","",w) for w in b1])
    words2 = array([re.sub("-","",w) for w in b2])
    n,m = len(b1[0]),len(b2[0])
    dp = zeros((n+1,m+1))
    pointers = zeros((n+1,m+1),int)
    pointers[0,1:] = 2
    pointers[1:,0] = 1
    for i in xrange(1,n+1):
        for j in xrange(1,m+1):
            insert = dp[i-1,j]
            delet = dp[i,j-1]
            match = dp[i-1,j-1] + sum([0 if '-' in [gs1[i-1],gs2[j-1]] else lib[w1,w2][pos(gs1,i-1),pos(gs2,j-1)]
            for (w1,gs1) in zip(words1,b1)
            for (w2,gs2) in zip(words2,b2)])
            dp[i,j] = max([match,insert,delet])
            pointers[i,j] = argmax([match,insert,delet])
    al1 = transpose(array([list(w) for w in b1]))
    al2 = transpose(array([list(w) for w in b2]))
    alCombined = []
    while max(i,j)>0:
        p = pointers[i,j]
        if p == 0:
            alCombined = [list(al1[i-1]) + list(al2[j-1])] + alCombined
            i-=1
            j-=1
        elif p==1:
            alCombined = [list(al1[i-1])+['-']*len(b2)] + alCombined
            i-=1
        else: 
            alCombined = [['-']*len(b1)+list(al2[j-1])] + alCombined
            j-=1
    return array([''.join(x) for x in array(alCombined).T]),dp[-1,-1]


if __name__ == '__main__':
    pass