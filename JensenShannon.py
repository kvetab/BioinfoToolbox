import math
import numpy as np
import MSAparser
import pandas as pd

AA = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P','S','T', 'W', 'Y', 'V']
aadict = pd.read_csv("dict.txt", sep=" ", header=0)
lamb = 0.5
q = np.full((20),1/20)


def RelativeEntropy(pc, q):
    res = 0
    for aa in AA:
        res += pc[aa]*math.log((pc[aa]/q[aa]), 2)
    return res

def JSD(pc,q):
    r = lamb*pc + (1-lamb)*q
    res = lamb * RelativeEntropy(pc, r) + (1 - lamb) * RelativeEntropy(q, r)
    return res

#returns weights for given position in MSA
def CalcColFreq(col, leng):
    aas = np.zeros(20)      #first stores number of occurences of aa, then weighted
    used = np.zeros(20)     #marks each aa used in given column to get total number of different aas
    seqs = np.zeros(leng)
    for aa in col:
        if aa != "-":
            used[aadict[aa]] = 1
            aas[aadict[aa]] += 1
    n = used.sum()
    aas = aas*n
    gap = True
    for i in range(leng):
        if col[i] != "-":
            gap = False
            residue = aadict.loc[0,col[i]]
            seqs[i] = 1/aas[residue]
    return seqs, gap


def CalcWeights(reader):
    #reader = MSAparser.MSA(msa)
    alignment = reader.GetAl()
    num = reader.num        #number of sequences in MSA
    weights = np.zeros(num)
    count = 0
    for i in range(reader.length):      #length of each sequence
        freq, gap = CalcColFreq(reader.GetColumn(i), num)
        weights += freq
        if not gap:
            count += 1
    weights = weights/count
    print(weights)
    return weights

def CalcDistr(parser, weights):
    pseudocount = 0.000001
    distributions = np.full((20, parser.length), pseudocount)
    k = -1
    for sequence in parser.alignment:
        k += 1
        weight = weights[k]
        for i in range(len(sequence)):
            if sequence[i] != "-":
                distributions[aadict.loc[0,sequence[i]], i] += weight
    distributions = distributions/(1+(20*0.000001))
    return distributions


if __name__ == "__main__":
    reader = MSAparser.MSA("MSAsimple.txt")
    weights = CalcWeights(reader)
    dist = CalcDistr(reader, weights)
    res = JSD(dist[:,51],q)
    print(res)

