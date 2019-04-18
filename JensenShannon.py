import math
import numpy as np
import MSAparser
import pandas as pd


class JSD:
    def __init__(self, reader):
        self.AA = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
        self.AAs = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P','S','T', 'W', 'Y', 'V']
        self.aadict = pd.read_csv("dict.txt", sep=" ", header=0)
        self.lamb = 0.5
        self.q = np.full((20),1/20)     # background AA distribution (here uniform)
        self.reader = reader

        self.CalcWeights()
        self.CalcDistr()


    def RelativeEntropy(self, pc, q):
        res = 0
        for aa in self.AA:
            res += pc[aa]*math.log((pc[aa]/q[aa]), 2)
        return res

    # calculates conservation score of given position
    # first argument is distribution of AAs in given column
    def JSD(self, pos, q):
        pc = self.dist[:,pos]
        r = self.lamb*pc + (1-self.lamb)*q
        res = self.lamb * self.RelativeEntropy(pc, r) + (1 - self.lamb) * self.RelativeEntropy(q, r)
        return res

    #returns weights for given position in MSA
    def CalcColFreq(self, col, leng):
        aas = np.zeros(20)      #first stores number of occurences of aa, then weighted
        used = np.zeros(20)     #marks each aa used in given column to get total number of different aas
        seqs = np.zeros(leng)
        for aa in col:
            if aa != "-":
                used[self.aadict[aa]] = 1
                aas[self.aadict[aa]] += 1
        n = used.sum()
        aas = aas*n
        gap = True
        for i in range(leng):
            if col[i] != "-":
                gap = False
                residue = self.aadict.loc[0,col[i]]
                seqs[i] = 1/aas[residue]
        self.seqs = seqs
        self.gap = gap
        return seqs, gap


    def CalcWeights(self):
        #reader = MSAparser.MSA(msa)
        alignment = reader.GetAl()
        num = reader.num        #number of sequences in MSA
        weightsk = np.zeros(num)
        count = 0
        for i in range(self.reader.length):      #length of each sequence
            freq, gap = self.CalcColFreq(reader.GetColumn(i), num)
            weightsk += freq
            if not gap:
                count += 1
        weightsk = weightsk/count
        self.weights = weightsk
        #print(weightsk)
        #return weights

# calculates distributions of AAs in each column
    def CalcDistr(self):
        pseudocount = 0.000001
        distributions = np.full((20, self.reader.length), pseudocount)
        k = -1
        for sequence in self.reader.alignment:
            k += 1
            weight = self.weights[k]
            for i in range(len(sequence)):
                if sequence[i] != "-":
                    distributions[self.aadict.loc[0,sequence[i]], i] += weight
        distributions = distributions/(1+(20*0.000001))
        self.dist = distributions
        return distributions


if __name__ == "__main__":
    reader = MSAparser.MSA("MSAsimple.txt")
    jsd = JSD(reader)
    #weights =
    #jsd.CalcWeights()
    #dist = jsd.CalcDistr()
    res = jsd.JSD(51, jsd.q)
    print(res)

