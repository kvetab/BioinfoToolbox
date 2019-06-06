import numpy as np
import math

class EditDist:
    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.m = None

    def CalcMin(self, i, j, m):
        res1 = m[i-1, j, 0] + 1
        res2 = m[i, j-1, 0] + 1
        res3 = m[i-1, j-1, 0]
        if self.seq1[i-1] != self.seq2[j-1]:
            res3 += 1
        res = min(res1, res2, res3)
        if res == res1:
            parentx = i-1
            parenty = j
        elif res == res2:
            parentx = i
            parenty = j-1
        else:
            parentx = i-1
            parenty = j-1
        return res, parentx, parenty

# returns the edit distance of two sequences
    def GetDist(self):
        len1 = len(self.seq1)
        len2 = len(self.seq2)
        m = np.zeros((len1+1, len2+1, 3), dtype=int)
        for i in range(len1):
            m[i+1,0,0] = i
        for j in range(len2):
            m[0,j+1,0] = j
        for k in range(len1):
            i = k+1
            for l in range(len2):
                j = l + 1
                res, x, y = EditDist.CalcMin(self, i, j, m)
                m[i,j,0] = res
                m[i,j,1] = x
                m[i, j, 2] = y
        self.m = m
        return  m[len1, len2, 0]

    # returns one of the optimal alignments
    def GetAlignment(self):
        if self.m == None:
            self.GetDist()
        len1 = len(self.seq1)
        len2 = len(self.seq2)
        stack = []
        tup = (len1, len2)
        while tup != (0,0):
            stack.append(tup)
            tup = (self.m[tup[0], tup[1], 1], self.m[tup[0], tup[1], 2])
            # tup = m[tup[0], tup[1], 1]
        lasti = 0
        lastj = 0
        al1 = ""
        al2 = ""
        while stack != []:
            tup = stack.pop()
            i = tup[0]
            j = tup[1]
            if i - lasti > 0:
                al1 = al1 + self.seq1[i-1]
            else:
                al1 = al1 + "-"
            if j - lastj >0:
                al2 = al2 + self.seq2[j-1]
            else:
                al2 = al2 + "-"
            lastj = j
            lasti = i
        # print(al1)
        # print(al2)
        return al1, al2

if __name__ == "__main__":
    seq1 = "AACTGTCT"
    seq2 = "ATCTGGCCT"
    edit = EditDist(seq1, seq2)
    #edit.GetAlignment()
    print(edit.GetDist())
