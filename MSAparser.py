from Bio import AlignIO
import pandas as pd

class MSA:
    def __init__(self, filename):
        self.alignment = AlignIO.read(filename, "clustal")
        self.num = len(self.alignment)
        self.length = len(self.alignment[0])


    def GetSeqById(self, id):
        if id < self.num:
            return self.alignment[id]
        else:
            print("Index out of bounds. There are "+str(self.num)+" sequences.")

    def GetAl(self):
        return self.alignment
    def GetNumOfSeq(self):
        return self.num
    def GetLenOfAl(self):
        return self.length

    def GetColumn(self, id):
        if id < self.length:
            return self.alignment[:, id]
        else:
            print("Index out of bounds. There are " + str(self.length) + " columns.")

    def LoadMatrix(self, filename):
        self.matrix = pd.read_csv(filename, sep=";", header=0, index_col=0)

    def SumOfPairsCol(self, col):
        if col < len(self.alignment[0]):
            score = 0
            column = self.alignment[:,col]
            for i in range(0, len(column)):
                for j in range(i+1, len(column)):
                    score += self.matrix.loc[column[i], column[j]]
            #print(score)
            return score
        else:
            print("Column index out of range.")

    def SumOfPairsAl(self):
        score = 0
        for i in range(len(self.alignment[0])):
            score += self.SumOfPairsCol(i)
        print(score)
        return score

if __name__ == "__main__":
    msa = MSA("MSAsimple.txt")
    msa.LoadMatrix("BLOSUM62_.csv")
    sop = msa.SumOfPairsCol(50)
    print(sop)
    msa.SumOfPairsAl()