from Bio import AlignIO
import numpy as np
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

    def GetColumn(self, id):
        if id < self.length:
            return self.alignment[:, id]
        else:
            print("Index out of bounds. There are " + str(self.length) + " columns.")

    def LoadMatrix(self, filename):
        self.matrix = np.loadtxt(filename,dtype=int)
        df = pd.read_csv('BLOSUM62.csv', sep=";", header=0, index_col=0)
