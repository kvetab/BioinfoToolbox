from Bio import SeqIO
import Hamming
import EditDist

class FASTAparser:
    def __init__(self, filename):
        self.records = list(SeqIO.parse(filename, "fasta"))

    def checkNum(self, n):
        if n > 0 and n <= len(self.records):
            return True
        else:
            print("The number is not valid.")
            return False

    #returns description of molecule number n, counting from 1
    def getDescr(self, n):
        if FASTAparser.checkNum(self, n):
            print(self.records[n-1].description)
            return self.records[n-1].description

    #returns sequence of molecule n
    def getSeq(self, n):
        if FASTAparser.checkNum(self, n):
            print(self.records[n-1].seq)
            return self.records[n - 1].seq

    #returns sequence length for molecule n
    def getSeqLen(self, n):
        if FASTAparser.checkNum(self, n):
            print(len(self.records[n-1].seq))
            return len(self.records[n-1].seq)

    #returns subsequence of molecule n, indexed from 1, including beginning and end inices
    def getSubseq(self, beg, end, n):
        if FASTAparser.checkNum(self, n) and beg <= end:
            print(self.records[n-1].seq[beg-1:end])
            return self.records[n-1].seq[beg-1:end]

def main():
    parser = FASTAparser("fasta_simple.txt")
    seq = parser.getSeq(1)
    seq2 = parser.getSeq(2)
    #descr = parser.getDescr(1)
    #len = parser.getSeqLen(1)
    #subs = parser.getSubseq(2,8,1)
    #hamming = Hamming.HammingDist()
    #dist = hamming.GetDist(seq, parser.getSeq(2))
    #edit = EditDist.EditDist(seq, seq2)
    #dist, matrix = edit.GetDist()
    #print(dist)
    #al1, al2 = edit.GetAlignment(matrix, 21, 21)
    #print(al1)
    #print(al2)
    #print(seq)



if __name__ == "__main__":
  main()