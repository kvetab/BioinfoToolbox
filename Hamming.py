import FASTAparser

class HammingDist:

    # calculates the Hamming distance of two sequences of equal length
    def Hamming(self, seq1, seq2):
        dist = 0
        str1 = str(seq1)
        str2 = str(seq2)
        for i in range(len(str1)):
            if str1[i] != str2[i]:
                dist += 1
        return dist

  # returns distance of two sequences or an error message, if their lengths differ.
    def GetDist(self, seq1, seq2):
        if len(seq1) == len(seq2):
            return self.Hamming(seq1, seq2)
        else:
            print("The lengths of the sequences don't match.")

if __name__ == "__main__":
    #seq1 = "AACTGTCA"
    #seq2 = "ATCTTTC"
    parser = FASTAparser.FASTAparser("fasta_simple.txt")
    seq1 = parser.getSeq(1)
    seq2 = parser.getSeq(2)
    dist = HammingDist()
    print(dist.GetDist(seq1, seq2))