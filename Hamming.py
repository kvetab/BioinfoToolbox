
class HammingDist:

    def Hamming(self, seq1, seq2):
        dist = 0
        str1 = str(seq1)
        str2 = str(seq2)
        for i in range(len(str1)):
            if str1[i] != str2[i]:
                dist += 1
        return dist
    def GetDist(self, seq1, seq2):
        if len(seq1) == len(seq2):
            return self.Hamming(seq1, seq2)
        else:
            print("The lengths of the sequences don't match.")

