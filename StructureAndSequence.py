import PDBfile
import MSAparser
import JensenShannon
from Bio import PDB

# assumes one model in file
class StructAndSeq:
    # inputs: file with multiple sequence alignment, file with pdb structure, index of sequence in the file belonging to structure
    def __init__(self, msafile, structfile, sequence):
        self.PDBparser = PDBfile.PDBparser(structfile)
        self.MSAparser = MSAparser.MSA(msafile)
        #self.hetlist = self.PDBparser.GetHetResidues(0)     # list of non-amino acid residues
        # for testing!!
        self.hetlist = PDB.Selection.unfold_entities(self.PDBparser.GetModel(0), "R")[23:28]
        print(len(self.hetlist))
        resnum = input("Which residue is the ligand? ")
        res = self.hetlist[int(resnum)]
        # self.alnum = input("Which sequence in the alignment belongs to the structure?")
        atm = res.get_list()[0]
        self.activesite = self.PDBparser.GetResInRadius(atm, 2.5)
        # the sequence corresponding to the structure
        self.seq = self.MSAparser.GetSeqById(sequence)

    def LoadNewStruct(self, msafile, structfile, seq):
        self.__init__(msafile, structfile, seq)


# assumes that residues are in the same order in the structure as in the alignment!!
    def CalcConsInActiveS(self):
        jsd = JensenShannon.JSD(self.MSAparser)
        scores = []
        positions = self.RecalculateResPositions()
        for pos in positions:
            score = jsd.JSD(pos, jsd.q)
            scores.append(score)
        return scores

    def GetRelCons(self, actscores):
        jsd = JensenShannon.JSD(self.MSAparser)
        #scores = []
        count = 0
        scores = 0
        for i in range(len(self.seq)):
            score = jsd.JSD(i, jsd.q)
            scores += score
            count += 1
        avg = scores / count
        actscore = 0
        for num in actscores:
            actscore += num
        avgact = actscore / len(actscores)
        print(avgact)
        print(avg)
        return avgact/ avg


    def RecalculateResPositions(self):
        gaps = 0
        residues = PDB.Selection.unfold_entities(self.PDBparser.GetModel(0), "R")   # list of all resiues in the model
        site_index = -1
        resnum = -1
        seqpos = -1
        poslist = []
        for res in residues:
            resnum += 1
            seqpos += 1
            while self.seq[seqpos] == "-":
                seqpos += 1
                gaps += 1
            if res == self.activesite[site_index + 1]:
                site_index += 1
                poslist.append(seqpos)
        return poslist


if __name__ == "__main__":
    tool = StructAndSeq("epomsa.txt", "1buy.pdb", 0)
    scores = tool.CalcConsInActiveS()
    rel = tool.GetRelCons(scores)
    print(rel)
