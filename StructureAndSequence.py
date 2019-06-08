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
        self.hetlist = self.PDBparser.GetHetResidues(0)     # list of non-amino acid residues
        self.activesite = []
        # for testing!!
        #self.hetlist = PDB.Selection.unfold_entities(self.PDBparser.GetModel(0), "R")[23:28]
        if len(self.hetlist) > 0:
            print("There are " + str(len(self.hetlist)) + "non-aa residues.")
            resnum = input("Which residue is the ligand? ")
            res = self.hetlist[int(resnum)]
            atm = res.get_list()[0]
            self.activesite = self.PDBparser.GetResInRadius(atm, 2.5)
        self.seq = self.MSAparser.GetSeqById(sequence)      # the sequence corresponding to the structure

# changes the loaded structure and MSA for a new one
    def LoadNewStruct(self, msafile, structfile, seq):
        self.__init__(msafile, structfile, seq)


# assumes that residues are in the same order in the structure as in the alignment!!
# calculates the conservation scores of residues in the active site, returns them as a list
    def CalcConsInActiveS(self):
        jsd = JensenShannon.JSD(self.MSAparser)
        scores = []
        positions = self.RecalculateResPositions()
        for pos in positions:
            score = jsd.JSDivergence(pos, jsd.q)
            scores.append(score)
        return scores

# calculates conservation scores for all residues and compares their average to the average scores of the active site
# returns ratio of conservation in active site vs. in the whole molecule
    def GetRelCons(self, actscores):
        jsd = JensenShannon.JSD(self.MSAparser)
        #scores = []
        count = 0
        scores = 0
        for i in range(len(self.seq)):
            score = jsd.JSDivergence(i, jsd.q)
            scores += score
            count += 1
        avg = scores / count
        actscore = sum(actscores)
        #for num in actscores:
        #    actscore += num
        if len(actscores) > 0:
            avgact = actscore / len(actscores)
        else:
            avgact = 0
        print(avgact)
        print(avg)
        return avgact/ avg

# determines which residues (positions in the sequence) belong to the active site, stores them as a list
    def RecalculateResPositions(self):
        gaps = 0
        residues = PDB.Selection.unfold_entities(self.PDBparser.GetModel(0), "R")   # list of all residues in the model
        newactsite = sorted(self.activesite, key=lambda x: residues.index(x))
        self.activesite = newactsite
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
            if site_index < len(self.activesite)-1 and res == self.activesite[site_index + 1]:
                site_index += 1
                poslist.append(seqpos)
        return poslist


if __name__ == "__main__":
    tool = StructAndSeq("epomsa.txt", "1buy.pdb", 0)
    scores = tool.CalcConsInActiveS()
    rel = tool.GetRelCons(scores)
    print(rel)
