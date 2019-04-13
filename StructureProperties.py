import PDBfile
from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt

class StructAnalysis:
    def __init__(self, model, RAD):
        self.hse = PDB.HSExposureCB(model, radius = RAD)
        self.hse = PDB.HSExposureCA(model, radius=RAD)
        self.hse = PDB.ExposureCN(model, radius=RAD)
        self.model = model
        self.calculated = False
        filename = "polarAA.txt"
        self.polarAAs = {}
        with open(filename) as fh:
            for line in fh:
                aa, polar = line.strip().split(' ')
                polar = int(polar)
                self.polarAAs[aa] = polar

# what is a good measure of if a residue is buried or exposed?
# HSE supposedly works, bu I could not find any threshold values.
# HSE beta Up is supposed to be the most informative according to Hamelryck (2005) -> used here
    def GetHSE(self, residue):
        CN = None
        bU = None
        try:
            # Contact number
            CN = residue.xtra["EXP_CN"]
            # HSE alpha up
            # aU = residue.xtra["EXP_HSE_A_U"]
            # HSE alpha down
            # aD = residue.xtra["EXP_HSE_A_D"]
            # HSE beta up
            bU = residue.xtra["EXP_HSE_B_U"]
            # HSE beta down
            # bD = residue.xtra["EXP_HSE_B_D"]
        except:
            pass
        return CN, bU

    def GetDiam(self, parser):
        return parser.GetWidth()

    def GetCNs(self):
        residues = self.model.get_residues()
        buried = 0
        exposed = 0
        other = 0
        self.cns = []
        self.bus = []
        for res in residues:
            cn, bu = self.GetHSE(res)
            if cn is not None:
                self.cns.append(cn)
                self.bus.append(bu)
            if cn is not None and bu < 13 and cn < 30:
                exposed += 1
            elif cn is not None and bu > 22 and cn > 30:
                buried += 1
            else:
                other += 1
        if buried > 0:
            ratio = exposed / buried
        else:
            ratio = exposed
        self.calculated = True
        return ratio

    def GetHistogram(self):
        if self.calculated == False:
            self.GetCNs()
        arrc = np.array(self.cns)
        arrb = np.array(self.bus)
        histc = np.histogram(arrc)
        plt.hist(arrc)
        plt.show()
        histb = np.histogram(arrb)
        plt.hist(arrb)
        plt.show()
        return histb, histc

# mode is an integer
# 0 for ratio of polar aa in the core and all aa in the core
# 1 is like 0 but for surface aa
# 2 for a ratio of polar aa in the core and polar aa on the surface
    def PortionPolar(self, mode):
        res = self.model.get_residues()
        totalcore = 0
        totalout = 0
        polarcore = 0
        polarout = 0
        for aa in res:
            if PDB.is_aa(aa):
                if aa.xtra["EXP_HSE_B_U"] <13:
                    polarout += self.polarAAs[aa.get_resname()]
                    totalout += 1
                elif aa.xtra["EXP_HSE_B_U"] > 13:
                    polarcore += self.polarAAs[aa.get_resname()]
                    totalcore += 1
        if mode == 0:
            return polarcore/totalcore
        elif mode == 1:
            return polarout/totalout
        else:
            return polarcore/polarout




"""
 Use the parser you have implemented, to process a PDB structure and given its structure:

    Compute the diameter of the protein. ++
    Compute the ratio of surface and buried amino acids. ++
    Output a histogram (or at least data for a histogram) of amino acids composition of buried and exposed amino acids. ++
    Quantify portion of polar amino acids in the core and on the surface of the protein.
    Use the structures of A2a and caffeine receptor and hemoglobin (1b0b) to quantify the portion surface and buried amino acids and also the portion of polar amino acids from previous task. Do you see any significant differences? If so, try to interpret them.

"""

parser = PDBfile.PDBparser("5y41.pdb")
mod = parser.GetModel(0)
structure = StructAnalysis(mod, 12)
print(structure.PortionPolar(2))
#ratio = structure.GetCNs()
#print(ratio)
#structure.GetHistogram()