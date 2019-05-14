from Bio import PDB

class PDBparser:
    def __init__(self, filename):
        parser = PDB.PDBParser(PERMISSIVE=1)
        self.struct =parser.get_structure("structure", filename)
        self.nr_models = len(self.struct)

    def GetModel(self, id):
        if id < self.nr_models:
            return self.struct[id]
        else:
            print("Index out of range. There are " + str(self.nr_models) + " models.")

    #returns the chain with given id in given model
    def GetChainById(self, mod, id):
        return mod[id]

    def GetChain(self, mod, id):
        chain_list = PDB.Selection.unfold_entities(mod, "C")
        nr = len(chain_list)
        if id < nr:
            return chain_list[id]
        else:
            print("Index out of range. This model has " + str(nr) + " chains.")

    def GetResById(self, chain, id, het = " ", ins = " "):
        return chain[(het, id, ins)]

    def GetResByPos(self, chain, pos):
        nr = len(chain)
        if pos < nr:
            reslist = chain.get_list()
            return reslist[pos]
        else:
            print("Index out of range. This chain has " + str(nr) + " residues.")

# returns a list of non-amino acid residues
    def GetHetResidues(self, modnum):
        residues = self.GetModel(modnum).get_residues()
        hetlist = []
        for residue in residues:
            residue_id = residue.get_id()
            hetfield = residue_id[0]
            if hetfield[0] == "H":
                hetlist.append(residue)

    def GetAtomByPos(self, res, pos):
        nr = len(res)
        if pos < nr:
            return res[pos]
        else:
            print("Index out of range. This residue has " + str(nr) + " atoms.")

    def GetAtomById(self, res, id):
        return res[id]

    def GetNumOdChildren(self, ent):
        return len(ent)
    def GetNumOfMods(self):
        return  self.nr_models
    def GetNumOfChains(self):
        n = 0
        for mod in self.struct:
            n += len(mod)
        return n
    def GetNumOfResidues(self):
        n = 0
        for mod in self.struct:
            for chain in mod:
                n += len(chain)
        return n
    def GetNumOfAtoms(self):
        n = 0
        for mod in self.struct:
            for chain in mod:
                for res in chain:
                    n += len(res)
        return n

# returns the maximal width of the structure
    def GetWidth(self):
        max = 0
        atoms = self.struct.get_atoms()
        for atm in atoms:
            for atom in atoms:
                if atm - atom > max:
                    max = atm - atom
        return max

# returns all atoms (in a list) within given distance from a certain ligand (atom)
    def GetAtmsInRadius(self, hetatm, radius):
        list = []
        atoms = self.struct.get_atoms()
        for atom in atoms:
            if atom - hetatm <= radius:
                list.append(atom)
        return  list

# returns all residues (in a list) within given distance from a certain ligand (atom)
    def GetResInRadius(self, hetatm, radius):
        list = self.GetAtmsInRadius(hetatm, radius)
        reslist = PDB.Selection.unfold_entities(list, "R")
        return reslist

if __name__ == "__main__":
    parser = PDBparser("5y41.pdb")
    mod = parser.GetModel(0)
    chain = parser.GetChain(mod, 0)
    res = parser.GetResByPos(chain, 3)


