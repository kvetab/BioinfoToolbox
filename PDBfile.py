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

    #returns the chain at position id in given model
    def GetChain(self, mod, id):
        nr = len(mod)
        if id < nr:
            return mod[id]
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

    def GetWidth(self):
        max = 0
        atoms = self.struct.get_atoms()
        for atm in atoms:
            for atom in atoms:
                if atm - atom > max:
                    max = atm - atom
        return max

    def GetAtmsInRadius(self, hetatm, radius):
        list = []
        atoms = self.struct.get_atoms()
        for atom in atoms:
            if atom - hetatm <= radius:
                list.append(atom)
        return  list

    def GetResInRadius(self, hetatm, radius):
        list = self.GetAtmsInRadius(self, hetatm, radius)
        reslist = PDB.Selection.unfold_entities(list, "R")
        return reslist


