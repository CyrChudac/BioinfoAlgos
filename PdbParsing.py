import Task
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.vectors import Vector
from Bio.PDB.Selection import unfold_entities

class PdbParser(Task.Task):
    def __init__(self):
        self.name = "pdb parsing"
    def help(self, pre: str, mult: int):
        print(pre*mult + "[file]")
        print(pre*(mult+1) + "'file' is a path to the pdb data file")
        print()
        print(pre*(mult+1) + "-s --show")
        print(pre*(mult+2) + "displays all residues and their atoms")
        print(pre*(mult+1) + "-a --atoms")
        print(pre*(mult+2) + "displays all atoms")
        print(pre*(mult+1) + "-r --residues")
        print(pre*(mult+2) + "displays all residues")
        print(pre*(mult+1) + "-q --quantities")
        print(pre*(mult+2) + "displays quantity of models, chains, residues and atoms")
        print(pre*(mult+1) + "-w --width")
        print(pre*(mult+2) + "displays the distance of the two far most atoms")
        print(pre*(mult+1) + "-la --ligand-atoms [ligand_name] [radius]")
        print(pre*(mult+2) + "displays all atoms in radius around ligand given by it's tri-letter name")
        print(pre*(mult+1) + "-lr --ligand-residues [ligand_name] [radius]")
        print(pre*(mult+2) + "displays all residues in radius around ligand given by it's tri-letter name")
        print(pre*(mult+1) + "-gl --get-ligands")
        print(pre*(mult+2) + "displays all ligands present in the structure")
    def averageVector(self, resi):
        result = Vector(0,0,0)
        for v in resi:
            result += v.get_vector()
        result /= len(resi)
        return result.get_array()
    def getClose(self, center_name: str, structure, radius: float, layer: str, predicate, toPrint):
        resi = None
        for m in structure:
            for c in m:
                for r in c:
                    if r.get_full_id()[3][0] == center_name:
                        resi = r
        if resi != None:
            ns = NeighborSearch(unfold_entities(structure, 'A'))
            center = self.averageVector(resi)
            atoms = ns.search(center, radius, layer)
            atoms = [i for i in atoms if predicate(i, resi)]
            return Task.Result(atoms, Task.Result.printList_advanced(toPrint))
        else:
            print("ligand of name " + center_name[2:] + " not found in given file")
            return None
    def run(self, params):
        if len(params) < 2:
            print("wrong number of parameters:")
            Task.Task.help(self)
        else:
            if not Task.isFile(params[0]):
                print("one of given paths is not a path")
                return None
            structure = PDBParser().get_structure("s", params[0])
            if params[1] == "-s" or params[1] == "--show":
                return Task.Result(structure, PdbParser.showStructure) 
            elif params[1] == "-a" or params[1] == "--atoms":
                r = ""
                for atom in structure.get_atoms():
                    r = r + atom.get_name() + " "
                return Task.Result(r)
            elif params[1] == "-r" or params[1] == "--residues":
                r = ""
                for model in structure:
                    for resi in model.get_residues():
                        r = r + resi.get_resname() + " "
                return Task.Result(r)
            elif params[1] == "-q" or params[1] == "--quantities":
                m = len(structure)
                c = 0
                r = 0
                a = 0
                for model in structure:
                    c += len(model)
                    for chain in model:
                        r += len(chain)
                        for resi in chain:
                            a += len(resi)
                return Task.Result((m,c,r,a), PdbParser.showQuantities)
            elif params[1] == "-w" or params[1] == "--width":
                m = 0
                for a1 in structure.get_atoms():
                    for a2 in structure.get_atoms():
                        m = max(m, a1 - a2)
                return Task.Result(m)
            elif params[1] == "-la" or params[1] == "--ligand-atoms":
                m = 0
                if len(params) < 4:
                    Task.Task.help(self)
                else:
                    name = "H_" + params[2]
                    return self.getClose(
                        name, structure,
                        float(params[3]),
                        'A',
                        lambda a,r: a.get_parent() != r,
                        lambda a: a.get_name() + "\t" + str(a.get_coord()))
            elif params[1] == "-lr" or params[1] == "--ligand-residues":
                m = 0
                if len(params) < 4:
                    Task.Task.help(self)
                else:
                    name = "H_" + params[2]
                    return self.getClose(
                        name, 
                        structure, 
                        float(params[3]),
                        'R',
                        lambda r1,r2: r1 != r2,
                        lambda r: r.get_resname() + "\t" + str(self.averageVector(r)))
            elif params[1] == "-gl" or params[1] == "--get-ligands":
                ligs = []
                for m in structure:
                    for c in m:
                        for r in c:
                            if r.get_full_id()[3][0][0:2] == "H_":
                                ligs.append(r)
                return Task.Result(ligs, Task.Result.printList_advanced(
                    lambda l: l.get_resname() + "\t" + str(self.averageVector(l))))
            else:
                print("argument " + params[1] + " is not a valid argument:")
                Task.Task.help(self)

    def showQuantities(quadr):
        m,c,r,a = quadr
        d = "\t"
        print("m" + d + "c" + d + "r" + d + "a")
        print(str(m) + d + str(c) + d + str(r) + d + str(a))
    def showStructure(structure):
        for model in structure:
            for chain in model:
                for resi in chain:
                    print(resi.get_resname())
                    r = "\t"
                    for atom in resi:
                        r = r + atom.get_name() + " "
                    print(r)
                print("--------------------------")