import Task
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.vectors import Vector
from Bio.PDB.Selection import unfold_entities

class PDBParsing(Task.Task):
    class NotALigandError(Exception):
        """the chosen  ligand is not present"""
        pass
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
            for a in atoms:
                if predicate(a, resi):
                    print(toPrint(a))
            return atoms
        else:
            raise self.NotALigandError()
    def run(self, params):
        if len(params) <= 2:
            Task.Task.help(self)
        else:
            structure = PDBParser().get_structure("s", params[1])
            if params[2] == "-s" or params[2] == "--show":
                for model in structure:
                    for chain in model:
                        for resi in chain:
                            print(resi.get_resname())
                            r = "\t"
                            for atom in resi:
                                r = r + atom.get_name() + " "
                            print(r)
                        print("--------------------------")
                return structure
            elif params[2] == "-a" or params[2] == "--atoms":
                r = ""
                for atom in structure.get_atoms():
                    r = r + atom.get_name() + " "
                print(r)
                return r
            elif params[2] == "-r" or params[2] == "--residues":
                r = ""
                for model in structure:
                    for resi in model.get_residues():
                        r = r + resi.get_resname() + " "
                print(r)
                return r
            elif params[2] == "-q" or params[2] == "--quantities":
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
                d = "\t"
                print("m" + d + "c" + d + "r" + d + "a")
                print(str(m) + d + str(c) + d + str(r) + d + str(a))
                return {"models": m,
                        "chains": c,
                        "residues": r,
                        "atoms": a}
            elif params[2] == "-w" or params[2] == "--width":
                m = 0
                for a1 in structure.get_atoms():
                    for a2 in structure.get_atoms():
                        m = max(m, a1 - a2)
                print(m)
                return m
            elif params[2] == "-la" or params[2] == "--ligand-atoms":
                m = 0
                if len(params) < 4:
                    Task.Task.help(self)
                else:
                    name = "H_" + params[3]
                    return self.getClose(
                        name, structure,
                        int(params[4]),
                        'A',
                        lambda a,r: a.get_parent() != r,
                        lambda a: a.get_fullname() + "\t" + str(a.get_coord()))
            elif params[2] == "-lr" or params[2] == "--ligand-residues":
                m = 0
                if len(params) < 4:
                    Task.Task.help(self)
                else:
                    name = "H_" + params[3]
                    return self.getClose(
                        name, 
                        structure, 
                        int(params[4]), 'R',
                        lambda r1,r2: r1 != r2,
                        lambda r: r.get_resname() + "\t" + str(self.averageVector(r)))
            elif params[2] == "-gl" or params[2] == "--get-ligands":
                ligs = []
                for m in structure:
                    for c in m:
                        for r in c:
                            if r.get_full_id()[3][0][0:2] == "H_":
                                ligs.append(r)
                for l in ligs:
                    print(l.get_resname() + "\t" + str(self.averageVector(l)))
                return ligs
            else:
                print("argument " + params[2] + " is not a valid argument:")
                Task.Task.help(self)