
import Task
from Bio import AlignIO
from Bio.PDB.Polypeptide import PPBuilder
import PdbParsing
import MsaParsing
import positionConservation

class MSAStructure(Task.Task):
    radius = 2.5
    def __init__(self):
        self.name = "msa structure"
    def help(self, pre: str, mult: int):
        print(pre*mult + "[alignment_file] [pdb_file] [ligand_name] [scoring_matrix] [threshold] [window_half] [radius]")
        print(pre*(mult+1) + "'alignment_file' is a path to the alignment file")
        print(pre*(mult+1) + "'pdb_file' is a path to the pdb data file")
        print(pre*(mult+1) + "'ligand_name' is the tri-letter name of the ligand")
        print(pre*(mult+1) + "'scoring_matrix' is a path to the scoring matrix file")
        print(pre*(mult+1) + "'threshold' is a number between 0 and 1 indicating")
        print(pre*(mult+2) + "how much is the sequence supposed to be conservated")
        print(pre*(mult+1) + "'window_half' indicates how many residues are supposed to be conservated")
        print(pre*(mult+2) + "to each side from the active residue")
        print(pre*(mult+1) + "'radius' radius in which is a residue considered to be interacting with the ligand")
        print(pre*(mult+2) + "it's optional parameter with default of " + str(MSAStructure.radius))
        print()
        print(pre*(mult+1) + "finds out, if the sequence in the active site is conservated regarding the alignment")
    def findSeq(self, align, seq):
        for i in range(len(align)):
            alignSeq = align[i].seq.replace("-","")
            if alignSeq == seq:
                return i
        return -1
    def findResi(self, seq, num):
        found = 0
        curr = 0
        while found < num and curr < len(seq):
            if seq[curr] != "-":
                found += 1
            curr += 1
        if found == num:
            return curr
        return -1
    def run(self, params):
        if len(params) < 6:
            Task.task.help(self)
        else:
            pdbParser = PdbParsing.PdbParser()
            structure = pdbParser.run([params[1],"-s"])
            if structure == None:
                print("your pdb file can't be read")
                return None
            structure = structure.val
            radius = str(MSAStructure.radius)
            if len(params) > 6:
                radius = params[6]
            around = pdbParser.run([params[1],"--ligand-residues", params[2], radius]).val
            if around == None:
                print("en error while finding the active site in the pdb file")
                return None
            if len(around) == 0:
                return Task.Result(around, MSAStructure.noActive)
            around = [r for r in around if r.get_full_id()[3][0] == " "]
            msaParser = MsaParsing.MSAParser()
            align = msaParser.run([params[0],"-s"])
            if align == None:
                print("your alignment file can't be read")
                return None
            align = align.val
            ppb=PPBuilder()
            present = False
            isConserv = True
            for resi in around:
                id = resi.get_full_id()
                seq = ppb.build_peptides(structure[id[1]][id[2]])
                resiNum = id[3][2]
                seqNum = self.findSeq(align,seq)
                if seqNum != -1:
                    inAlignNum = self.findResi(align[seqNum].seq, resiNum)
                    if inAlignNum >= 0:
                        present = True
                        pc = positionConservation.PositionConservation()
                        isConserv = isConserv and pc.run([params[0],params[3],"-ic","-t",params[4],"-w",params[5]])
            if present:
                return Task.Result(isConserv)
            else:
                print("no active residue was found in both alignment and sequence")
    def noActive(x):
        print("no residue in the active radius around the ligand")