# imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB, pairwise2
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
import argparse, os, sys
import random , copy
import string
from Bio.PDB import PDBIO
from Bio.PDB.PDBIO import Select


chains = []
pdb_file = "1G65.pdb"
parser = PDB.PDBParser()
structure = parser.get_structure(pdb_file[:-4],pdb_file)
for chain in structure.get_chains():
    chains.append(chain)
# chains.append(structure.get_chains())
backbone = {"CA", "C1\'"}
all_interact = []
for i in range(len(chains)):
    fixed_chain = []
    fixed_chain = [atom for atom in chains[i].get_atoms()]  # Gets only the backbone atoms
    # print(fixed_chain)
    ns = PDB.NeighborSearch(fixed_chain)  # Generates a neigbour search tree to speed up distance calculations
    for m in range(i):
        compare_chain = []
        compare_chain = [atom for atom in chains[m].get_atoms()]
        clashes = 0
        for atom in compare_chain:
            clashes += bool(ns.search(atom.coord, 3.5))  # If this atom shows clashes, add 1 to the clashes counter
        # print(clashes)
        if clashes/len(fixed_chain) >= 0.01 or clashes/len(fixed_chain) >= 0.01:
            all_interact.append((chains[i], chains[m]))
        else:  # Otherwise return no
            continue# print("Doesn't")

ppb=PPBuilder()

all_interact_copy = copy.deepcopy(all_interact)
for i in range(len(all_interact)):
    for m in range(i):

        same_AC = False
        same_AD = False
        same_BC = False
        same_BD = False

        peptideA = ppb.build_peptides(all_interact[i][0])[0]
        peptideB = ppb.build_peptides(all_interact[i][1])[0]
        peptideC = ppb.build_peptides(all_interact[m][0])[0]
        peptideD = ppb.build_peptides(all_interact[m][1])[0]

        alignmentAC = pairwise2.align.globalxx(peptideA.get_sequence(), peptideC.get_sequence())[0]
        aln_seq_A = alignmentAC[0]
        aln_seq_C = alignmentAC[1]
        al_length = len(alignmentAC[0])
        ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_A, aln_seq_C))
        if ident/al_length >= 0.95:
            same_AC = True

        alignmentAD = pairwise2.align.globalxx(peptideA.get_sequence(), peptideD.get_sequence())[0]
        aln_seq_A = alignmentAD[0]
        aln_seq_D = alignmentAD[1]
        al_length = len(alignmentAD[0])
        ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_A, aln_seq_D))
        if ident/al_length >= 0.95:
            same_AD= True

        alignmentBC = pairwise2.align.globalxx(peptideB.get_sequence(), peptideC.get_sequence())[0]
        aln_seq_B = alignmentBC[0]
        aln_seq_C = alignmentBC[1]
        al_length = len(alignmentBC[0])
        ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_B, aln_seq_C))
        if ident/al_length >= 0.95:
            same_BC = True

        alignmentBD = pairwise2.align.globalxx(peptideB.get_sequence(), peptideD.get_sequence())[0]
        aln_seq_B = alignmentBD[0]
        aln_seq_D = alignmentBD[1]
        al_length = len(alignmentBD[0])
        ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_B, aln_seq_D))
        if ident/al_length >= 0.95:
            same_BD = True

        if (same_AC and same_BD) or (same_AD and same_BC):
            if all_interact[m] in all_interact_copy:
                all_interact_copy.remove((all_interact[m]))

# flat_list = [item for sublist in all_interact for item in sublist]
print(all_interact_copy)

io = PDBIO()
io.set_structure(structure)
cwd = os.getcwd()
path = os.path.join(cwd, structure.get_id())
os.mkdir(path)
for interact in all_interact_copy:
    class ChainSelect(Select):
        def accept_chain(self, chain):
            if chain.get_id()==interact[0].get_id():
                return 1
            elif chain.get_id()==interact[1].get_id():
                return 1
            else:
                return 0
    io.save(path + '/' +interact[0].get_id()+interact[1].get_id()+'.pdb', ChainSelect())
