# imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB, pairwise2
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB import PDBIO
from Bio.PDB.PDBIO import Select
import os

""" produces interaction files to create examples to test the program """

chains = []
pdb_file = "6om3.pdb"
parser = PDB.PDBParser()
structure = parser.get_structure("6om3",pdb_file)
for chain in structure.get_chains():
    chains.append(chain)
# chains.append(structure.get_chains())
print(chains)
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
        #print(clashes)
        if clashes >= 1:
            #print(clashes)
            all_interact.append((chains[i].get_id(), chains[m].get_id()))
        else:  # Otherwise return no
            continue# print("Doesn't")

flat_list = [item for sublist in all_interact for item in sublist]

for chain in chains:
    if chain.get_id() not in flat_list:
        print(f"{chain.get_id()} not in list")
io = PDBIO()
io.set_structure(structure)
cwd = os.getcwd()
path = os.path.join(cwd, structure.get_id())
os.mkdir(path)
for interact in all_interact:
    class ChainSelect(Select):
        def accept_chain(self, chain):
            if chain.get_id()==interact[0]:
                return 1
            elif chain.get_id()==interact[1]:
                return 1
            else:
                return 0
    io.save(path + '/' +interact[0]+interact[1]+'.pdb', ChainSelect())
