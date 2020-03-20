# imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB, pairwise2
from Bio.PDB.Polypeptide import PPBuilder
import argparse, os, sys, UserInteraction
import random

class Chain(object):
    """ DESCRIPTION """

    def __init__(self, sequence, file_index, interactions):
        self.__sequence = sequence
        self.__file_index = file_index
        self.__interactions = interactions

    def get_sequence(self):
        return self.__sequence

    def get_file_index(self):
        return self.__file_index

    def get_interaction(self):
        return self.__interactions

#main function that is called when running the script
if __name__ == "__main__":
    """ Macrocomplex builder based on structure superimposition."""
    
    # obtaining fasta and pdb files 
    
    fasta_files, pdb_files = UserInteraction.getUserInput()

# PARSING OF DATA

# TODO: insert case of empty fasta file 
    seq_record_list = []
    for seq in fasta_files:
        for seq_record in SeqIO.parse(seq, "fasta"):
            seq_record_list.append(seq_record)

# TODO: insert case of empty pdb-file 
    parser = PDB.PDBParser()
    interact_structure = []
    for pdb_struct in pdb_files:
        interact_structure.append(parser.get_structure(pdb_struct,pdb_struct))
    print("interact structure", interact_structure)
    ppb = PPBuilder()
    pdb_seq = []
    for i in range(len(interact_structure)):
            peptide1 = ppb.build_peptides(interact_structure[i])[0]
            # saves the record as a chain object with pdb-file index and sequence
            peptide2 = ppb.build_peptides(interact_structure[i])[1]
            pdb_seq.append(Chain(peptide1.get_sequence(), i, peptide2))
            pdb_seq.append(Chain(peptide2.get_sequence(), i, peptide1))
    for i in range(len(pdb_seq)):
        print(pdb_seq[i])
    
    # for i in range(len(pdb_files)):
    #     for record in SeqIO.parse(pdb_files[i], "pdb-seqres"):
    #         # saves the record together with the index of the pdb file
    #         pdb_seq = pdb_seq.append([record,i])

# SEQUENCE ALIGNMENTS

# find the sequences that occur multiple times in pdb files and save all proteins for each structural aln in a separate list
sequences = []
for i in range(len(pdb_seq)):
    for m in range(i):
        print("i:", i)
        print("m:", m)
        # just check sequence alignments if sequences are not in the same pair
        if (pdb_seq[i].get_file_index() != pdb_seq[m].get_file_index()):
            print("Seq1:", pdb_seq[i].get_file_index())
            print("Seq2:", pdb_seq[m].get_file_index())
            # find the best alignment for two sequences (first element of list of alignments)
            alignment = pairwise2.align.globalxx(pdb_seq[i].get_sequence(), pdb_seq[m].get_sequence())[0]
            aln_seq_1 = alignment[0]
            aln_seq_2 = alignment[1]
            al_length = len(alignment[0])
            ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_1, aln_seq_2))
            if ident/al_length >= 0.95:
                print("is identical")
                inserted = True
                for similar_seq in sequences:
                    if pdb_seq[i] in similar_seq:
                        if pdb_seq[m] not in similar_seq:
                            similar_seq.append(pdb_seq[m])
                            inserted = False
                            break
                    if pdb_seq[m] in similar_seq:
                        if pdb_seq[i] not in similar_seq:
                            similar_seq.append(pdb_seq[i])
                            inserted = False
                            break
                    if pdb_seq[m] in similar_seq and pdb_seq[i] in similar_seq:
                        inserted = False
                        break
                if inserted:
                    sequences.append([pdb_seq[i], pdb_seq[m]])
                print(sequences)
for i in range(len(sequences)):
    for el in sequences[i]:
        print("elem", i, el.get_file_index())

# BUILDING UP THE COMPLEX

# start with pdb-file with the most interactions
chain_to_superimpose = random.choice(max(len(similar_seq)))
starting_complex = pdb_files[chain_to_superimpose.get_file_index]
#print()
create_macrocomplex(starting_complex, chain_to_superimpose)


def get_superimpose_options(chain_to_superimpose):
    for similar_seq in sequences:
        if chain_to_superimpose in similar_seq:
            # return the list of similar sequences without the chain itself??
            return similar_seq

# TODO: check both chains of starting complex and combine them to complete complex

final_complexes = []
def create_macrocomplex(current_complex, superimpose_chain, threshold):
    # reach threshold
    if (threshold == 0):
        # append to list of final complexes
        return final_complexes.append(current_complex) 
    superimpose_options = get_superimpose_options(superimpose_chain)
    # no other superimposition options
    if not superimpose_options:
        # append to list of final complexes
        return final_complexes.append(current_complex) 
    for chain_option in superimpose_options:
        # TODO: combine multiple superimpose options
        created_complex = superimpose(current_complex, chain_option)
        threshold = 0
        if created_complex.rms > threshold:
            # implement: checking if complex works (surface stuff)
            create_macrocomplex(created_complex, chain_option.get_interaction(),threshold-1)


def superimpose(chain_a, chain_b):
    superimp = PDB.Superimposer()
    superimp.set_atoms(chain_a, chain_b) 
    return superimp.apply(chain_b.get_atoms())

# Feasible complexes - No clashes when superimposed
# Backbone of the model is not interacting with the rest of the complex (threshold example 2 Amstrons)
# We have to obtain the atoms thanks to the PDB
# 1. Obtain the backnone atoms of the complex and the model
# 2. use Bio.PDB.NeighborSearch(complex_atoms) --> Generates a neigbour search tree to speed up distance calculations
# for atom in chain_atoms:
# clashes += bool(ns.search(atom.coord, 2))  # If this atom shows clashes, add 1 to the clashes counter
# 3. threshold of % clashes (to discard the chain) --> Pau 3%
