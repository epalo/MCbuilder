# imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB, pairwise2
import argparse, os, sys, UserInteraction


#main function that is called when running the script
if __name__ == "__main__":
    """ Macrocomplex builder based on structure superimposition."""
    
    fasta_files, pdb_files = UserInteraction.getUserInput()

    seq_record_list = []
    for seq in fasta_files:
        for seq_record in SeqIO.parse(seq, "fasta"):
            seq_record_list.append(seq_record)

    parser = PDB.PDBParser()
    interact_structure = []
    for pdb_struct in pdb_files:
        interact_structure.append(parser.get_structure(pdb_struct,pdb_struct))

    pdb_seq = []
    for pdb in pdb_files:
        for record in SeqIO.parse(pdb, "pdb-seqres"):
            pdb_seq.append(record)

    for i in range(len(pdb_seq)):
        for m in range(i):
            alignment = pairwise2.align.globalxx(pdb_seq[i].seq, pdb_seq[m].seq)[0]
            aln_seq_1 = alignment[0]
            aln_seq_2 = alignment[1]
            al_length = len(alignment[0])
            ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_1, aln_seq_2))
    #     if ident/al_length >= 0.95:  # If 95% of identity, return known sequence
    #         return k_seq
