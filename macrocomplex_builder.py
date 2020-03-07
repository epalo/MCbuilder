# imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB
import argparse, os, sys, UserInteraction


#main function that is called when running the script
if __name__ == "__main__":
    fasta_files, pdb_files = UserInteraction.getUserInput()
    seq_record_list = []
    for seq in fasta_files:
        for seq_record in SeqIO.parse(seq, "fasta"):
            seq_record_list.append(seq_record)
    print (seq_record_list)
