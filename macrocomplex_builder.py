# imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB
import argparse, os, sys, UserInteraction


#main function that is called when running the script
if __name__ == "__main__":
    fasta_files, pdb_files = UserInteraction.getUserInput()