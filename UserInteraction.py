# imports
import argparse
import os
import sys


def getUserInput():
    # flags used when running in terminal
    parser = argparse.ArgumentParser(description=" ")

    parser.add_argument('-i', '--input',
                        dest="infile",
                        action="store",
                        default=None,
                        nargs="*",
                        help="Input fasta and PDB files, or directory containing these")

    parser.add_argument('-o', '--output',
                        dest="outfile",
                        action="store",
                        default=None,
                        help=" ")

    parser.add_argument('-v', '--verbose',
                        dest="verbose",
                        action="store_true",
                        default=False,
                        help="Print progression log to standard error")

    parser.add_argument('-p', '--pattern',
                        dest="pattern",
                        action="store",
                        default=None,
                        help="Regular expression pattern to search in the translated sequences")

    parser.add_argument('-r', '--random',
                        dest="random_output_num",
                        action="store",
                        default=None,
                        type=int,
                        help="Random how many sequence to print")

    options = parser.parse_args()

    # getting input files in fasta and pdb format
    input_list = options.infile
    fasta_files = []
    pdb_files = []
    if len(input_list) == 0: 
        input_list = os.getcwd()
    for input in input_list:
        if os.path.isdir(input):
            fasta_files.append([os.path.join(input, f) for f in os.listdir(input) if f.endswith(".fa") or f.endswith(".fasta")])
            pdb_files.append([os.path.join(input, f) for f in os.listdir(input) if f.endswith(".pdb")])
        elif os.path.isfile(input) and (input.endswith(".fa") or input.endswith(".fasta")):
            fasta_files.append(input)
        elif os.path.isfile(input) and (input.endswith(".pdb")):
            pdb_files.append(input)
        

    if len(fasta_files) == 0 and len(pdb_files) == 0:
        raise Exception("No fasta or pdb files were found")
    if len(fasta_files) == 0 or len(pdb_files) == 0:
        raise Exception("fasta or pdb file is missing")
    
    return (fasta_files, pdb_files)
