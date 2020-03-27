# imports
import argparse
import re
# import os
# import sys

# def getUserInput():
""" Import user data.

Requires input of fasta file(.fasta or .fa) containg protein to be modelled along with PDB files (.pdb) containing relevant interactions.
Only the mentioned file types will be accepted.
"""
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

parser.add_argument('-s', '--stoichiometry',
                    dest="stoich",
                    action="store",
                    default=None,
                    help="Set stoichiometry for the macrocomplex. Input in standard form (e.g. A1B4C6)")

parser.add_argument('-r', '--random',
                    dest="random_output_num",
                    action="store",
                    default=None,
                    type=int,
                    help="Random how many sequence to print")

options = parser.parse_args()

def getUserInput():
    return options.infile

def getVerboseOption():
    return options.verbose

def getOutputDirectory():
    return options.outfile

def getStoichiometry():
    if options.stoich:
        regex = re.compile("([A-Z]+[0-9]+)", re.IGNORECASE)
        stoich = {}
        stoich_all = regex.findall(options.stoich)
        for chain in stoich_all:
            stoich[chain[0]] = int(chain[1:])
        return stoich
    else:
        return options.stoich


    # getting input files in fasta and pdb format
    # input_list = options.infile
    # fasta_files = []
    # pdb_files = []
    # if len(input_list) == 0:
    #     input_list = [os.getcwd()]
    # for input in input_list:
    #     if os.path.isdir(input):
    #         print(os.listdir(input))
    #         fasta_files = [os.path.join(input, f) for f in os.listdir(input) if f.endswith(".fa") or f.endswith(".fasta")]
    #         pdb_files = [os.path.join(input, f) for f in os.listdir(input) if f.endswith(".pdb")]
    #     elif os.path.isfile(input) and (input.endswith(".fa") or input.endswith(".fasta")):
    #         fasta_files.append(input)
    #     elif os.path.isfile(input) and (input.endswith(".pdb")):
    #         pdb_files.append(input)
    #
    #
    # if len(fasta_files) == 0 and len(pdb_files) == 0:
    #     raise Exception("No fasta or pdb files were found")
    # if len(fasta_files) == 0 or len(pdb_files) == 0:
    #     raise Exception("fasta or pdb file is missing")

    # return (fasta_files, pdb_files)
