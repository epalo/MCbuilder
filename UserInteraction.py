# imports
import argparse
import re
import os
import sys
# import macrocomplex_builder
import logging
from Bio.PDB import PDBIO


# def getUserInput():
# """ Import user data.
#
# Requires input of fasta file(.fasta or .fa) containg protein to be modelled along with PDB files (.pdb) containing relevant interactions.
# Only the mentioned file types will be accepted.
# """
# flags used when running in terminal

parser = argparse.ArgumentParser(description="Macrocomplex builder, creates protein macrocomplex from individual PDB files and FASTA files.")

parser.add_argument('-i', '--input',
                    dest="infile",
                    action="store",
                    default=None,
                    nargs="*",
                    help="Input fasta and PDB files, or directory containing these")

parser.add_argument('-o', '--output',
                    dest="outfile",
                    action="store",
                    default="macrocomplex.pdb",
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

parser.add_argument('-l', '--limit',
                    dest="limit",
                    action="store",
                    default=None,
                    type=int,
                    help="Limit number of chains in the protein")

parser.add_argument('-c', '--complete',
                    dest="complete",
                    action="store_true",
                    default=False,
                    help="Run full recursion. Explores all possibilities")

parser.add_argument('-r', '--rmsd',
                    dest="rmsd_threshold",
                    action="store",
                    type=float,
                    default= 0.5,
                    help="Set up rmsd threshold. Default is set to 0.5")


options = parser.parse_args()


def get_userinput():
    return options.infile

def get_verbose_option():
    return options.verbose

def get_output_directory():
    return options.outfile

def get_runtype_option():
    return options.complete

def get_rmsd_threshold():
    return options.get_rmsd_threshold

def get_stoichiometry():
    if options.stoich:
        regex = re.compile("([A-Z]+[0-9]+)", re.IGNORECASE)
        stoich = {}
        stoich_all = regex.findall(options.stoich)
        for chain in stoich_all:
            stoich[chain[0]] = int(chain[1:])
        return stoich
    else:
        return options.stoich

def get_protein_limit():
    return options.limit



verbose =  get_verbose_option()
# create logger
def create_logger():
    log = logging.getLogger(__name__)
    log.setLevel(level=logging.INFO)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(message)s')

    if verbose:
            # create console handler for logger.
            soh = logging.StreamHandler()
            soh.setLevel(level=logging.INFO)
            soh.setFormatter(formatter)
    # create file handler for logger.
    fh = logging.FileHandler('mbuilder.log')
    fh.setLevel(level=logging.WARNING)
    fh.setFormatter(formatter)

    # add handlers to logger.
    if verbose:
        log.addHandler(soh)

    log.addHandler(fh)

    return log


def process_input():
    """ Read pdb and FASTA files. """
    log = create_logger()
    log.info("Initialized")
    input_list = get_userinput()
    fasta_files = []
    pdb_files = []
    if len(input_list) == 0:
        input_list = [os.getcwd()]
    for input in input_list:
        if os.path.isdir(input):
            print(os.listdir(input))
            fasta_files = [os.path.join(input, f) for f in os.listdir(input) if f.endswith(".fa") or f.endswith(".fasta")]
            pdb_files = [os.path.join(input, f) for f in os.listdir(input) if f.endswith(".pdb")]
        elif os.path.isfile(input) and (input.endswith(".fa") or input.endswith(".fasta")):
            fasta_files.append(input)
        elif os.path.isfile(input) and (input.endswith(".pdb")):
            pdb_files.append(input)


    if len(fasta_files) == 0 and len(pdb_files) == 0:
        raise Exception("No fasta or pdb files were found")
    if len(fasta_files) == 0 or len(pdb_files) == 0:
        raise Exception("fasta or pdb file is missing")

    log.info(f"{len(fasta_files)} FASTA file(s) and {len(pdb_files)} PDB file(s) were processed")

    return (fasta_files, pdb_files, log)

def create_output_PDB(best_complex):
    """ Create a PDB with the final complex """
    io = PDBIO()
    io.set_structure(best_complex.get_model())
    io.save(get_output_directory())
