# imports
import argparse
import os
import sys
import loggingSetup
import UserInteraction


def processInput():
    """ Read pdb and FASTA files. """

    log = loggingSetup.createLogger()
    log.info("Initialized")
    input_list = UserInteraction.getUserInput()
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
