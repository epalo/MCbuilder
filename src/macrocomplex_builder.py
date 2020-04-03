from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB, pairwise2
from Bio.PDB.Polypeptide import PPBuilder, CaPPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
import argparse, os, sys, UserInteraction
import random , copy
from InteractingChain import InteractingChain
from Complex import *
from Interaction import Interaction
import string

###  Macrocomplex builder based on structural superimposition  ###

##################################################################

dna = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T'}
rna = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}

def all_residues_in_dict(chain, type_dict):
    for residue in chain.get_residues():
        if not residue.resname.strip() in type_dict:
            return False
    return True

def build_sequence(chain, type_dict):
    sequence = ""
    for residue in chain.get_residues():
        sequence = sequence + type_dict[residue.resname.strip(" ")]
    return sequence

def get_sequence_for_chain(chain):
    if all_residues_in_dict(chain, dna):
        # all residues are dna
        return build_sequence(chain, dna)
    elif all_residues_in_dict(chain, rna):
        # all residues are rna
        return build_sequence(chain,rna)
    else:
        # building the peptide with PPBuilder
        ppb=PPBuilder()
        peptide = ppb.build_peptides(chain)
        if peptide:
            return peptide[0].get_sequence()
        else:
            # trying to build the peptide with CaPPBuilder
            ca_ppb = CaPPBuilder()
            peptide = ca_ppb.build_peptides(chain)
            if peptide:
                return peptide[0].get_sequence()
            else:
                # peptide couldn't be build
                log.warning("This program does not support heteroatoms")
                return None

# function gets a list of interacting chains and returns a list of lists with homologous chains
def find_homologous_chains(chains):
    homologous_chains = []
    for i in range(len(chains)):
        for j in range(i):
            # just do sequence alignments if chains are not in the same interacting pair
            if (chains[i].get_file_index() != chains[j].get_file_index()):
                # find the best alignment for two homo_chains (first element of list of alignments)
                alignment = pairwise2.align.globalxx(chains[i].get_sequence(), chains[j].get_sequence())[0]
                if alignment[2]/alignment[4] >= 0.95:
                    inserted = True
                    log.info(f"{chains[i].get_biopy_chain().get_id()} and {chains[j].get_biopy_chain().get_id()} have 95% or more sequence identity")
                    for similar_seq in homologous_chains:
                        if chains[i] in similar_seq:
                            if chains[j] not in similar_seq:
                                similar_seq.append(chains[j])
                                inserted = False
                                break
                        if chains[j] in similar_seq:
                            if chains[i] not in similar_seq:
                                similar_seq.append(chains[i])
                                inserted = False
                                break
                        if chains[j] in similar_seq and chains[i] in similar_seq:
                            inserted = False
                            break
                    if inserted:
                        homologous_chains.append([chains[i], chains[j]])
    return homologous_chains

def get_most_interacting_interaction(interaction_list, homo_chain_list):
    best_interaction = None
    interaction_sum = 0
    for interaction in interaction_list:
        sum = len(interaction.get_chain_a().get_homo_chains(homo_chain_list)) + len(interaction.get_chain_b().get_homo_chains(homo_chain_list))
        if sum > interaction_sum :
            best_interaction = interaction
            interaction_sum = sum
    return best_interaction

if __name__ == "__main__":
    """ main function that is called when running the application """

   ##################################################################
    #  USER INPUT
    # obtaining fasta and pdb files
    fasta_files, pdb_files, log = UserInteraction.process_input()
    protein_limit = UserInteraction.get_protein_limit()

    ##################################################################
    #  PARSING DATA

    # TODO: insert case of empty fasta file
    seq_record_list = []
    for seq in fasta_files:
        for seq_record in SeqIO.parse(seq, "fasta"):
            seq_record_list.append(seq_record)

    # TODO: insert case of empty pdb-file
    parser = PDB.PDBParser()
    interactions = []
    chains = []


  # iterate through all pdb files and return a list of interaction objects
    for i in range(len(pdb_files)):
        model = parser.get_structure(pdb_files[i],pdb_files[i])[0]
        # get sequences
        chain_a, chain_b = model.get_chains()
        sequence_a = get_sequence_for_chain(chain_a)
        sequence_b = get_sequence_for_chain(chain_b)

        if sequence_a == None or sequence_b == None:
            continue
        # build up the list of chains for each interaction
        biopy_chain_a, biopy_chain_b = model.get_chains()
        interacting_a = InteractingChain(biopy_chain_a, i, sequence_a, biopy_chain_b)
        interacting_b = InteractingChain(biopy_chain_b, i, sequence_b,biopy_chain_a)
        interacting_a.set_interacting_chain(interacting_b)
        interacting_b.set_interacting_chain(interacting_a)
        interactions.append(Interaction(model, interacting_a, interacting_b))
        # create a list that contains all the interacting chains
        chains.append(interacting_a)
        chains.append(interacting_b)
        log.info("PDB interactions processed")

    ##################################################################
    # SEQUENCE ALIGNMENTS
    homo_chains = find_homologous_chains(chains)

    log.info(f"{len(homo_chains)} homologous chains found")

    ##################################################################
    # SETTING THE STOICHIOMETRY
    stoichiometry = UserInteraction.get_stoichiometry()
    log.info(f"Stoichiometry has been set.{stoichiometry}")

    ##################################################################
    #SETTING THE STARTING COMPLEX

    #the starting_model is the interaction with the most homologous chains to both chains
    starting_interaction= get_most_interacting_interaction(interactions, homo_chains)

    starting_complex = Complex(starting_interaction.get_model(), [starting_interaction.get_chain_a(), starting_interaction.get_chain_b()], log)
    if stoichiometry:
        starting_complex.set_stoich_complex(stoichiometry)
        starting_complex.add_to_stoich(starting_interaction.get_chain_a(), homo_chains)
        starting_complex.add_to_stoich(starting_interaction.get_chain_b(), homo_chains)

    ##################################################################
    # BUILDING THE COMPLEX

    number_list = list(range(0,10000))
    dict_of_chains = { i : 0 for i in chains }
    dict_of_chains[starting_interaction.get_chain_a()] = True
    dict_of_chains[starting_interaction.get_chain_b()] = True
    if UserInteraction.get_runtype_option():
        log.info("Creating Macrocomplex with full version.")
        starting_complex.create_macrocomplex_full(homo_chains,protein_limit, stoichiometry, number_list, dict_of_chains, interactions)
        macrocomplex_list = final_complexes
        # TODO: create output pdb's
        exit(1)
    else:
        log.info("Creating Macrocomplex with simple version.")
        best_complex = starting_complex.create_macrocomplex(homo_chains,protein_limit, stoichiometry, number_list, dict_of_chains, interactions)
        new_id_list = list(string.ascii_letters)
        for chain in best_complex.get_model().get_chains():
            chain.id = random.choice(new_id_list)
            new_id_list.remove(chain.id)
        UserInteraction.create_output_PDB(best_complex, log)
        exit(1)
