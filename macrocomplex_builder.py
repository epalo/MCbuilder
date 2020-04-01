# imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB, pairwise2
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
import argparse, os, sys, UserInteraction
import random , copy
from InteractingChain import InteractingChain
from Complex import Complex
from Interaction import Interaction
import string

#main function that is called when running the script
if __name__ == "__main__":
    """ Macrocomplex builder based on structure superimposition."""

    ######################################################
    #  USER INPUT
    # obtaining fasta and pdb files
    fasta_files, pdb_files, log = UserInteraction.process_input()
    protein_limit = UserInteraction.get_protein_limit()
    
    ######################################################
    #  PARSING DATA

    # TODO: insert case of empty fasta file
    # fasta files still needed????
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

    ######################################################
    # SEQUENCE ALIGNMENTS
    homo_chains = find_homologous_chains(chains)
	
    log.info(f"{len(homo_chains)} homologous chains found")

    ######################################################
    # SETTING THE STOICHIOMETRY
    stoichiometry = UserInteraction.get_stoichiometry()
    print(stoichiometry)
    log.info(f"Stoichiometry has been set.{stoichiometry}")

    ######################################################
    #SETTING THE STARTING COMPLEX

    #the starting_model is the interaction with the most homologous chains to both chains
    starting_model= None
    interaction_sum = 0
    for interaction in interactions:
        sum = len(interaction.get_chain_a().get_homo_chains(homo_chains)) + len(interaction.get_chain_b().get_homo_chains(homo_chains))
        if sum > interaction_sum :
            starting_interaction = interaction
            interaction_sum = sum
            
    starting_complex = Complex(starting_interaction.get_model(), [starting_interaction.get_chain_a(), starting_interaction.get_chain_b()], log)
    if stoichiometry:
        starting_complex.set_stoich_complex(stoichiometry)
        starting_complex.add_to_stoich(starting_interaction.get_chain_a(), homo_chains)
        starting_complex.add_to_stoich(starting_interaction.get_chain_b(), homo_chains)
    print("Start",starting_interaction.get_chain_a().get_biopy_chain().get_id())
    print("Start",starting_interaction.get_chain_b().get_biopy_chain().get_id())

    # BUILDING THE MODEL WITH THE STARTING COMPLEX
    # number_list = list(range(0,10000))
    # best_complex = starting_complex.create_macrocomplex(homo_chains,protein_limit, stoichiometry, number_list, homo_chains)
    # print("FINAL COMPLEX:",best_complex)
    # print("NUMBER OF CHAINS:", len(best_complex.get_chains()))
    # for chain in best_complex.get_model().get_chains():
    #     print("CHAIN IDS", chain.get_id())

    ######################################################
    # BUILDING THE COMPLEX

    number_list = list(range(0,10000))
    if UserInteraction.get_runtype_option():
        best_complex = starting_complex.create_macrocomplex_full(homo_chains,protein_limit, stoichiometry, number_list, homo_chains)
    else:
        best_complex = starting_complex.create_macrocomplex(homo_chains,protein_limit, stoichiometry, number_list, homo_chains)
    print("ready to print")
    print("best", best_complex.get_chains())
    print("starting", starting_complex.get_chains())
    new_id_list = list(string.ascii_letters)
    for chain in best_complex.get_model().get_chains():
        chain.id = random.choice(new_id_list)
        new_id_list.remove(chain.id)
    UserInteraction.create_output_PDB(best_complex)
    exit(1)

######################################################
# HELPER FUNCTIONS

protein = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', "UNK": "X"}
dna = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T'}
rna = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}

def all_residues_in_dict(chain, type_dict):
    for residue in chain.get_residues():
        if not residue.resname in type_dict:
            return False
    return True
    
def build_sequence(chain, type_dict):
    sequence = ""
    for residue in chain.get_residues():
        sequence = sequence + type_dict[residue.resname.strip(" ")]
    return sequence

def get_sequence_for_chain(chain):
    if all_residues_in_dict(chain, protein):
        # all residues are aminoacids
        ppb=PPBuilder()
        peptide = ppb.build_peptides(chain)
        return peptide[0].get_sequence()
    if all_residues_in_dict(chain, dna):
        # all residues are dna
        return build_sequence(chain, dna)
    if all_residues_in_dict(chain, rna):
        # all residues are rna
        return build_sequence(chain,rna)
    else:
        print("exception!")

# function gets a list of interacting chains and returns a list of lists with homologous chains
def find_homologous_chains(chains):
    homologous_chains = []
    for i in range(len(chains)):
        for j in range(i):
            # just do sequence alignments if chains are not in the same interacting pair
            if (chains[i].get_file_index() != chains[j].get_file_index()):
                # find the best alignment for two homo_chains (first element of list of alignments)
                alignment = pairwise2.align.globalxx(chains[i].get_sequence(), chains[j].get_sequence())[0]
                aln_seq_1 = alignment[0]
                aln_seq_2 = alignment[1]
                #al_length = len(alignment[0])
                #ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_1, aln_seq_2))
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
