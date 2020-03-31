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

#main function that is called when running the script
if __name__ == "__main__":
    """ Macrocomplex builder based on structure superimposition."""

    # USER INPUT
    # obtaining fasta and pdb files
    fasta_files, pdb_files, log = UserInteraction.process_input()

    # PARSING OF DATA

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
        ppb=PPBuilder()
        # build the peptide to obtain the sequences of the chains
        peptide = ppb.build_peptides(model)
        sequence_a = peptide[0].get_sequence()
        sequence_b = peptide[1].get_sequence()
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

    
    # SEQUENCE ALIGNMENTS 

    # find the sequences that occur multiple times in pdb files and save all proteins for each structural aln in a separate list
    clashes_dict = {}
    homo_chains = []
    for i in range(len(chains)):
        for m in range(i):
            # just check sequence alignments if homo_chains are not in the same pair
            if (chains[i].get_file_index() != chains[m].get_file_index()):
                # find the best alignment for two homo_chains (first element of list of alignments)
                alignment = pairwise2.align.globalxx(chains[i].get_sequence(), chains[m].get_sequence())[0]
                aln_seq_1 = alignment[0]
                aln_seq_2 = alignment[1]
                al_length = len(alignment[0])
                ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_1, aln_seq_2))
                if ident/al_length >= 0.95:
                    inserted = True
                    log.info(f"{chains[i].get_biopy_chain().get_id()} and {chains[m].get_biopy_chain().get_id()} have 95% or more sequence identity")
                    for similar_seq in homo_chains:
                        if chains[i] in similar_seq:
                            if chains[m] not in similar_seq:
                                similar_seq.append(chains[m])
                                inserted = False
                                break
                        if chains[m] in similar_seq:
                            if chains[i] not in similar_seq:
                                similar_seq.append(chains[i])
                                inserted = False
                                break
                        if chains[m] in similar_seq and chains[i] in similar_seq:
                            inserted = False
                            break
                    if inserted:
                        homo_chains.append([chains[i], chains[m]])

    # for elem in homo_chains:
    #     print("Next list:")
    #     for el in elem:
    #         print(el.get_biopy_chain().get_id())
    log.info(f"{len(homo_chains)} homologous chains found")
    
    # SETTING THE STOICHIOMETRY
    stoichiometry = UserInteraction.get_stoichiometry()
    print(stoichiometry)
    log.info(f"Stoichiometry has been set.{stoichiometry}")

    # the starting_model is the interaction with the most homologous chains to both chains
    starting_model= None
    interaction_sum = 0
    # find interaction with the most homo_chains
    for interaction in interactions:
        sum = len(interaction.get_chain_a().get_homo_chains(homo_chains)) + len(interaction.get_chain_b().get_homo_chains(homo_chains))
        if sum > interaction_sum :
            starting_interaction = interaction
            interaction_sum = sum
        # set starting complex with stoichometry
    starting_complex = Complex(starting_interaction.get_model(), [starting_interaction.get_chain_a(), starting_interaction.get_chain_b()], log)
    starting_complex.add_to_stoich(starting_interaction.get_chain_a(), homo_chains)
    starting_complex.add_to_stoich(starting_interaction.get_chain_b(), homo_chains)
    if stoichiometry:
        starting_complex.set_stoich_complex(stoichiometry)
    print("Start",starting_interaction.get_chain_a().get_biopy_chain().get_id())
    print("Start",starting_interaction.get_chain_b().get_biopy_chain().get_id())
    
    # BUILDING THE MODEL WITH THE STARTING COMPLEX
    best_complex = starting_complex.create_macrocomplex(homo_chains, 20, stoichiometry)
    UserInteraction.create_output_PDB(best_complex)



