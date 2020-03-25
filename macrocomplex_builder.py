# imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB, pairwise2
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
import argparse, os, sys, UserInteraction
import processInputFiles
import random
import loggingSetup
#import chimera
#from DetectClash import detectClash

# maybe transform into an extend of the actual pdb chain class and add function to retrieve sequence 
class Interacting_Chain():
    """ DESCRIPTION """

    def __init__(self, chain, file_index, sequence, interacting_chain):
        self.__chain = chain 
        self.__file_index = file_index
        self.__structure = structure
        self.__sequence = sequence
        self.__interacting_chain = interacting_chain
    
    def get_chain(self):
        return self.__chain 

    def get_file_index(self):
        return self.__file_index
    
    def get_sequence(self):
        return self.__sequence

    def get_interacting_chain(self):
        return self.__interacting_chain

    def __len__(self):
        return len(self.__sequence)

class Complex(object):
    """ DESCRIPTION """

# chain attribute needed?
    def __init__(self, structure, chains, pdb_files=False):
        self.__structure = structure
        self.__chains = chains
        self.__pdb_files = pdb_files

    def get_structure(self):
        return self.__structure

    def get_chains(self):
        return self.__chains
    
    def get_pdb_files(self):
        return self.__pdb_files
    
    def add_chain(self, chain):
        return Complex(structure.get_model().add(chain.get_chain()),chains.append(chain))
    
    def calc_z_score(self):
        # how to calculate z_score? 
        return 

class Interaction():
    def __init__(self, structure, chain_a, chain_b):
        self.__structure = structure
        self.__chain_a = chain_a
        self.__chain_b = chain_b
    
    def get_structure(self):
        return self.__structure 

    def get_chain_a(self):
        return self.__chain_a
    
    def get_chain_b(self):
        return self.__chain_b
    

#main function that is called when running the script
if __name__ == "__main__":
    """ Macrocomplex builder based on structure superimposition."""

    # obtaining fasta and pdb files

    fasta_files, pdb_files, log = processInputFiles.processInput()


# PARSING OF DATA
# TODO: insert case of empty fasta file
    seq_record_list = []
    for seq in fasta_files:
        for seq_record in SeqIO.parse(seq, "fasta"):
            seq_record_list.append(seq_record)

# TODO: insert case of empty pdb-file
    parser = PDB.PDBParser()
    interactions = []
    # iterate through all pdb files and return a list of interaction objects
    for i in range(len(pdb_files)):
        structure = parser.get_structure(pdb_files[i],pdb_files[i])
        ppb=PPBuilder()
        # build the peptide to obtain the sequences of the chains
        peptide = ppb.build_peptides(structure)
        sequence_a = peptide[0].get_sequence()
        sequence_b = peptide[1].get_sequence()
        # build up the list of chains for each interaction 
        chain_a, chain_b = structure.get_chains()
        interacting_a = Interacting_Chain(chain_a, i, sequence_a, chain_b)
        interacting_b = Interacting_Chain(chain_b, i, sequence_b,chain_a)
        interactions.append(Interaction(structure, interacting_a, interacting_b))


    # get all the chains of a list of interactions
    chains = []
    for interaction in interactions:
        chains.append(interaction.get_chain_a())
        chains.append(interaction.get_chain_b())
    log.info("PDB interactions processed")

# SEQUENCE ALIGNMENTS

# find the sequences that occur multiple times in pdb files and save all proteins for each structural aln in a separate list
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
    log.info(f"{len(homo_chains)} homologous chains found")


    # HELPER FUNCTIONS
    # returns a list of chains out of a list of chains that are similar to the input chain
    def get_homo_chains(chain):
        for lst in homo_chains:
            if chain in lst:
                return lst
        else:
            return []

    # returns a list with all possible chains that can be added to a current complex
    def get_superimpose_options(current_complex):
        superimpose_options = []
        for chain in current_complex.get_chains():  
            similar_chains = get_homo_chains(chain)
            superimpose_options = superimpose_options + similar_chains
        print("Superimpose options:", superimpose_options)
        return superimpose_options
    # TODO: check both chains of starting complex and combine them to complete complex

    # final_complexes = []
    # def create_macrocomplex(current_complex, superimpose_chain, threshold):
    #     superimpose_options = get_superimpose_options(superimpose_chain)
    #     # no other superimposition options -> reached a leaf of the tree
    #     if not superimpose_options:
    #         # append to list of final complexes
    #         return superimpose(superimpose_chain,) 
    #      # reach threshold
    #     if (threshold == 0):
    #         # append to list of final complexes
    #         return 
    #     for chain_option in superimpose_options:
    #         # TODO: combine multiple superimpose options
    #         superimposition = superimpose(superimpose_chain, chain_option)
    #         created_complex = 
    #         current_rmsd = superimpose(superimpose_chain, chain_option)
    #         rmsd_threshold = 0
    #         if created_complex.rms < rmsd_threshold:
    #             if not is_clashing(current_complex, chain_option):
    #             #if not is_clashing2(created_complex):
    #                 create_macrocomplex(created_complex, chain_option.get_interaction(),threshold-1)


    # def create_macrocomplex(current_complex, next_chain, threshold):
    #     # gives us the new complex and the according rmsd
    #     superimposition = superimpose(superimpose_chain, next_chain)
        
    #     superimpose_options = get_superimpose_options(superimposition[0],chains)
    #     # new complex has no other superimposition options -> reached a leaf of the tree
    #     if not superimpose_options:
    #         # append to list of final complexes
    #         return superimposition 
    #      # reach threshold
    #     if (threshold == 0):
    #         # append to list of final complexes
    #         return superimposition
    #     ## add here case for stoichometry
    #     if FALSE:
    #         return superimposition
        
    #     for option in superimpose_options:
    #         created_complex = create_macrocomplex(created_complex, option ,threshold-1)
    #         current_structure = created_complex[]
    #         current_rmsd = created_complex[1]
                    


    def create_macrocomplex(current_complex, threshold):
        superimpose_options = get_superimpose_options(current_complex)
        best_complex = current_complex
        # starting complex has no superimposition options
        if not superimpose_options:
            # then just return the starting complex
            print("no options!")
            return best_complex
        else:
            for option in superimpose_options:
                option_complex = superimpose(current_complex, option)
                print("Option Complex", option_complex)
                print("Current Structure:", option_complex.get_structure())
                # no other superimposition options for the complex available (leaf)
                # or reached threshold
                # or TODO: ADD STOICHOMETRY option
                if not get_superimpose_options(option_complex) or \
                    (threshold == 0) or \
                        False:
                    # if Z-Score for option complex is lower than for the current best complex replace it
                    if option_complex.calc_z_score < best_complex.calc_z_score:
                        best_complex = option_complex
                # reach threshold
                else:
                    # if we didn't reach the leaf yet, recursive call
                    print("recursion!")
                    create_macrocomplex(option_complex ,threshold-1)
        return best_complex
            
    
    
    def is_clashing(current_complex, chain_to_superimpose):
        backbone = {"CA", "C1\'"}
        model_atoms = [atom for atom in current_complex.get_atoms() if atom.id in backbone]
        chain_atoms = [atom for atom in chain_to_superimpose if atom.id in backbone]
        n_search = PDB.NeighborSearch(model_atoms) # Generates a neigbour search tree
        clashes = 0
        for atom in chain_atoms:
            clashes += bool(n_search.search(atom.coord, 1.7))  # If this atom shows clashes, add 1 to the clashes counter
        if clashes/len(chain_atoms) >= 0.03:  # If more than 3% of atoms show clashes return yes
            return True
        else:  # Otherwise return no
            return False

    # def is_clashing2(created_complex):
    #     clashes = detectClash(created_complex.get_atoms())
    #     sum = 0
    #     for clashDict in clashes.values():
    #         sum += len(clashDict)
    #     numClashes = sum/2
    #     if numClashes/len() >= 0.03:  # If more than 3% of atoms show clashes return yes
    #         return True
    #     else:  # Otherwise return no
    #         return False

    # return all the chains of a current complex where a chain_b can possibly be superimposed
    def get_superimpose_positions(current_complex, chain_b):
        superimpose_positions = []
        print("Chain id:", chain_b.get_chain().get_id())
        homos_chain_b = get_homo_chains(chain_b)
        print("Homos chain b:",homos_chain_b)
        print("Current complex:",current_complex.get_chains())
        for chain in current_complex.get_chains():
            if chain in homos_chain_b:
                superimpose_positions.append(chain)
        return superimpose_positions

    def superimpose(current_complex, chain_b):
        # TODO: check for clashing 
        # TODO: only with backbone
        # TODO: check all different positions ???
        superimposition_positions = get_superimpose_positions(current_complex, chain_b)
        print("Superimpose positions",superimposition_positions)
        superimp = PDB.Superimposer()
        best_superimposition = None
        best_rmsd = 10
        for chain in superimposition_positions:
            atoms_a = []
            atoms_b = []
            for elem in chain.get_chain().get_atoms():
                atoms_a.append(elem)
            print("atoms a:",atoms_a)
            for elem in chain_b.get_chain().get_atoms():
                atoms_b.append(elem)
            print("atoms b:",atoms_b)
            # setting fixed and moving atoms 
            superimp.set_atoms(atoms_a, atoms_b)
            # apply the superimposition, TODO: think about copy??
            copy_chain_b = chain_b.get_chain()
            superimp.apply(copy_chain_b)
            rmsd = superimp.rms
            if rmsd < best_rmsd:
                best_superimposition = copy_chain_b
                best_rmsd = rmsd
        print("Best superimposition:",copy_chain_b)
        created_complex = current_complex.add_chain(best_superimposition)
        return created_complex

    # BUILDING UP THE COMPLEX

    starting_interaction = None
    interaction_sum = 0
    # find pdb-file with the most interactions
    for interaction in interactions:
        sum = len(get_homo_chains(interaction.get_chain_a())) + len(get_homo_chains(interaction.get_chain_b()))
        if sum > interaction_sum :
            starting_interaction = interaction
            interaction_sum = sum
    starting_complex = Complex(starting_interaction, [starting_interaction.get_chain_a(), starting_interaction.get_chain_b()])
    create_macrocomplex(starting_complex, 100)

    # TODO: check for DNA