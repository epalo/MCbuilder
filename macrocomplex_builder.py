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
import random , copy
#import chimera
#from DetectClash import detectClash

# maybe transform into an extend of the actual pdb chain class and add function to retrieve sequence
class Interacting_Chain():
    """ DESCRIPTION """

    def __init__(self, biopy_chain, file_index, sequence, interacting_chain=None):
        self.__biopy_chain = biopy_chain
        self.__file_index = file_index
        self.__sequence = sequence
        self.__interacting_chain = interacting_chain

    def get_biopy_chain(self):
        return self.__biopy_chain

    def get_file_index(self):
        return self.__file_index

    def get_sequence(self):
        return self.__sequence

    def get_interacting_chain(self):
        return self.__interacting_chain

    def set_biopy_chain(self, biopy_chain):
        self.__biopy_chain = biopy_chain

    def set_interacting_chain(self, interacting_chain):
        self.__interacting_chain = interacting_chain

    def __len__(self):
        return len(self.__sequence)

class Complex(object):
    """ DESCRIPTION """

# chain attribute needed?
    def __init__(self, model, chains, pdb_files=False, stoichiometry=None):
        self.__model = model
        self.__chains = chains
        self.__pdb_files = pdb_files
        self.__stoichiometry = stoichiometry

    def get_model(self):
        return self.__model

    def get_chains(self):
        return self.__chains

    def get_pdb_files(self):
        return self.__pdb_files

    def get_stoichiometry(self):
        return self.__stoichiometry

    def set_stoichiometry(self):
        self.__stoichiometry = stoichiometry

    def set_model(self, model):
        self.__model = model

    def set_chains(self,chains):
        self.__chains = chains

    def add_chain(self, chain):
        self.__model.add(chain.get_biopy_chain())
        self.__chains.append(chain)

    def calc_z_score(self):
        # how to calculate z_score?
        return

# an interaction is a model out of two chains
class Interaction():
    def __init__(self, model, chain_a, chain_b):
        self.__model = model
        self.__chain_a = chain_a
        self.__chain_b = chain_b

    def get_model(self):
        return self.__model

    def get_chain_a(self):
        return self.__chain_a

    def get_chain_b(self):
        return self.__chain_b

#TODO : Update UserInter stoichiometry if chain id number
#TODO : How to join homologous chains so that under one name

#main function that is called when running the script
if __name__ == "__main__":
    """ Macrocomplex builder based on structure superimposition."""

    # obtaining fasta and pdb files
    number_list = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
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
        model = parser.get_structure(pdb_files[i],pdb_files[i])[0]
        ppb=PPBuilder()
        # build the peptide to obtain the sequences of the chains
        peptide = ppb.build_peptides(model)
        sequence_a = peptide[0].get_sequence()
        sequence_b = peptide[1].get_sequence()
        # build up the list of chains for each interaction
        biopy_chain_a, biopy_chain_b = model.get_chains()
        interacting_a = Interacting_Chain(biopy_chain_a, i, sequence_a, biopy_chain_b)
        interacting_b = Interacting_Chain(biopy_chain_b, i, sequence_b,biopy_chain_a)
        interacting_a.set_interacting_chain(interacting_b)
        interacting_b.set_interacting_chain(interacting_a)
        interactions.append(Interaction(model, interacting_a, interacting_b))
        # if stoichiometry:
        #     stoich_complex.setdefault(biopy_chain_a.get_id(), 0)
        #     stoich_complex.setdefault(biopy_chain_b.get_id(), 0)
    # get all the chains of a list of interactions
    chains = []
    for interaction in interactions:
        chains.append(interaction.get_chain_a())
        chains.append(interaction.get_chain_b())
    log.info("PDB interactions processed")

    stoichiometry_temp = UserInteraction.getStoichiometry()
    stoichiometry = {}
    stoich_complex = {}
    if stoichiometry_temp:
        for stoich_chain in stoichiometry_temp:
            for elem in chains:
                if stoich_chain == elem.get_biopy_chain().get_id() and elem.get_biopy_chain() not in stoichiometry_temp:
                    print(stoich_chain)
                    print(elem.get_biopy_chain().get_full_id())

                    stoichiometry[elem.get_biopy_chain()] = stoichiometry_temp[stoich_chain]
        log.info(f"Stoichiometry has been set.{stoichiometry}")


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
    log.info(f"{len(homo_chains)} homologous chains found")


    # HELPER FUNCTIONS
    # returns a list of chains out of a list of chains that are similar to the input chain
    def get_homo_chains(chain):
        for lst in homo_chains:
            if chain in lst:
                log.info(f"Found homologs for chain {chain.get_biopy_chain().get_id()}")
                return lst
        else:
            log.info("No homologous chains found")
            return []

    # returns a list with all possible chains that can be added to a current complex
    def get_superimpose_options(current_complex):
        superimpose_options = []
        superimpose_options_verbose = ''
        for chain in current_complex.get_chains():
            similar_chains = get_homo_chains(chain)
            superimpose_options = superimpose_options + similar_chains
        for option in superimpose_options:
            superimpose_options_verbose = superimpose_options_verbose + " " + option.get_biopy_chain().get_id()
        log.info(f"The following chains are homologous to those currently in the complex:{superimpose_options_verbose}")
        # print("Superimpose options:", superimpose_options)
        return superimpose_options

    def create_stoich_complex(chain, stoich_complex):
        similar_chains = get_homo_chains(chain.get_chain_a())
        for chain in similar_chains:
            chain

    # TODO: check both chains of starting complex and combine them to complete complex
    def create_macrocomplex(current_complex, threshold):
        superimpose_options = get_superimpose_options(current_complex)
        best_complex = current_complex
        # print("best_complex",best_complex)
        # starting complex has no superimposition options
        if not superimpose_options:
            # then just return the starting complex
            log.info("There are no superimposition options available")
            # print("no options!")
            return best_complex
        else:
            for option in superimpose_options:
                log.info(f"Attempting to superimpose chain {option.get_biopy_chain().get_id()}")
                option_complex = superimpose(current_complex, option)
                if (option_complex == None):
                    log.warning("The current option could not be added!")
                else:
                    log.info("Option complex was be found!")
                    # print("Option Complex", option_complex)
                    # print("Current model:", option_complex.get_model())
                    # no other superimposition options for the complex available (leaf)
                    # or reached threshold
                    # or TODO: ADD STOICHOMETRY option

                    if not get_superimpose_options(option_complex) or \
                        (threshold == 0):
                        return best_complex
                        # if Z-Score for option complex is lower than for the current best complex replace it
                        if option_complex.calc_z_score < best_complex.calc_z_score:
                            best_complex = option_complex
                    # reach threshold
                    elif stoichiometry:
                        print(option_complex.get_stoichiometry())
                        if option_complex.get_stoichiometry()[option.get_interacting_chain().get_biopy_chain().get_id()] >= \
                                    stoichiometry[option.get_interacting_chain().get_biopy_chain().get_id()]:
                            log.info(f"Chain {option.get_interacting_chain().get_biopy_chain().get_id()} has reached its stoichiometry limit.")
                            # if Z-Score for option complex is lower than for the current best complex replace it
                            if option_complex.calc_z_score < best_complex.calc_z_score:
                                best_complex = option_complex
                    else:
                        # if we didn't reach the leaf yet, recursive call
                        print("recursion!")
                        if stoichiometry:
                            stoich_complex[option.get_interacting_chain().get_biopy_chain().get_id()]+= 1
                        create_macrocomplex(option_complex ,threshold-1)
        return best_complex

    def is_clashing(current_complex, atom_list):
        backbone = {"CA", "C1\'"}
        model_atoms = [atom for atom in current_complex.get_model().get_atoms() if atom.id in backbone]
        chain_atoms = [atom for atom in atom_list if atom.id in backbone]
        # for atom in model_atoms:
        #     print(atom.get_coord())
        clashes_list = []
        chain_list = []
        # for atom in chain_atoms:
        n_search = PDB.NeighborSearch(model_atoms) # Generates a neigbour search tree
        clashes = 0
        for atom in chain_atoms:
            #print(atom)
            clashes += bool(n_search.search(atom.coord, 1.7))  # If this atom shows clashes, add 1 to the clashes counter
            # clashes_list.append(n_search.search(atom.coord, 1.7))
        # print("Num of clashes:", clashes)
        if clashes/len(chain_atoms) >= 0.03:  # If more than 3% of atoms show clashes return yes
            # for item in clashes_list:
            #     chain_list.append(PDB.Selection.unfold_entities(item, 'C'))
            #
            # print(clashes_list)
            log.info(f"Leads to clashes! {chain_list}")
            # clashes_dict[chain] = [chain for chain in chain.get_coords()]
            # print("Clashes dict: ", clashes_dict)
            return True
        else:  # Otherwise return no
            return False

    # return all the chains of a current complex where a chain_b can possibly be superimposed
    def get_superimpose_positions(current_complex, chain_b):
        superimpose_positions = []
        homos_chain_b = get_homo_chains(chain_b)
        for chain in current_complex.get_chains():
            if chain in homos_chain_b:
                superimpose_positions.append(chain)
        return superimpose_positions

    def superimpose(current_complex, chain_b):
        # TODO: only with backbone
        # TODO: check all different positions ???

        # if no complex can be created with the requested chain it returns None
        created_complex = None
        superimposition_positions = get_superimpose_positions(current_complex, chain_b)
        superimp = PDB.Superimposer()
        best_superimposition_matrix = None
        best_rmsd = 10
        for chain in superimposition_positions:
            atoms_a = []
            atoms_b = []
            for elem in chain.get_biopy_chain().get_atoms():
                atoms_a.append(elem)
            # print("atoms a:",atoms_a)
            for elem in chain_b.get_biopy_chain().get_atoms():
                atoms_b.append(elem)
            # print("atoms b:",atoms_b)
            # setting fixed and moving atoms, calculate the superimposition matrix
            superimp.set_atoms(atoms_a, atoms_b)
            rmsd = superimp.rms
            # update the best superimposition according to its rmsd
            if rmsd < best_rmsd:
                # check if the superimposition leads to clashes
                log.info(f"Checking whether {chain_b.get_interacting_chain().get_biopy_chain().get_id()} has any clashes")
                chain_to_try = copy.copy(chain_b.get_interacting_chain())
                superimp.apply(chain_to_try.get_biopy_chain())
                chain_sorted = sorted(chain_to_try.get_biopy_chain().get_atoms())
                if not (is_clashing(current_complex, chain_sorted)):
                    log.info(f"Chain {chain_b.get_interacting_chain().get_biopy_chain().get_id()} did not have any clashes. Feasible addition.")
                    # print("feasible addition")
                    best_superimposition_matrix = superimp
                    best_rmsd = rmsd
                    best_chain_position = chain_to_try
                    print("best chain",best_chain_position)

        # backbone = {"CA", "C1\'"}
        # chain_atoms1 = [atom for atom in chain_b.get_biopy_chain().get_atoms() if atom.id in backbone]
        # for elem in chain_atoms1:
        #     print("atoms b:",elem.get_coord())

        # apply the superimposition matrix to chain_b and its interacting chain
        if not (best_superimposition_matrix == None):
            created_complex = current_complex
            # best_superimposition_matrix.apply(chain_b.get_interacting_chain().get_biopy_chain())

            # chain_atoms2 = [atom for atom in chain_b.get_biopy_chain().get_atoms() if atom.id in backbone]
            # for elem in chain_atoms2:
            #     print("atoms b:",elem.get_coord())
            try:
                created_complex.add_chain(best_chain_position)
            except PDB.PDBExceptions.PDBConstructionException:
                log.warning(f"ID twice error current id {best_chain_position.get_biopy_chain().get_id()}.")
                # chain_to_add = copy.copy(chain_b.get_interacting_chain())
                # log.warning(f"ID twice error current id {chain_to_add.get_biopy_chain().get_id()}.")
                best_chain_position.get_biopy_chain().id = random.choice(number_list)
                log.warning(f"ID twice error new id {best_chain_position.get_biopy_chain().get_id()}.")
                created_complex.add_chain(best_chain_position)
            # created_complex.add_chain(chain_b.get_interacting_chain())
            # print(created_complex)
        return created_complex

    # BUILDING UP THE COMPLEX

    # the starting_model is the interaction with the most homologous chains to both chains
    starting_model= None
    interaction_sum = 0
    # find interaction with the most homo_chains
    for interaction in interactions:
        if stoichiometry:
            create_stoich_complex(interaction, stoich_complex)
        sum = len(get_homo_chains(interaction.get_chain_a())) + len(get_homo_chains(interaction.get_chain_b()))
        if sum > interaction_sum :
            starting_interaction = interaction
            interaction_sum = sum
    if stoichiometry:
        log.info(f"This is the stoich complex {stoich_complex}")
        starting_complex = Complex(starting_interaction.get_model(), [starting_interaction.get_chain_a(), starting_interaction.get_chain_b()],stoich_complex)
    else:
        starting_complex = Complex(starting_interaction.get_model(), [starting_interaction.get_chain_a(), starting_interaction.get_chain_b()])
    best_complex = create_macrocomplex(starting_complex, 20)
    processInputFiles.createPDB(best_complex)
    # TODO: check for DNA
