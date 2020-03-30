# imports
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO, PDB, pairwise2
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
import argparse, os, sys, user_interaction
import random , copy


# maybe transform into an extend of the actual pdb chain class and add function to retrieve sequence
class InteractingChain():
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

    def get_ca_atoms(self):
        backbone = {"CA", "C1\'"}
        model_atoms = [atom for atom in self.get_biopy_chain().get_atoms() if atom.id in backbone]
        return model_atoms

class Complex(object):
    """ DESCRIPTION """

    def __init__(self, model, chains, pdb_files=False, stoich_complex=None):
        self.__model = model
        self.__chains = chains
        self.__pdb_files = pdb_files
        self.__stoich_complex = stoich_complex

    def get_model(self):
        return self.__model

    def get_chains(self):
        return self.__chains

    def get_pdb_files(self):
        return self.__pdb_files

    def get_stoich_complex(self):
        return self.__stoich_complex

    def set_stoich_complex(self):
        self.__stoich_complex = stoich_complex

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

    # each entry in the stoichiometry is one representation for all homo-chains, so if a homologous chain occurs also the counter has to be set up 
    def add_to_stoich(self, chain):
        if not self.__stoich_complex == None:
            # get all the homo-chains for the chain to add
            homo_chain_ids = []
            for elem in get_homo_chains(chain):
                homo_chain_ids.append(elem.get_biopy_chain().get_id())
            # remove duplicates
            homo_chain_ids = list(dict.fromkeys(homo_chain_ids))
            print("homo chain ids:", homo_chain_ids)
            for id in homo_chain_ids:
                if id in self.__stoich_complex:
                    self.__stoich_complex[id] += 1
    
    def stoich_is_complete(self):
        is_complete = False
        if stoichiometry:
            # check if stoichiometry has reached its limit for all chains (the two dictionaries are equal)
            for id in stoichiometry:
                if not self.__stoich_complex[id] == stoichiometry[id]:
                    is_complete = False
                    break
                is_complete = True
        return is_complete
    
    def stoich_is_overfull(self):
        is_overfull = False
        if stoichiometry:
            for id in stoichiometry:
                if not self.__stoich_complex[id] <= stoichiometry[id]:
                    is_overfull = True
                    break 
                is_overfull = False
        return is_overfull 


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

#main function that is called when running the script
if __name__ == "__main__":
    """ Macrocomplex builder based on structure superimposition."""

    # USER INPUT
    # obtaining fasta and pdb files
    number_list = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    fasta_files, pdb_files, log = user_interaction.process_input()
    protein_limit = user_interaction.get_protein_limit()
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

    stoichiometry = user_interaction.get_stoichiometry()
    print(stoichiometry)
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

    # for elem in homo_chains:
    #     print("Next list:")
    #     for el in elem:
    #         print(el.get_biopy_chain().get_id())
    log.info(f"{len(homo_chains)} homologous chains found")

    # returns a list of chains out of a list of chains that are similar to the input chain
    def get_homo_chains(chain):
        to_return = []
        for lst in homo_chains:
            if chain in lst:
                log.info(f"Found homologs for chain {chain.get_biopy_chain().get_id()}")
                return lst


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
        return superimpose_options

    def create_stoich_complex(stoichiometry):
        stoich_complex = {}
        for id in stoichiometry:
            stoich_complex[id] = 0
        return stoich_complex


    def update_homo_chains(original, best_chain_position):
        i=0
        for list in homo_chains:
            chain_names = [chain.get_biopy_chain() for chain in list]
            if original.get_biopy_chain() in chain_names:
                break
            i+=1
        homo_chains[i].append(best_chain_position)

    def create_macrocomplex(current_complex):
        # superimpose_options = get_superimpose_options(current_complex)
        # # print("best_complex",best_complex)
        # # starting complex has no superimposition options
        # if not superimpose_options:
        #     # then just return the starting complex
        #     log.info("There are no superimposition options available")
        #     # print("no options!")
        #     return best_complex
        # else:
        best_complex = current_complex
        for option in get_superimpose_options(current_complex):
            print("option from complex", option)
            log.info(f"Attempting to superimpose chain {option.get_biopy_chain().get_id()}")
            print("stoich of current complex before superimposition:", current_complex.get_stoich_complex())
            option_complex = superimpose(current_complex, option)
            # don't go into recursion of there is no option-complex found 
            if (option_complex == None):
                log.warning("The current option could not be added!")
            else:
                log.info("Option complex was be found!")
                # no other superimposition options for the complex available (leaf)
                # or reached threshold
                # or reached stoichiometry
                # print("stoich of current complex after superimposition:", option_complex.get_stoich_complex())
                if len(option_complex.get_chains()) == protein_limit or \
                        option_complex.stoich_is_complete():
                    print("returning the final complex!")
                    user_interaction.create_output_PDB(best_complex)
                    return best_complex
                    break
                    # if Z-Score for option complex is lower than for the current best complex replace it
                    # if option_complex.calc_z_score < best_complex.calc_z_score:
                    #     best_complex = option_complex    
                else:
                    # if we didn't reach the leaf yet, recursive call
                    currently = [chain for chain in option_complex.get_model().get_chains()]
                    log.warning(f"Currently in complex: {currently}")
                    log.warning("recursion!")
                    print("recursion!")
                    create_macrocomplex(option_complex)
        return best_complex

    def is_clashing(current_complex, chain):
        backbone = {"CA", "C1\'"}
        model_atoms = [atom for atom in current_complex.get_model().get_atoms() if atom.id in backbone]
        chain_atoms = chain.get_ca_atoms()

        chain_list = []
        # for atom in chain_atoms:
        n_search = PDB.NeighborSearch(model_atoms) # Generates a neigbour search tree
        clashes = 0
        for atom in chain_atoms:
            clashes += bool(n_search.search(atom.coord, 1.7))  # If this atom shows clashes, add 1 to the clashes counter
        if clashes/len(chain_atoms) >= 0.03:  # If more than 3% of atoms show clashes return yes
            log.info(f"Leads to clashes! {chain_list}")
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


    def superimpose(current_complex, chain_to_superimp):
        # if no complex can be created with the requested chain it returns None
        created_complex = None
        superimposition_options = get_homo_chains(chain_to_superimp)
        superimp = PDB.Superimposer()
        best_chain_position = None
        best_rmsd = 10
        if superimposition_options == []:
            print("CHAIN THATS NOT WORKING:", chain_to_superimp.get_biopy_chain().get_id())
        for chain in superimposition_options:
            atoms_a = []
            atoms_b = []
            atoms_a = chain.get_ca_atoms()
            # for elem in chain.get_biopy_chain().get_atoms():
            #     atoms_a.append(elem)
            # print("atoms a:",atoms_a)
            atoms_b = chain_to_superimp.get_ca_atoms()
            # for elem in chain_to_superimp.get_biopy_chain().get_atoms():
            #     atoms_b.append(elem)
            # print("atoms b:",atoms_b)
            # setting fixed and moving atoms, calculate the superimposition matrix
            superimp.set_atoms(atoms_a, atoms_b)
            rmsd = superimp.rms
            # update the best superimposition according to its rmsd
            if rmsd < best_rmsd:
                # check if the superimposition leads to clashes
                log.info(f"Checking whether {chain.get_interacting_chain().get_biopy_chain().get_id()} has any clashes")
                chain_to_try = copy.copy(chain.get_interacting_chain())
                superimp.apply(chain_to_try.get_biopy_chain())
                if not (is_clashing(current_complex, chain_to_try)):
                    log.info(f"Chain {chain.get_interacting_chain().get_biopy_chain().get_id()} did not have any clashes. Feasible addition.")
                    original = chain.get_interacting_chain()
                    best_rmsd = rmsd
                    best_chain_position = chain_to_try


        # apply the superimposition matrix to chain_b and its interacting chain
        if not (best_chain_position == None):
            new_id = random.choice(number_list)
            best_chain_position.get_biopy_chain().id = new_id
            number_list.remove(new_id)
            update_homo_chains(original, best_chain_position)
            created_complex = copy.deepcopy(current_complex)

            try:
                created_complex.add_chain(best_chain_position)
                # print(created_complex)
                # print("Stoich before adding:",created_complex.get_stoich_complex())
                # print("checking for id:", best_chain_position.get_biopy_chain().get_id())
                # print("chain object:", best_chain_position)
                # if the added chain is specified in the stoichiometry change the counter of the added chain
                created_complex.add_to_stoich(best_chain_position)
                # print("Stoich after adding:",created_complex.get_stoich_complex())
                # print(best_chain_position.get_biopy_chain().get_id())
                # print(best_chain_position.get_interacting_chain().get_biopy_chain().get_id())
                # print(best_rmsd)

            except PDB.PDBExceptions.PDBConstructionException:
                created_complex.add_chain(best_chain_position)
            # created_complex.add_chain(chain_b.get_interacting_chain())
            # print(created_complex)
            # if stoichiometry limits are overfull set the option complex to None
            if created_complex.stoich_is_overfull():
                print("stoich is overfull!!!")
                created_complex = None
        print("returned complex:", created_complex)
        return created_complex

    # BUILDING UP THE COMPLEX

    # the starting_model is the interaction with the most homologous chains to both chains
    starting_model= None
    interaction_sum = 0
    # find interaction with the most homo_chains
    for interaction in interactions:
        sum = len(get_homo_chains(interaction.get_chain_a())) + len(get_homo_chains(interaction.get_chain_b()))
        if sum > interaction_sum :
            starting_interaction = interaction
            interaction_sum = sum
    if stoichiometry:
        # set starting complex with stoichometry
        stoich_complex = create_stoich_complex(stoichiometry)
        starting_complex = Complex(starting_interaction.get_model(), [starting_interaction.get_chain_a(), starting_interaction.get_chain_b()], stoich_complex=stoich_complex)
        starting_complex.add_to_stoich(starting_interaction.get_chain_a())
        starting_complex.add_to_stoich(starting_interaction.get_chain_b())
    else:
        # starting complex without stoichiometry
        print("set again")
        starting_complex = Complex(starting_interaction.get_model(), [starting_interaction.get_chain_a(), starting_interaction.get_chain_b()])
    # print("Start",starting_interaction.get_chain_a().get_biopy_chain().get_id())
    # print("Start",starting_interaction.get_chain_b().get_biopy_chain().get_id())
    best_complex = create_macrocomplex(starting_complex)
    print("ready to print")
    user_interaction.create_output_PDB(best_complex)
    exit(1)
    # TODO: check for DNA
