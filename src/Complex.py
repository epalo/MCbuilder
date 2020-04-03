from Bio import SeqIO, PDB, pairwise2
import copy, random
from InteractingChain import InteractingChain
import UserInteraction
import string

class Complex(object):

    """ A Complex object stores a Bio.PDB.Model.Model, a list of InteractingChains,
    a logging.getLogger(__name__) object, and (optionally) a dictionary for stoichiometry """

    def __init__(self, model, chains, logger, stoich_complex=None):
        self.__model = model
        self.__chains = chains
        self.__logger = logger
        self.__stoich_complex = stoich_complex

    def get_model(self):
        return self.__model

    def get_chains(self):
        return self.__chains

    def get_stoich_complex(self):
        return self.__stoich_complex

    def set_model(self, model):
        self.__model = model

    def set_chains(self,chains):
        self.__chains = chains

    def add_chain(self, chain):
        self.__model.add(chain.get_biopy_chain())
        self.__chains.append(chain)

    def set_stoich_complex(self, stoich):
        stoich_complex = {}
        for id in stoich:
          stoich_complex[id] = 0
        self.__stoich_complex = stoich_complex

    def calc_z_score(self):
        # how to calculate z_score?
        return

    # each entry in the stoichiometry is one representation for all homo-chains, so if a homologous chain occurs also the counter has to be set up
    def add_to_stoich(self, chain, chain_list):
        if not self.__stoich_complex == None:
            # get all the homo-chains for the chain to add
            homo_chain_ids = []
            for elem in chain.get_homo_chains(chain_list):
                homo_chain_ids.append(elem.get_biopy_chain().get_id())
            # remove duplicates
            homo_chain_ids = list(dict.fromkeys(homo_chain_ids))
            # print("homo chain ids:", homo_chain_ids)
            for id in homo_chain_ids:
                if id in self.__stoich_complex:
                    self.__stoich_complex[id] += 1

    def stoich_is_complete(self, stoich):
        is_complete = False
        if stoich:
            # check if stoichiometry has reached its limit for all chains (the two dictionaries are equal)
            for id in stoich:
                if not self.__stoich_complex[id] == stoich[id]:
                    is_complete = False
                    break
                is_complete = True
        return is_complete

    def stoich_is_overfull(self, stoich):
        is_overfull = False
        if stoich:
            for id in stoich:
                if not self.__stoich_complex[id] <= stoich[id]:
                    is_overfull = True
                    break
                is_overfull = False
        return is_overfull

    # returns a list with all possible chains that can be added to a current complex
    def get_superimpose_options(self, chain_list):
        superimpose_options = []
        for chain in self.__chains:
            similar_chains = chain.get_homo_chains(chain_list)
            superimpose_options = superimpose_options + similar_chains

        return superimpose_options

    def each_chain_occurs_in_list(self, chain_list):
        for chain in self.__chains:
            if not chain in chain_list:
                return False
        return True

    def get_most_interacting_chain(self, list_of_chains, homo_chain_list):
        best_chain = list_of_chains[0]
        most_homo_chains = len(best_chain.get_homo_chains(homo_chain_list))
        for chain in list_of_chains:
            num_homo_chains = len(chain.get_homo_chains(homo_chain_list))
            if (num_homo_chains > most_homo_chains):
                best_chain = chain
                most_homo_chains = num_homo_chains
        return best_chain
    

    def create_new_subunit(self,homo_chain_list, protein_limit, stoich, number_list, initial_chains, interaction_files, version):
        # get remaining interactions
        remaining_chains = [chain for chain,value in initial_chains.items() if value == 0]
        print("Remaining chain:",remaining_chains)
        # get next interaction with the most interactions
        new_start_chain = self.get_most_interacting_chain(remaining_chains, homo_chain_list)
        # TODO: set coordinates for the next subunit
        # add new_start_chain to optioncomplex and run again in recursive call (creation of new subunit)
        if new_start_chain in initial_chains:
            initial_chains[new_start_chain] = True
        number_list = new_start_chain.set_numeric_id(number_list)
        self.add_chain(new_start_chain)
        if(version == "full"):
            return self.create_macrocomplex_full(homo_chain_list, protein_limit, stoich, number_list, initial_chains, interaction_files)
        else:
            return self.create_macrocomplex(homo_chain_list, protein_limit, stoich, number_list, initial_chains, interaction_files)

    # simple version of algorithm:
    def create_macrocomplex(self, homo_chain_list, protein_limit, stoich, number_list, initial_chains, interaction_files):
        # adds each time the function is called all possibly addable chains AT ONCE
        for chain in self.__chains:
            option_complex, updated_numbers = self.superimpose(chain, homo_chain_list, stoich, number_list, initial_chains)
            if (option_complex == None):
                self.__logger.warning("The current option could not be added!")
            else:
                self.__logger.info("Option complex was be found!")
        # if end-options are reached, don't go into recursion
        # endoptions: - all chains are clashing, -protein-limit reached, -complete stoichiometry
        if  len(option_complex.get_chains()) == protein_limit or \
                option_complex.stoich_is_complete(stoich) or \
                len(option_complex.get_chains()) == len(self.__chains):
            # check if all pdb-files were used at least once
            if all(initial_chains.values()) or protein_limit or option_complex.stoich_is_complete(stoich):
                print("FINAL COMPLEX FOUND")
                print("Initial Values are",initial_chains.values())
                print("optionComplex is:",option_complex)
                return option_complex
            else: # if not all pdb-files were used at least once but further adding leads to clashes --> creation of new subunit
                return option_complex.create_new_subunit(homo_chain_list, protein_limit, stoich, number_list, initial_chains, interaction_files, "simple")
        else:
            # recursively add chains that can still be added to the current complex
            currently = [chain for chain in option_complex.get_model().get_chains()]
            self.__logger.warning(f"Currently in complex: {currently}")
            self.__logger.warning("recursion!")
            print("recursion!")
            return option_complex.create_macrocomplex(homo_chain_list, protein_limit, stoich, updated_numbers, initial_chains, interaction_files)

    # full version of algorithm:
    # in each recursive step one possible chain is added and from there the whole recursive tree is searched through for possible complexes
    def create_macrocomplex_full(self, homo_chain_list, protein_limit, stoich, number_list, initial_chains, interaction_files):
        for option in self.get_superimpose_options(homo_chain_list):
            # superimpose the option-chain to the current complex
            self.__logger.info(f"Attempting to superimpose chain {option.get_biopy_chain().get_id()}")
            option_complex, updated_numbers = self.superimpose(option, homo_chain_list, stoich, number_list, initial_chains)

            if (option_complex == None): # if no option complex was found don't go into the recursive call
                self.__logger.warning("The current option could not be added!")
            else:
                self.__logger.info("Option complex was be found!")
                # if end-options are reached, don't go into recursion
                # endoptions: - all chains are clashing, -protein-limit reached, -complete stoichiometry
                if  len(option_complex.get_chains()) == protein_limit or \
                        option_complex.stoich_is_complete(stoich) or \
                        len(option_complex.get_chains()) == len(self.__chains):
                    # check if all pdb-files were used at least once
                    if all(initial_chains.values()):
                        print(initial_chains.values())
                        print(option_complex)
                        print("COMPLEX FOUND")
                        # add option_complex to list of final complexes
                        return option_complex
                    else: # if not all pdb-files were used at least once but further adding leads to clashes --> creation of new subunit
                        return option_complex.create_new_subunit(homo_chain_list, protein_limit, stoich, number_list, initial_chains, interaction_files, "full")
                else:
                    # if we didn't reach the leaf yet, recursive call
                    currently = [chain for chain in option_complex.get_model().get_chains()]
                    self.__logger.warning(f"Currently in complex: {currently}")
                    self.__logger.warning("recursion!")
                    print("full recursion!")
                    return option_complex.create_macrocomplex_full(homo_chain_list, protein_limit, stoich, updated_numbers, initial_chains, interaction_files)



    def superimpose(self, chain_to_superimp, homo_chain_list, stoich, number_list, initial_chains):
        # if no complex can be created with the requested chain it returns None
        created_complex = None

        superimposition_options = [chain for chain in chain_to_superimp.get_homo_chains(homo_chain_list) if chain in initial_chains]

        original = None
        superimp = PDB.Superimposer()
        best_chain_position = None
        best_rmsd = UserInteraction.get_rmsd_threshold()

        for chain in superimposition_options:
            atoms_a = []
            atoms_b = []
            atoms_a = chain_to_superimp.get_ca_atoms()
            atoms_b = chain.get_ca_atoms()
            if len(atoms_a) > len(atoms_b):
                diff = len(atoms_a) - len(atoms_b)
                if diff/len(atoms_a) >= 0.1:
                    continue
                atoms_a = atoms_a[:-diff]
            elif len(atoms_b) > len(atoms_a):
                diff = len(atoms_b) - len(atoms_a)
                if diff/len(atoms_b) >= 0.1:
                    continue
                atoms_b = atoms_b[:-diff]
            # setting fixed and moving atoms, calculate the superimposition matrix
            superimp.set_atoms(atoms_a, atoms_b)
            rmsd = superimp.rms
            # update the best superimposition according to its rmsd
            if rmsd < best_rmsd:
                # check if the superimposition leads to clashes
                self.__logger.info(f"Checking whether {chain.get_interacting_chain().get_biopy_chain().get_id()} has any clashes")
                chain_to_try = copy.deepcopy(chain.get_interacting_chain())
                superimp.apply(chain_to_try.get_biopy_chain())
                if not (self.is_clashing(chain_to_try)):
                    self.__logger.info(f"Chain {chain.get_interacting_chain().get_biopy_chain().get_id()} did not have any clashes. Feasible addition.")
                    original = chain.get_interacting_chain()
                    best_rmsd = rmsd
                    best_chain_position = chain_to_try


        # apply the superimposition matrix to chain_b and its interacting chain
        if not (best_chain_position == None):
            #new_id_list = list(string.ascii_letters)
            new_id = random.choice(number_list)
            best_chain_position.get_biopy_chain().id = new_id
            number_list.remove(new_id)
            self.update_homo_chains(original, best_chain_position, homo_chain_list)
            created_complex = copy.deepcopy(self)

            created_complex.add_chain(best_chain_position)
            # set chain to True if its in the initial chain dictionary
            if original in initial_chains:
                initial_chains[original] = True
            print(initial_chains)
            # if the added chain is specified in the stoichiometry change the counter of the added chain
            created_complex.add_to_stoich(best_chain_position,homo_chain_list)


            # if stoichiometry limits are overfull set the option complex to None
            if created_complex.stoich_is_overfull(stoich):
                # print("stoich is overfull!!!")
                created_complex = None
        else:
            created_complex = self
        return created_complex, number_list

    def is_clashing(self, chain):
        backbone = {"CA", "C1\'"}
        model_atoms = [atom for atom in self.__model.get_atoms() if atom.id in backbone]
        chain_atoms = chain.get_ca_atoms()

        chain_list = []
        # for atom in chain_atoms:
        n_search = PDB.NeighborSearch(model_atoms) # Generates a neigbour search tree
        clashes = 0
        for atom in chain_atoms:
            clashes += bool(n_search.search(atom.coord, 1.9))  # If this atom shows clashes, add 1 to the clashes counter
        if clashes/len(chain_atoms) >= 0.03:  # If more than 3% of atoms show clashes return yes
            self.__logger.info(f"Leads to clashes! {chain_list}")
            return True
        else:  # Otherwise return no
            return False

    # return all the chains of a current complex where a chain_b can possibly be superimposed
    def get_superimpose_positions(self, chain_b):
        superimpose_positions = []
        homos_chain_b = chain_b.get_homo_chains()
        for chain in self.__chains():
            if chain in homos_chain_b:
                superimpose_positions.append(chain)
        return superimpose_positions

    def update_homo_chains(self, original, best_chain_position, chain_list):
        for list in chain_list:
            chain_names = [chain.get_biopy_chain() for chain in list]
            if original.get_biopy_chain() in chain_names:
                list.append(best_chain_position)
