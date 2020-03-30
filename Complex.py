import copy
from InteractingChain import InteractingChain

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
            print("homo chain ids:", homo_chain_ids)
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
        superimpose_options_verbose = ''
        for chain in self.__chains():
            similar_chains = chain.get_homo_chains(chain_list)
            superimpose_options = superimpose_options + similar_chains
        for option in superimpose_options:
            superimpose_options_verbose = superimpose_options_verbose + " " + option.get_biopy_chain().get_id()
        log.info(f"The following chains are homologous to those currently in the complex:{superimpose_options_verbose}")
        return superimpose_options

    def create_macrocomplex(self, chain_list, threshold):
        # superimpose_options = get_superimpose_options(current_complex)
        # # print("best_complex",best_complex)
        # # starting complex has no superimposition options
        # if not superimpose_options:
        #     # then just return the starting complex
        #     log.info("There are no superimposition options available")
        #     # print("no options!")
        #     return best_complex
        # else:
        best_complex = self
        for option in self.get_superimpose_options(chain_list):
            print("option from complex", option)
            log.info(f"Attempting to superimpose chain {option.get_biopy_chain().get_id()}")
            print("stoich of current complex before superimposition:", self.__stoich_complex())
            option_complex = self.superimpose(option)
            # don't go into recursion of there is no option-complex found 
            if (option_complex == None):
                log.warning("The current option could not be added!")
            else:
                log.info("Option complex was be found!")
                # no other superimposition options for the complex available (leaf)
                # or reached threshold
                # or reached stoichiometry 
                print("stoich of current complex after superimposition:", option_complex.get_stoich_complex())
                if not self.get_superimpose_options(chain_list) or \
                    (threshold == 0) or \
                        option_complex.stoich_is_complete():
                    print("returning the final complex!")
                    return best_complex
                    # if Z-Score for option complex is lower than for the current best complex replace it
                    # if option_complex.calc_z_score < best_complex.calc_z_score:
                    #     best_complex = option_complex    
                else:
                    # if we didn't reach the leaf yet, recursive call
                    currently = [chain for chain in option_complex.get_model().get_chains()]
                    log.warning(f"Currently in complex: {currently}")
                    log.warning("recursion!")
                    print("recursion!")
                    option_complex.create_macrocomplex(threshold-1)
        return best_complex


    def superimpose(self, chain_to_superimp):
        # if no complex can be created with the requested chain it returns None
        created_complex = None
        superimposition_options = chain_to_superimp.get_homo_chains()
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
                if not (self.is_clashing(chain_to_try)):
                    log.info(f"Chain {chain.get_interacting_chain().get_biopy_chain().get_id()} did not have any clashes. Feasible addition.")
                    original = chain.get_interacting_chain()
                    best_rmsd = rmsd
                    best_chain_position = chain_to_try


        # apply the superimposition matrix to chain_b and its interacting chain
        if not (best_chain_position == None):
            self.update_homo_chains(original, best_chain_position)
            created_complex = copy.deepcopy(self)

            try:
                created_complex.add_chain(best_chain_position)
                print(created_complex)
                print("Stoich before adding:",created_complex.get_stoich_complex())
                print("checking for id:", best_chain_position.get_biopy_chain().get_id())
                print("chain object:", best_chain_position)
                # if the added chain is specified in the stoichiometry change the counter of the added chain
                created_complex.add_to_stoich(best_chain_position)
                print("Stoich after adding:",created_complex.get_stoich_complex())
                print(best_chain_position.get_biopy_chain().get_id())
                print(best_chain_position.get_interacting_chain().get_biopy_chain().get_id())
                print(best_rmsd)

            except PDB.PDBExceptions.PDBConstructionException:
                number_list = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
                original = best_chain_position
                new_id = random.choice(number_list)
                best_chain_position.get_biopy_chain().id = new_id
                number_list.remove(new_id)
                created_complex.add_chain(best_chain_position)
            # created_complex.add_chain(chain_b.get_interacting_chain())
            # print(created_complex)
            # if stoichiometry limits are overfull set the option complex to None
            if created_complex.stoich_is_overfull():
                print("stoich is overfull!!!")
                created_complex = None
        print("returned complex:", created_complex)
        return created_complex
    
    def is_clashing(self, chain):
        backbone = {"CA", "C1\'"}
        model_atoms = [atom for atom in self.__model().get_atoms() if atom.id in backbone]
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
    def get_superimpose_positions(self, chain_b):
        superimpose_positions = []
        homos_chain_b = chain_b.get_homo_chains()
        for chain in self.__chains():
            if chain in homos_chain_b:
                superimpose_positions.append(chain)
        return superimpose_positions

    def update_homo_chains(self, original, best_chain_position, chain_list):
        i=0
        for list in chain_list:
            chain_names = [chain.get_biopy_chain() for chain in list]
            if original.get_biopy_chain() in chain_names:
                break
            i+=1
        chain_list[i].append(best_chain_position)

