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
