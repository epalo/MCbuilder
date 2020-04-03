import random

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

    # returns a list of chains out of a list of chains that are similar to the input chain
    def get_homo_chains(self, list_of_chains):
        to_return = []
        for lst in list_of_chains:
            if self in lst:
                to_return = lst
                break
        return to_return

    def set_numeric_id(self, number_list):
        new_id = random.choice(number_list)
        self.__biopy_chain.id = new_id
        number_list.remove(new_id)
        return self.__biopy_chain, number_list
