class Interaction():

    """ An Interaction Object stores a Bio.PDB.Model.Model, a first InteractingChain
    and a second Interacting Chain  """

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