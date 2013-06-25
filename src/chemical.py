import sys
from rdkit import Chem


class Chemical:
    """Represents a chemical compound."""

    def __init__(self, smiles):
        self.smi = smiles.strip()
        try:
            self.mol = Chem.MolFromSmiles(self.smi)
        except Exception:
            print 'Error: invalid compound SMILES: {}.'.format(self.smi)
            sys.exit(1)

    def make_retrostep(self, transform):
        """Returns unique reactions obtained by retrosynthesis."""

        rxns = set([])
        product_list = []

        try:
            product_list = transform.formula.RunReactants((self.mol,))
        except Exception as error:
            # For now, ignore silently any errors
            #print 'Error: applying transform failed:', error
            pass

        if product_list:
            for products in product_list:
                reactant_smis = '.'.join(Chem.MolToSmiles(p) for p in products)
            rxns.add(reactant_smis + '>>' + self.smi)
        else:
            return None
        return rxns


if __name__ == '__main__':
    pass
