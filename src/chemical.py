import sys
from rdkit import Chem


class Chemical:
    """Represents a chemical compound."""

    def __init__(self, smiles):
        self.smiles = smiles.strip()
        try:
            self.mol = Chem.MolFromSmiles(self.smiles)
        except Exception:
            print 'Error: invalid compound SMILES: {}.'.format(self.smiles)
            sys.exit(1)

    def make_retrostep(self, transform):
        """Returns unique reactions obtained by retrosynthesis."""

        product_list = []
        try:
            product_list = transform.formula.RunReactants((self.mol,))
        except Exception:
            # For now, ignore silently any errors
            pass

        rxn_smis = set([])
        if product_list:
            for products in product_list:
                reactant_smis = '.'.join(Chem.MolToSmiles(p) for p in products)
            rxn_smis.add(reactant_smis + '>>' + self.smiles)
        return rxn_smis


if __name__ == '__main__':
    pass
