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
        self.functional_groups = None

    def make_retrostep(self, transform):
        """Returns unique reaction smiles obtained by retrosynthesis."""

        product_list = []
        try:
            product_list = transform.formula.RunReactants((self.mol,))
        except Exception:
            # For now, ignore silently any errors.
            pass

        rxn_smis = set([])
        if product_list:
            for products in product_list:
                if not products:
                    continue
                reactant_smis = '.'.join(Chem.MolToSmiles(p) for p in products)
            rxn_smis.add(reactant_smis + '>>' + self.smiles)
        return rxn_smis

    def find_groups(self, groups):
        """Finds which functional groups are present on the chemical.

        Argument 'groups' is expected to be a dictionary which keys are
        functional gropus SMARTS and values are RDKit's query mols representing
        a gienb group.
        """
        self.functional_groups = {}
        for smarts, pattern in groups.items():
            if self.mol.HasSubstructMatch(pattern):
                self.functional_groups.setdefault(smarts, 0)
                self.functional_groups[smarts] += \
                    len(self.mol.GetSubstructMatches(pattern))


if __name__ == '__main__':
    groups = {}
    with open('functional_groups.txt') as f:
        for line in f:
            gid, smarts = line.split()
            try:
                pattern = Chem.MolFromSmarts(smarts)
            except Exception:
                continue
            groups[smarts] = pattern

    chem = Chemical('NC(=O)c1ccc(C(N)=O)cc1')
    chem.find_groups(groups)
    print chem.smiles, chem.functional_groups
