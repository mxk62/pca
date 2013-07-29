import sys
from rdkit.Chem import AllChem
from chemical import Chemical


class Reaction:
    """Represents a reaction."""

    def __init__(self, smiles):
        """Initialize entry."""
        self.smiles = smiles.strip()
        self.react_smis, self.prod_smis = [s.split('.')
                                           for s in self.smiles.split('>>')]
        self.reactants = [Chemical(smi) for smi in self.react_smis]
        self.products = [Chemical(smi) for smi in self.prod_smis]
        self.is_published = None

    def get_descriptors(self):
        """Calculate and returns a list of reaction descriptors."""

        # Count carbons in reaction SMILES.
        carbons = self.smiles.count('C') + self.smiles.count('c')

        # Esitmate a total mass of reactants.
        mass = 0
        for chem in self.reactants:
            if chem.mol is not None:
                for a in chem.mol.GetAtoms():
                    mass += a.GetMass()

        return [carbons, mass]

    def get_group_descriptor(self, groups):
        """Return descriptor based on functional group count.

        Function returns a vector which elemnts indicates how many functional
        group of a given type are present in reaction's reactants.
        """

        group_count = {}
        for chem in self.reactants:
            for group, count in chem.functional_groups.items():
                group_count.setdefault(group, 0.0)
                group_count[group] += count

        descriptor = []
        for smarts in sorted(groups.keys()):
            descriptor.append(group_count.get(smarts, 0))
        return descriptor


class Transform:
    """Represents retrosynthetic transform."""

    def __init__(self, smarts):
        self.smarts = smarts
        self.retrons, self.synthons = [s.split('.')
                                       for s in smarts.split('>>')]
        try:
            self.formula = AllChem.ReactionFromSmarts(smarts)
        except Exception:
            print 'Error: invalid transform SMARTS: {}'.format(self.smarts)
            sys.exit(1)


if __name__ == '__main__':
    pass
