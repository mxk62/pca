import sys
from rdkit import Chem
from rdkit.Chem import AllChem


class Reaction:
    """Represents a reaction."""

    def __init__(self, smiles):
        """Initialize entry."""
        self.smi = smiles.strip()
        self.react_smis, self.prod_smis = [s.split('.')
                                           for s in self.smi.split('>>')]
        self.is_published = None

    def get_descriptors(self):
        """Calculate and returns a list of reaction descriptors."""

        # Count carbons in reaction SMILES.
        carbons = self.smi.count('C') + self.smi.count('c')

        # Esitmate a total mass of reactants.
        mass = 0
        for smi in self.react_smis:
            m = Chem.MolFromSmiles(smi)
            if m is not None:
                for a in m.GetAtoms():
                    mass += a.GetMass()

        return [carbons, mass]

    def is_published(self):
        """Returns True if a reactions is published."""
        pass


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
