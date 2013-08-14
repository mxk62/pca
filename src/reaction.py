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

        # Calculate number of atoms.
        r_atoms = sum(chem.count_atoms() for chem in self.reactants)
        p_atoms = sum(chem.count_atoms() for chem in self.products)

        # Calculate number of bonds.
        r_bonds = sum(chem.count_bonds() for chem in self.reactants)
        p_bonds = sum(chem.count_bonds() for chem in self.products)

        # Calculate number of rings.
        r_rings = sum(chem.count_rings() for chem in self.reactants)
        p_rings = sum(chem.count_rings() for chem in self.products)

        # Calculate masses.
        r_mass = sum(chem.get_weight() for chem in self.reactants)
        p_mass = sum(chem.get_weight() for chem in self.products)

        # Calculate Randic indices.
        r_randic = sum(chem.get_randic() for chem in self.reactants)
        p_randic = sum(chem.get_randic() for chem in self.products)

        # Calculate Balban J indices.
        r_balaban = sum(chem.get_balaban() for chem in self.reactants)
        p_balaban = sum(chem.get_balaban() for chem in self.products)

        # Calculate Bert indices.
        r_bertz = sum(chem.get_bertz() for chem in self.reactants)
        p_bertz = sum(chem.get_bertz() for chem in self.products)

        # Calculate Kier flexibility indices.
        r_kier_flex = sum(chem.get_kier_flex() for chem in self.reactants)
        p_kier_flex = sum(chem.get_kier_flex() for chem in self.products)

        return [r_atoms, p_atoms, r_bonds, p_bonds, r_rings, p_rings,
                r_mass, p_mass, r_balaban, p_balaban, r_bertz, p_bertz,
                r_kier_flex, p_kier_flex, r_randic, p_randic]

    def get_group_descriptor(self, groups):
        """Return descriptor based on functional group count.

        Function returns a vector which elemnts indicates how many functional
        group of a given type are present in reaction's reactants.
        """
        group_count = {}
        for chem in self.reactants + self.products:
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
