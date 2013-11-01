from rdkit.Chem import AllChem
from chemical import Chemical
from popularity import POPULARITY


class Reaction:
    """Represents a reaction."""

    def __init__(self, smiles, popularity=None, tid=None, rxnid=None, year=None):
        """Initialize entry."""
        self.popularity = popularity
        self.tid = tid
        self.rxnid = rxnid
        self.year = year
        self.smiles = smiles.strip()
        self.react_smis, self.prod_smis = [s.split('.')
                                           for s in self.smiles.split('>>')]
        try:
            self.reactants = [Chemical(smi) for smi in self.react_smis]
            self.products = [Chemical(smi) for smi in self.prod_smis]
        except ValueError:
            raise ValueError('invalid substrates and/or products')

    def equals(self, other):
        """Return True if two reactions are considered to be equal.

        Reactions B is considered to be equal to reaction A, if its both
        reactant and product SMILES are respective subsets of reaction A's
        reactants and products SMILES.
        """
        return (set(other.react_smis).issubset(self.react_smis) and
                set(other.prod_smis).issubset(self.prod_smis))


    def get_descriptors(self):
        """Calculates and returns a list of reaction descriptors.

        Currently, we are using purely molecular descriptors. In an attempt
        use them to characterize a reaction, for each available descriptor,
        we are calculating a total sum its values for all reactants and
        products separately. The resulting pair of numbers is our reaction
        descriptor corresponding to a given molecular descriptor.

        Simple compounds with a single-atom H-depleted molecular graph (e.g.
        water, 'O') are *excluded* from the count as both  adjacency and
        distance matrices on which many topological descriptors rely are
        not well defined (they are just 0) in such a case.
        """

        # If any descriptors should be ignored, put their names in the list
        # below.
        blacklist = []

        results = []
        for name, descriptor in Chemical:
            if name in blacklist:
                continue
            results.append(sum(descriptor(chem) for chem in self.reactants
                               if chem.a.size > 1))
            results.append(sum(descriptor(chem) for chem in self.products
                               if chem.a.size > 1))
        return results

    def get_group_descriptor(self, groups):
        """Return descriptor based on functional group count.

        Function returns a vector which elements indicates how many functional
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

    def __init__(self, smarts, dbid=None):
        self.id = dbid
        self.smarts = smarts
        self.retrons, self.synthons = [s.split('.')
                                       for s in smarts.split('>>')]
        try:
            self.formula = AllChem.ReactionFromSmarts(smarts)
        except Exception:
            raise ValueError('invalid transform SMARTS')

        try:
            self.popularity = POPULARITY[dbid]
        except KeyError:
            self.popularity = -1


if __name__ == '__main__':
    # An example reaction with single atom reactant.
    smi = 'CCCCCCC(=C=O)CC(N)C(=O)O.O>>CCCCCC[C@H](C[C@H](N)C(=O)O)C(=O)O'
    rxn = Reaction(smi)
    rxn.get_descriptors()
