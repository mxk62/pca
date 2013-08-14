import sys
from rdkit import Chem
from rdkit.Chem import Descriptors


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

    def count_atoms(self):
        return len(self.mol.GetAtoms())

    def count_bonds(self):
        return len(self.mol.GetBonds())

    def count_rings(self):
        return self.mol.GetRingInfo().NumRings()

    def get_weight(self):
        return sum(a.GetMass() for a in self.mol.GetAtoms())

    def get_randic(self):
        kappa = 0
        for b in self.mol.GetBonds():
            di = b.GetBeginAtom().GetDegree()
            dj = b.GetEndAtom().GetDegree()
            kappa += 1.0 / (di * dj)
        return kappa

    def get_balaban(self):
        return Descriptors.BalabanJ(self.mol)

    def get_bertz(self):
        return Descriptors.BertzCT(self.mol)

    def get_ipc(self):
        return Descriptors.IPC(self.mol)

    def get_kier_flex(self):
        k1 = Descriptors.Kappa1(self.mol)
        k2 = Descriptors.Kappa2(self.mol)
        return k1 * k2 / len(self.mol.GetAtoms())

    def make_retrostep(self, transform):
        """Returns unique reaction smiles obtained by retrosynthesis."""

        rxn_smis = set([])

        product_list = []
        try:
            product_list = transform.formula.RunReactants((self.mol,))
        except Exception:
            return rxn_smis

        if product_list:
            for products in product_list:
                try:
                    [Chem.SanitizeMol(m) for m in products]
                except Exception:
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
