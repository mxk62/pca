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
    
    def get_fraction_CSP3(self):
        return Descriptors.FractionCSP3(self.mol)
    
    def get_num_aliphatic_carbocycles(self):
        return Descriptors.NumAliphaticCarbocycles(self.mol)
    
    def get_num_aliphatic_heterocycles(self):
        return Descriptors.NumAliphaticHeterocycles(self.mol)
    
    def get_num_saturated_rings(self):
        return Descriptors.NumSaturatedRings(self.mol)
    
    def get_num_saturated_heterocycles(self):
        return Descriptors.NumSaturatedHeterocycles(self.mol)
    
    def get_num_saturated_carbocycles(self):
        return Descriptors.NumSaturatedCarbocycles(self.mol)
    
    def get_num_aromatic_rings(self):
        return Descriptors.NumAromaticRings(self.mol)
    
    def get_num_aromatic_heterocycles(self):
        return Descriptors.NumAromaticHeterocycles(self.mol)
    
    def get_num_aromatic_carbocycles(self):
        return Descriptors.NumAromaticCarbocycles(self.mol)
    
    def get_num_aliphatic_rings(self):
        return Descriptors.NumAliphaticRings(self.mol)
    
    def get_num_H_donors(self):
        return Descriptors.NumHDonors(self.mol)
    
    def get_num_H_acceptors(self):
        return Descriptors.NumHAcceptors(self.mol)
    
    def get_hetero_atoms_mass(self):
        mass = 0
        for a in self.mol.GetAtoms():
            if a.GetAtomicNum() != 6 and a.GetAtomicNum() != 1:
                mass += a.GetMass()
        return mass
    
    def get_fraction_hetero_atoms(self):
        total_atom = len(self.mol.GetAtoms())
        hetero_atom = Descriptors.NumHeteroatoms(self.mol)
        fraction = hetero_atom/total_atom
        return hetero_atom
    
    def get_fraction_rotatable_bonds(self):
        return Descriptors.NumRotatableBonds(self.mol)
    
    def get_num_valance_electrons(self):
        return Descriptors.NumValenceElectrons(self.mol)
    
    def get_num_radical_electrons(self):
        return Descriptors.NumRadicalElectrons(self.mol)
    
    def get_heavy_atom_mol_wt(self):
        return Descriptors.HeavyAtomMolWt(self.mol)

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
