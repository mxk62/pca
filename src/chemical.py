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

    def count_H_acceptors(self):
        return Descriptors.NumHAcceptors(self.mol)

    def count_H_donors(self):
        return Descriptors.NumHDonors(self.mol)

    def count_atoms(self):
        return len(self.mol.GetAtoms())

    def count_hetero_atoms(self):
        return Descriptors.NumHeteroatoms(self.mol)

    def count_bonds(self):
        return len(self.mol.GetBonds())

    def get_fraction_rotatable_bonds(self):
        return Descriptors.NumRotatableBonds(self.mol)

    def count_rings(self):
        return self.mol.GetRingInfo().NumRings()

    def count_aliphatic_rings(self):
        return Descriptors.NumAliphaticRings(self.mol)

    def count_aromatic_rings(self):
        return Descriptors.NumAromaticRings(self.mol)

    def count_saturated_rings(self):
        return Descriptors.NumSaturatedRings(self.mol)

    def count_aliphatic_carbocycles(self):
        return Descriptors.NumAliphaticCarbocycles(self.mol)

    def count_aliphatic_heterocycles(self):
        return Descriptors.NumAliphaticHeterocycles(self.mol)

    def count_aromatic_carbocycles(self):
        return Descriptors.NumAromaticCarbocycles(self.mol)

    def count_aromatic_heterocycles(self):
        return Descriptors.NumAromaticHeterocycles(self.mol)

    def count_saturated_carbocycles(self):
        return Descriptors.NumSaturatedCarbocycles(self.mol)

    def count_saturated_heterocycles(self):
        return Descriptors.NumSaturatedHeterocycles(self.mol)

    def count_valance_electrons(self):
        return Descriptors.NumValenceElectrons(self.mol)

    def count_radical_electrons(self):
        return Descriptors.NumRadicalElectrons(self.mol)

    def get_mass(self):
        return sum(a.GetMass() for a in self.mol.GetAtoms())

    def get_heavy_atom_mass(self):
        return Descriptors.HeavyAtomMolWt(self.mol)

    def get_hetero_atom_mass(self):
        mass = 0
        for a in self.mol.GetAtoms():
            if a.GetAtomicNum() != 6 and a.GetAtomicNum() != 1:
                mass += a.GetMass()
        return mass

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

    def get_first_order_kappa(self):
        return Descriptors.Kappa1(self.mol)

    def get_second_order_kappa(self):
        return Descriptors.Kappa2(self.mol)

    def get_third_order_kappa(self):
        return Descriptors.Kappa3(self.mol)

    def get_kier_flex(self):
        k1 = Descriptors.Kappa1(self.mol)
        k2 = Descriptors.Kappa2(self.mol)
        return k1 * k2 / len(self.mol.GetAtoms())

    def get_first_zagreb(self):
        """Returns first Zagreb index.

        First Zagreb index is a topological index base on vertex degree
        $\delta$:
        \[
            M_{1} = \sum_{i = 1}{A} \delta_{i}^{2}.
        \]
        """
        return sum(pow(len(a.GetNeighbors()), 2) for a in self.mol.GetAtoms())

    def get_second_zagreb(self):
        """Returns second Zagreb index.

        Second Zagreb index is a topological index base on vertex degree
        $\delta$:
        \[
            M_{2} =
                \sum_{i = 1}{A - 1}\sum_{j = i + 1}{A}
                    a_{ij} \delta_{i} \delta{j}.
        \]
        where $a_{ij}$ are elemensts of adjacency matrix.
        """
        a = Chem.GetAdjacencyMatrix(self.mol)

        natoms = len(self.mol.GetAtoms())
        summ = 0
        for i in range(natoms - 1):
            di = len(self.mol.GetAtomWithIdx(i).GetNeighbors())
            for j in range(i + 1, natoms):
                dj = len(self.mol.GetAtomWithIdx(j).GetNeighbors())
                summ += a[i][j] * di * dj
        return summ

    def get_eccentric_connectivity(self):
        """Returns eccentric connectivity index.

        Eccentric connectivity index of a H-depleted molecule graph is defined
        as
        \[
            \zeta^{c} = \sum_{i = 1}^{A} \eta_{i} \delta_{i}
        \]
        where $\eta_{i}$ is the eccentricity ($\max_{j}(d_{ij})$) of $i$-th
        atom and $\delta_{i}$ is the vertex degree of vertex $v_{i}$.
        """
        a = Chem.GetAdjacencyMatrix(self.mol)
        d = Chem.GetDistanceMatrix(self.mol)
        return sum(max(d[i]) * sum(a[i]) for i in range(a.shape[0]))

    def get_eccentric_distance_sum(self):
        """Returns eccentric distance sum index.

        Eccentric distance sum of a H-depleted molecule graph is defined
        as
        \[
            \zeta^{DS} = \sum_{i = 1}^{A} \eta_{i} \sigma_{i}
        \]
        where $\eta_{i}$ is the eccentricity ($\max_{j}(d_{ij})$) of $i$-th
        atom and $\sigma_{i}$ is the distance degree of vertex $v_{i}$.
        """
        d = Chem.GetDistanceMatrix(self.mol)
        return sum(max(d[i]) * sum(d[i]) for i in range(d.shape[0]))

    def get_adjacent_eccentric_distance_sum(self):
        """Returns adjacent eccentric distance sum index.

        Adjacent eccentric distance sum index of a H-depleted molecule graph is
        defined as
        \[
            \zeta^{SV} = \sum_{i = 1}^{A}
                \frac{\eta_{i} \sigma_{i}}{\delta_{i}}
        \]
        where $\eta_{i}$ is the eccentricity ($\max_{j}(d_{ij})$) of $i$-th
        atom, $\sigma_{i}$ and $\delta_{i}$ are the distance and is vertex
        degree of the vertex $v_{i}$ respectively.
        """
        a = Chem.GetAdjacencyMatrix(self.mol)
        d = Chem.GetDistanceMatrix(self.mol)
        return sum(max(d[i]) * sum(d[i]) / sum(a[i])
                   for i in range(a.shape[0]))

    def get_connective_eccentricity(self):
        """Returns connective eccentricity index.

        Connective eccentricity index of a H-depleted molecule graph is defined
        as
        \[
            C^{\zeta} = \sum_{i = 1}^{A} \frac{\delta_{i}}{\eta_{i}}
        \]
        where $\eta_{i}$ and $\delta_{i}$ are the eccentricity
        ($\max_{j}(d_{ij})$) and vertex degree of vertex $v_{i}$ respectively.
        """
        a = Chem.GetAdjacencyMatrix(self.mol)
        d = Chem.GetDistanceMatrix(self.mol)
        return sum(sum(a[i]) / max(d[i]) for i in range(a.shape[0]))

    def get_eccentric_adjacency(self):
        """Returns eccentric adjacency index.

        Eccentric adjacency index of a H-depleted molecule graph is defined
        as
        \[
            \zeta^{A} = \sum_{i = 1}^{A} \frac{EC_{i}^{1}}{\eta_{i}}
        \]
        where $\eta_{i}$ and $EC_{i}^{1}$ are the eccentricity
        ($\max_{j}(d_{ij})$) and vertex degree and extended connectivity of
        fisrt order of vertex $v_{i}$ respectively.

        Extended connectivity is defined as
        \[
            EC_{i}^{k + 1} = \sum_{j = 1}^{A} a_{ij} EC_{j}^{k}
        \]
        where $EC_{j}^{0}$ is the vertex degree of vertex $v_{j}$.
        """
        a = Chem.GetAdjacencyMatrix(self.mol)
        d = Chem.GetDistanceMatrix(self.mol)

        summ = 0.0
        for i in range(a.shape[0]):
            eci = sum(a[i][j] * sum(a[j]) for j in range(a.shape[0]))
            summ += eci / max(d[i])
        return summ

    def get_superadjacency(self):
        """Returns superadjacency index.

        Eccentric adjacency index of a H-depleted molecule graph is defined
        as
        \[
            \int^{A} = \sum_{i = 1}^{A} \frac{EC_{i}^{1} \delta_{i}}{\eta_{i}}
        \]
        where $\eta_{i}$, $EC_{i}^{1}$, and $delta_{I} are the eccentricity
        ($\max_{j}(d_{ij})$), extended connectivity of first orederm, and
        vertex degree of vertex $v_{i}$ respectively.

        Extended connectivity is defined as
        \[
            EC_{i}^{k + 1} = \sum_{j = 1}^{A} a_{ij} EC_{j}^{k}
        \]
        where $EC_{j}^{0}$ is the vertex degree of vertex $v_{j}$.
        """
        a = Chem.GetAdjacencyMatrix(self.mol)
        d = Chem.GetDistanceMatrix(self.mol)

        summ = 0.0
        for i in range(a.shape[0]):
            eci = sum(a[i][j] * sum(a[j]) for j in range(a.shape[0]))
            summ += sum(a[i]) * eci / max(d[i])
        return summ

    def get_augmented_eccentric_connectivity(self):
        """Returns augmented eccentric connectivity index.

        Eccentric adjacency index of a H-depleted molecule graph is defined
        as
        \[
            ^{A}\zeta^{c} = \sum_{i = 1}^{A} \frac{M_{i}^{1}}{\eta_{i}}
        \]
        where $\eta_{i}$ is the eccentricity ($\max_{j}(d_{ij})$) and
        \[
            $M_{i} = \prod_{j = 1}^{A} (\delta_{j})^{a_{ij}}
        \]
        with $\delta_{j}$ being the vertex degree of vertex $v_{j}$.
        """
        a = Chem.GetAdjacencyMatrix(self.mol)
        d = Chem.GetDistanceMatrix(self.mol)

        summ = 0.0
        for i in range(a.shape[0]):
            mi = 1
            for j in range(a.shape[0]):
                mi *= pow(sum(a[j]), a[i][j])
            summ += mi / max(d[i])
        return summ

    def get_fraction_CSP3(self):
        return Descriptors.FractionCSP3(self.mol)

    def get_logP(self):
        return Descriptors.MolLogP(self.mol)

    def get_mr(self):
        return Descriptors.MolMR(self.mol)

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
