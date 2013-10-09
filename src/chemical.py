import math
import numpy
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdchem import BondType


class Chemical:
    """Represents a chemical compound."""

    def __init__(self, smiles):
        self.smiles = smiles.strip()
        try:
            self.mol = Chem.MolFromSmiles(self.smiles)
        except Exception:
            print 'Error: invalid compound SMILES: {}.'.format(self.smiles)
            sys.exit(1)
        self.a = Chem.GetAdjacencyMatrix(self.mol)
        self.d = Chem.GetDistanceMatrix(self.mol)
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
        return Descriptors.Ipc(self.mol)

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
        deltas = numpy.sum(self.a, axis=(1,))
        natoms = len(self.a)

        summ = 0
        for i in range(natoms - 1):
            for j in range(i + 1, natoms):
                summ += self.a[i][j] * deltas[i] * deltas[j]
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
        return sum(self.d[i].max() * self.a[i].sum()
                   for i in range(self.a.shape[0]))

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
        return sum(self.d[i].max() * self.d[i].sum()
                   for i in range(self.d.shape[0]))

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
        return sum(self.d[i].max() * self.d[i].sum() / self.a[i].sum()
                   for i in range(self.a.shape[0]))

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
        return sum(self.a[i].sum() / self.d[i].max()
                   for i in range(self.a.shape[0]))

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
        natoms = len(self.a)

        summ = 0.0
        for i in range(natoms):
            eci = sum(self.a[i][j] * sum(self.a[j]) for j in range(natoms))
            summ += eci / max(self.d[i])
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
        natoms = len(self.a)

        summ = 0.0
        for i in range(natoms):
            eci = sum(self.a[i][j] * sum(self.a[j]) for j in range(natoms))
            summ += sum(self.a[i]) * eci / max(self.d[i])
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
        natoms = len(self.a)
        summ = 0.0
        for i in range(natoms):
            mi = 1
            for j in range(natoms):
                mi *= pow(sum(self.a[j]), self.a[i][j])
            summ += mi / max(self.d[i])
        return summ

    def get_fraction_CSP3(self):
        return Descriptors.FractionCSP3(self.mol)

    def get_logP(self):
        return Descriptors.MolLogP(self.mol)

    def get_mr(self):
        return Descriptors.MolMR(self.mol)

    def get_schultz(self):
        """Returns Schultz molecular topological index.

        Schultz molecular topological index is defined as
        \[
            MTI = \sum_{i = 1}^{A}
                [(\mathbf(A) + \mathbf(D)) \dot \mathbf{v}]_{i} =
                    \sum_{i = 1}^{A} t_{i}
        \]
        with $\mathbf{A}$, $\mathbf{D}$ being adjacency and distance matrix
        respectively, and $\mathbf{v}$ being  $A$-dimensional column vector
        constitued by the vertex degrees of the atoms in the H-depleted
        molecular graph.
        """
        v = numpy.sum(self.a, axis=(1,))
        return numpy.sum(numpy.dot(numpy.sum(self.a, self.d), v))

    def get_wiener(self):
        """Returns Wiener index.

        Wiener index is defined as
        \[
            W = \frac{1}{2} \sum_{i}^{A} \sum_{j}^{A} d_{ij}
        \]
        $d_{ij}$ are elements of the distance martix of H-depleted molecular
        graph.
        """
        return 0.5 * numpy.sum(self.d)

    def get_induction_parameter(self):
        """Returns induction parameter of a molecule.

        Induction paramter of a molecule is defined as
        \[
            q_{ind} = 1 - \frac{N_{d}}{A},
        \]
        where $N_{d}$ is the number of double bonds and $A$ is the number of
        atoms in a molecule.
        """
        n_atom = len(self.mol.GetAtoms())
        n_double = len([b for b in self.mol.GetBonds()
                        if b.GetBondType() == BondType.DOUBLE])
        return 1.0 - n_double / n_atom

    def get_TPSA(self):
        return Descriptors.TPSA(self.mol)

    def get_ring_fusion_density(self):
        """Returns ring fusion density of a molecule.

        Ring fusion density is defined as
        \[
            RF_{\delta} = 2 \frac{R_{b}}{A_{r}}
        \]
        where $R_{b}$ is number of ring bridges and $A_{r}$ is the total number
        of atoms belonging to ring systems.
        """
        ri = self.mol.GetRingInfo()

        # Calculate number of ring bridges.
        Rb = len([b for b in self.mol.GetBonds()
                  if ri.NumBondRings(b.GetIdx()) > 1])

        # Calculate number of atoms in ring systems.
        Ar = len([a for a in self.mol.GetAtoms() if a.IsInRing()])

        return 2 * Rb / Ar if Ar != 0 else 0

    def get_cyclized_degree(self):
        """Returns molecular cyclized degree.

        Molecular cyclized degree is defined as
        \[
            MCD = \frac{A_{R}}{A}
        \]
        where $A_{R}$ is the total number of taoms belonging to any ring
        system.
        """
        A = self.mol.GetNumAtoms()
        Ar = len([a for a in self.mol.GetAtoms() if a.IsInRing()])
        return Ar / float(A)

    def get_ring_complexity(self):
        """Returns ring complexity index.

        Ring complexity index is defined as
        \[
            C_{R} = \frac{R}{A_{R}}
        \]
        where $R$ is total ring size and $A_{r}$ is the total number of atoms
        belonging to any ring system. For isolated rings $C_{R} = 1$, for fused
        or bridged ring system $C_{R} > 1$, for molecules with no rings $C_{R}
        = 0$.
        """
        ri = self.mol.GetRingInfo()

        # Calculate number of atoms in ring systems.
        Ar = len([a for a in self.mol.GetAtoms() if a.IsInRing()])

        # Calculate total ring size.
        R = sum([ri.NumAtomRings(a.GetIdx()) for a in self.mol.GetAtoms()])

        return R / float(Ar) if Ar != 0 else 0

    def get_total_information_content(self):
        """Returns total information content on the adjacency equality.

        The total information content on the adjacency equality is definde as
        \[
            ^{V}I_{adj}^{E} = A^{2} \log_{2}A^{2} - 2 B \log_{2}2B -
                (A^{2} - 2 B) \log_{2}(A^2 - 2 B)
        \]
        """
        nv = len(self.a)
        ne = numpy.sum(self.a) / 2

        nvsqr = math.pow(nv, 2)
        return (nvsqr * math.log(nvsqr, 2) - 2 * ne * math.log(2 * ne, 2) -
                (nvsqr - 2 * ne) * math.log(nvsqr - 2 * ne, 2))

    def get_total_information_on_atomic_composition(self):
        """Returns total information index on atomic composition.

        Information index on atomic composition is defined as
        \[
            I_{AC} = A^{h} \log_{2}A^{h} - \sum_{g = 1}^{G} A_{g} \log_{2}A_{g}
        \]
        where $A^{h}$ is the total number of atoms (hydrogen included) in a
        molecule and $A_{g}$ is the number of atoms of chemical element of type
        $g$.
        """
        atoms = self.mol.GetAtoms()
        natoms = len(atoms)
        nhs = 0

        # Count atoms of different types (including hydrogens).
        atypes = {}
        for a in atoms:
            atype = a.GetSymbol().upper()
            atypes.setdefault(atype, 0)
            atypes[atype] += 1

            nhs += a.GetTotalNumHs()

        return ((natoms + nhs) * math.log(natoms + nhs, 2) -
                sum(atypes[key] * math.log(atypes[key], 2)
                    for key in atypes.keys()))

    def get_information_bond(self):
        """Returns information bond index.

        Information bond index is defined as
        \[
            I_{B} = B \log_{2}B - \sum_{g = 1}^{G} B_{g} \log_{2}B_{g}
        \]
        where $B$ is the number of bonds in a molecule and $B_{g}$ is the
        number of bonds of type $g$; summation goes over all $G$ different
        types of bond in the molecule (single, double, triple, aromatic).
        """
        bonds = self.mol.GetBonds()
        nbonds = len(bonds)

        # Count bonds of different types.
        btypes = {}
        for b in bonds:
            btype = b.GetBondType()
            btypes.setdefault(btype, 0)
            btypes[btype] += 1

        return (nbonds * math.log(nbonds, 2) -
                sum(btypes[key] * math.log(btypes[key], 2)
                    for key in btypes.keys()))

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
