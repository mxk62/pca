import sys
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

        descriptors = []

        # Calculate number of H donors.
        r_H_donors = sum(chem.count_H_donors() for chem in self.reactants)
        p_H_donors = sum(chem.count_H_donors() for chem in self.products)
        descriptors.extend([r_H_donors, p_H_donors])

        # Calculate number of H acceptors.
        r_H_acceptors = sum(chem.count_H_acceptors() for chem in self.reactants)
        p_H_acceptors = sum(chem.count_H_acceptors() for chem in self.products)
        descriptors.extend([r_H_acceptors, p_H_acceptors])

        # Calculate number of atoms.
        r_atoms = sum(chem.count_atoms() for chem in self.reactants)
        p_atoms = sum(chem.count_atoms() for chem in self.products)
        descriptors.extend([r_atoms, p_atoms])

        # Calculate number of hetero atoms.
        r_hetero_atoms = sum(chem.count_hetero_atoms() for chem in self.reactants)
        p_hetero_atoms = sum(chem.count_hetero_atoms() for chem in self.products)
        descriptors.extend([r_hetero_atoms, p_hetero_atoms])

        # Calculate number of bonds.
        r_bonds = sum(chem.count_bonds() for chem in self.reactants)
        p_bonds = sum(chem.count_bonds() for chem in self.products)
        descriptors.extend([r_bonds, p_bonds])

        # Calculate fraction of rotatable bonds.
        r_rotatable_bonds = sum(chem.get_fraction_rotatable_bonds() for chem in self.reactants)
        p_rotatable_bonds = sum(chem.get_fraction_rotatable_bonds() for chem in self.products)
        descriptors.extend([r_rotatable_bonds, p_rotatable_bonds])

        # Calculate number of rings.
        r_rings = sum(chem.count_rings() for chem in self.reactants)
        p_rings = sum(chem.count_rings() for chem in self.products)
        descriptors.extend([r_rings, p_rings])

        # Calculate aliphatic rings.
        r_aliphatic_rings = sum(chem.count_aliphatic_rings() for chem in self.reactants)
        p_aliphatic_rings = sum(chem.count_aliphatic_rings() for chem in self.products)
        descriptors.extend([r_aliphatic_rings, p_aliphatic_rings])

        # Calculate number of aromatic rings.
        r_aromatic_rings = sum(chem.count_aromatic_rings() for chem in self.reactants)
        p_aromatic_rings = sum(chem.count_aromatic_rings() for chem in self.products)
        descriptors.extend([r_aromatic_rings, p_aromatic_rings])

        # Calculate number of saturated rings.
        r_saturated_rings = sum(chem.count_saturated_rings() for chem in self.reactants)
        p_saturated_rings = sum(chem.count_saturated_rings() for chem in self.products)
        descriptors.extend([r_saturated_rings, p_saturated_rings])

        # Calculate aliphatic carbocycles.
        r_aliphatic_carbocycles = sum(chem.count_aliphatic_carbocycles() for chem in self.reactants)
        p_aliphatic_carbocycles = sum(chem.count_aliphatic_carbocycles() for chem in self.products)
        descriptors.extend([r_aliphatic_carbocycles, p_aliphatic_carbocycles])

        # Calculate aliphatic heterocycles.
        r_aliphatic_heterocycles = sum(chem.count_aliphatic_heterocycles() for chem in self.reactants)
        p_aliphatic_heterocycles = sum(chem.count_aliphatic_heterocycles() for chem in self.products)
        descriptors.extend([r_aliphatic_carbocycles, p_aliphatic_carbocycles])

        # Calculate aromatic carbocycles.
        r_aromatic_carbocycles = sum(chem.count_aromatic_carbocycles() for chem in self.reactants)
        p_aromatic_carbocycles = sum(chem.count_aromatic_carbocycles() for chem in self.products)
        descriptors.extend([r_aromatic_carbocycles, p_aromatic_carbocycles])

        # Calculate aromatic heterocycles.
        r_aromatic_heterocycles = sum(chem.count_aromatic_heterocycles() for chem in self.reactants)
        p_aromatic_heterocycles = sum(chem.count_aromatic_heterocycles() for chem in self.products)
        descriptors.extend([r_aromatic_heterocycles, p_aromatic_heterocycles])

        # Calculate saturated carbocycles.
        r_saturated_carbocycles = sum(chem.count_saturated_carbocycles() for chem in self.reactants)
        p_saturated_carbocycles = sum(chem.count_saturated_carbocycles() for chem in self.products)
        descriptors.extend([r_saturated_carbocycles, p_saturated_carbocycles])

        # Calculate saturated heterocycles.
        r_saturated_heterocycles = sum(chem.count_saturated_heterocycles() for chem in self.reactants)
        p_saturated_heterocycles = sum(chem.count_saturated_heterocycles() for chem in self.products)
        descriptors.extend([r_saturated_heterocycles, p_saturated_heterocycles])

        # Calculate number of valence electrons.
        r_valence_electrons = sum(chem.count_valance_electrons() for chem in self.reactants)
        p_valence_electrons = sum(chem.count_valance_electrons() for chem in self.products)
        descriptors.extend([r_valence_electrons, p_valence_electrons])

        # Calculate number of radical electrons.
        r_radical_electrons = sum(chem.count_radical_electrons() for chem in self.reactants)
        p_radical_electrons = sum(chem.count_radical_electrons() for chem in self.products)
        descriptors.extend([r_radical_electrons, p_radical_electrons])

        # Calculate masses.
        r_mass = sum(chem.get_mass() for chem in self.reactants)
        p_mass = sum(chem.get_mass() for chem in self.products)
        descriptors.extend([r_mass, p_mass])

        # Calculate average molecular weight.
        r_heavy_atom_mol_wt = sum(chem.get_heavy_atom_mass() for chem in self.reactants)
        p_heavy_atom_mol_wt = sum(chem.get_heavy_atom_mass() for chem in self.products)
        descriptors.extend([r_heavy_atom_mol_wt, p_heavy_atom_mol_wt])

        # Calculate total mass of hetero atoms.
        r_hetero_atom_mass = sum(chem.get_hetero_atom_mass() for chem in self.reactants)
        p_hetero_atom_mass = sum(chem.get_hetero_atom_mass() for chem in self.products)
        descriptors.extend([r_hetero_atom_mass, p_hetero_atom_mass])

        # Calculate information content of the coefficients of the
        #characteristic polynomial of adjancency matrix.
        r_ipc = sum(chem.get_ipc() for chem in self.reactants)
        p_ipc = sum(chem.get_ipc() for chem in self.products)
        descriptors.extend([r_ipc, p_ipc])
        
        # Calculate Randic indices.
        r_randic = sum(chem.get_randic() for chem in self.reactants)
        p_randic = sum(chem.get_randic() for chem in self.products)
        descriptors.extend([r_randic, p_randic])

        # Calculate Balban J indices.
        r_balaban = sum(chem.get_balaban() for chem in self.reactants)
        p_balaban = sum(chem.get_balaban() for chem in self.products)
        descriptors.extend([r_balaban, p_balaban])

        # Calculate Bertz indices.
        r_bertz = sum(chem.get_bertz() for chem in self.reactants)
        p_bertz = sum(chem.get_bertz() for chem in self.products)
        descriptors.extend([r_bertz, p_bertz])
        
        # Calculate Wiener indices.
        r_wiener = sum(chem.get_wiener() for chem in self.reactants)
        p_wiener = sum(chem.get_wiener() for chem in self.products)
        descriptors.extend([r_wiener, p_wiener])

        # Calculate Kier flexibility indices.
        r_kier_flex = sum(chem.get_kier_flex() for chem in self.reactants)
        p_kier_flex = sum(chem.get_kier_flex() for chem in self.products)
        descriptors.extend([r_kier_flex, p_kier_flex])

        # Calculate fraction of C atoms SP3 hybridized.
        r_CSP3 = sum(chem.get_fraction_CSP3() for chem in self.reactants)
        p_CSP3 = sum(chem.get_fraction_CSP3() for chem in self.products)
        descriptors.extend([r_CSP3, p_CSP3])

        # Calculate Wildman-Crippen logP values.
        r_log_p = sum(chem.get_logP() for chem in self.reactants)
        p_log_p = sum(chem.get_logP() for chem in self.products)
        descriptors.extend([r_log_p, p_log_p])

        # Calculate Wildman-Crippen mr values.
        r_mr = sum(chem.get_mr() for chem in self.reactants)
        p_mr = sum(chem.get_mr() for chem in self.products)
        descriptors.extend([r_mr, p_mr])
        
        # Calculate induction parameter of molecule.
        r_induction_parameter = sum(chem.get_induction_parameter() for chem in self.reactants)
        p_induction_parameter = sum(chem.get_induction_parameter() for chem in self.products)
        descriptors.extend([r_induction_parameter, p_induction_parameter])
        
        # Calculate Topological Polar Surface Area of molecule.
        r_TPSA = sum(chem.get_TPSA() for chem in self.reactants)
        p_TPSA = sum(chem.get_TPSA() for chem in self.products)
        descriptors.extend([r_TPSA, p_TPSA])
        
        # Calculate ring fusion density
        r_RF = sum(chem.get_RF_delta() for chem in self.reactants)
        p_RF = sum(chem.get_RF_delta() for chem in self.products)
        descriptors.extend([r_RF, p_RF])
        
        return descriptors

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
