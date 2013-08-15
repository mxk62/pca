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
        
        #Calculate fraction of C atoms SP3 hybridized.
        r_CSP3 = sum(chem.get_fraction_CSP3() for chem in self.reactants)
        p_CSP3 = sum(chem.get_fraction_CSP3() for chem in self.products)
        
        #Calculate aliphatic carbocycles.
        r_aliphatic_carbocycles = sum(chem.get_num_aliphatic_carbocycles() for chem in self.reactants)
        p_aliphatic_carbocycles = sum(chem.get_num_aliphatic_carbocycles() for chem in self.products)
        
        #Calculate aliphatic heterocycles.
        r_aliphatic_heterocycles = sum(chem.get_num_aliphatic_heterocycles() for chem in self.reactants)
        p_aliphatic_heterocycles = sum(chem.get_num_aliphatic_heterocycles() for chem in self.products)
        
        #Calculate saturated rings.
        r_saturated_rings = sum(chem.get_num_saturated_rings() for chem in self.reactants)
        p_saturated_rings = sum(chem.get_num_saturated_rings() for chem in self.products)
        
        #Calculate saturated heterocycles.
        r_saturated_heterocycles = sum(chem.get_num_saturated_heterocycles() for chem in self.reactants)
        p_saturated_heterocycles = sum(chem.get_num_saturated_heterocycles() for chem in self.products)
        
        #Calculate saturated carbocycles.
        r_saturated_carbocycles = sum(chem.get_num_saturated_carbocycles() for chem in self.reactants)
        p_saturated_carbocycles = sum(chem.get_num_saturated_carbocycles() for chem in self.products)
        
        #Calculate aromatic rings.
        r_aromatic_rings = sum(chem.get_num_aromatic_rings() for chem in self.reactants)
        p_aromatic_rings = sum(chem.get_num_aromatic_rings() for chem in self.products)
        
        #Calculate aromatic heterocycles.
        r_aromatic_heterocycles = sum(chem.get_num_aromatic_heterocycles() for chem in self.reactants)
        p_aromatic_heterocycles = sum(chem.get_num_aromatic_heterocycles() for chem in self.products)
        
        #Calculate aromatic carbocycles.
        r_aromatic_carbocycles = sum(chem.get_num_aromatic_carbocycles() for chem in self.reactants)
        p_aromatic_carbocycles = sum(chem.get_num_aromatic_carbocycles() for chem in self.products)
        
        #Calculate aliphatic rings.
        r_aliphatic_rings = sum(chem.get_num_aliphatic_rings() for chem in self.reactants)
        p_aliphatic_rings = sum(chem.get_num_aliphatic_rings() for chem in self.products)
        
        #Calculate Number of H donors.
        r_H_donors = sum(chem.get_num_H_donors() for chem in self.reactants)
        p_H_donors = sum(chem.get_num_H_donors() for chem in self.products)
        
        #Calculate Number of H acceptors.
        r_H_acceptors = sum(chem.get_num_H_acceptors() for chem in self.reactants)
        p_H_acceptors = sum(chem.get_num_H_acceptors() for chem in self.products)
        
        #Calculate fraction of hetero atoms.
        r_hetero_atoms = sum(chem.get_fraction_hetero_atoms() for chem in self.reactants)
        p_hetero_atoms = sum(chem.get_fraction_hetero_atoms() for chem in self.products)
        
        #Calculate hetero atoms' total mass.
        r_hetero_atoms_mass = sum(chem.get_hetero_atoms_mass() for chem in self.reactants)
        p_hetero_atoms_mass = sum(chem.get_hetero_atoms_mass() for chem in self.products)
        
        #Calculate fraction of rotatable bonds.
        r_rotatable_bonds = sum(chem.get_fraction_rotatable_bonds() for chem in self.reactants)
        p_rotatable_bonds = sum(chem.get_fraction_rotatable_bonds() for chem in self.products)
        
        #Calculate number of valence electrons.
        r_valence_electrons = sum(chem.get_num_valance_electrons() for chem in self.reactants)
        p_valence_electrons = sum(chem.get_num_valance_electrons() for chem in self.products)
        
        #Calculate number of radical electrons.
        r_radical_electrons = sum(chem.get_num_radical_electrons() for chem in self.reactants)
        p_radical_electrons = sum(chem.get_num_radical_electrons() for chem in self.products)
        
        #Calculate average molecular weight ignoring hydrogens.
        r_heavy_atom_mol_wt = sum(chem.get_heavy_atom_mol_wt() for chem in self.reactants)
        p_heavy_atom_mol_wt = sum(chem.get_heavy_atom_mol_wt() for chem in self.products)

        return [r_atoms, p_atoms, r_bonds, p_bonds, r_rings, p_rings,
                r_mass, p_mass, r_balaban, p_balaban, r_bertz, p_bertz,
                r_kier_flex, p_kier_flex, r_randic, p_randic, r_CSP3, p_CSP3,
                r_aliphatic_carbocycles, p_aliphatic_carbocycles, r_aliphatic_heterocycles,
                p_aliphatic_heterocycles, r_saturated_rings, p_saturated_rings,
                 r_saturated_heterocycles,  p_saturated_heterocycles, r_saturated_carbocycles,
                 p_saturated_carbocycles, r_aromatic_rings, p_aromatic_rings,
                  r_aromatic_heterocycles,  p_aromatic_heterocycles, r_aromatic_carbocycles,
                  p_aromatic_carbocycles, r_aliphatic_rings, p_aliphatic_rings, r_H_donors,
                  p_H_donors, r_H_acceptors, p_H_acceptors, r_hetero_atoms_mass,
                  p_hetero_atoms_mass, r_hetero_atoms, p_hetero_atoms, r_rotatable_bonds,
                  p_rotatable_bonds, r_valence_electrons, p_valence_electrons,
                  r_radical_electrons, p_radical_electrons, r_heavy_atom_mol_wt,
                  p_heavy_atom_mol_wt]

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
