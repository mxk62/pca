#!/usr/bin/env python2
"""
rxnsort.py

10/11/2013
Syeda Sabrina

This script returns two files after grouping reactions with same product.

"""
from rdkit import Chem
from rdkit.Chem import AllChem



def loadList(filename):
    """Returns list of reaction smiles from a file."""

    reaction_smiles = []
    for line in open(filename):
        smiles = line.strip().split('\t')
        reaction_smiles.append(smiles)
    return reaction_smiles

def save_array(a, filename):
    """Saves array in a file."""
    with open(filename, 'w') as f:
        for row in a:
            for e in row:
                f.write(str(e) + ' ')
            f.write('\n')

if __name__ == "__main__":
    
    # Load list of reaction smiles.
    rxn_smiles = loadList('../src/smiles.dat')
    
    # Get reactants and products smiles.
    products = {}
    productList = []
    for idx,rsmiles in enumerate(rxn_smiles):
        for s in rsmiles:
            reactantsmi, productsmi =[smi.split('.') for smi in s.split('>>')]
            products[idx+1] = productsmi
            for p in productsmi:
                productList.append(p)
    
    # Group reactions with same product/products.
    reactions_group = {}
    rsmiles_group = {}  # key of this dictionary is common products smiles
                        # value of this dictionary reaction smiles./row number of the reaction in training set
    for p in productList:
        reactions_group[p] = []
        rsmiles_group[p] = []
        for k in products.keys():
            if p in products[k]:
                reactions_group[p].append(k)
                rsmiles_group[p].append(rxn_smiles[k-1])
                
    # Save results.            
    groupedreactions =[]
    groupedsmiles = []
    for psmi,row in reactions_group.items():
       groupedreactions.append(row)
    for psmi,rsmiles in rsmiles_group.items():
        groupedsmiles.append(rsmiles)
   # In the training set each row represent a reaction. So, row number can be used as
   # a reaction ID. Hence,
   # In the following two files each row represent all reactions in training set
   # having same product. The reactions are represented by their smiles in 
   # groupedsmiles.dat file and in groupedreactions.dat file the reactions are
   # represented by their row in the training set at each row in training set
   # represent a reaction.
    save_array(groupedreactions,'groupedreactions.dat')
    save_array(groupedsmiles,'groupedsmiles.dat')

            
