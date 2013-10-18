#!/usr/bin/env python2
"""
rxnpair.py

10/11/2013
Syeda Sabrina

This script groups all reactions in training set with same product.
After grouping the reaction, it returns a list of reaction pairs for 
validation of scoring function. In each pair, one is published and
another is unpublished reaction. Both reactions in each pair produce
same product.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from random import shuffle


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
    
    # Load list of reaction smiles and their status.
    rxn_smiles = loadList('../src/smiles.dat')
    rxn_status = loadList('../src/status.dat')
    
    # Get reactants and products smiles.
    products = {}       # key: reaction row number and value: product smiles
    productList = []    # list of products.
    for idx,rsmiles in enumerate(rxn_smiles):
        # idx+1 is the row number of reaction in the training data.
        # rsmiles[0] is a reaction's smiles.
        reactantsmi, productsmi =[smi.split('.') for smi in rsmiles[0].split('>>')]
        # productsmi is list of products' smiles in the reaction
        products[idx+1] = productsmi    # key: reaction row number value: list of products smiles.
        for p in productsmi:
            productList.append(p)
    
    # Group reactions with same product/products.
    reactions_group = {} # value of this dictionary is: row number, that 
                         # represents the reaction in training set. 
    rsmiles_group = {}  # value of this dictionary reaction smiles
    for p in productList:
        reactions_group[p] = []
        rsmiles_group[p] = []
        for k in products.keys():
            if p in products[k]:
                reactions_group[p].append(k)
                rsmiles_group[p].append(rxn_smiles[k-1])
    
    # Make reaction pair 
    rxn_pair = {}
    for p in reactions_group.keys():
        rxns = reactions_group[p]
        l = len(rxns)
        rxn_pair[p] = []
        x = list(range(l))
        shuffle(x)
        for i in x:
            r = rxns[i]
            if int(rxn_status[r-1][0]) == 0:
               rxn_pair[p].append(r)
               break
        for i in x:
            r = rxns[i]
            if int(rxn_status[r-1][0]) == 1:
               rxn_pair[p].append(r)
               break
                
    for k in rxn_pair.keys():
        if len(rxn_pair[k]) == 1:
            del (rxn_pair[k])

 
    # Save results.            
    pairedrxn = []
    for psmi,row in rxn_pair.items():
        pairedrxn.append(row)
    #for psmi,row in reactions_group.items():
     #  groupedreactions.append(row)
    #for psmi,rsmiles in rsmiles_group.items():
     #   groupedsmiles.append(rsmiles)

   # In the training set each row represent a reaction. So, row number can be used as
   # a reaction ID. Hence,
   # In the following two files each row represent all reactions in training set
   # having same product. The reactions are represented by their smiles in 
   # groupedsmiles.dat file and in groupedreactions.dat file the reactions are
   # represented by their row in the training set at each row in training set
   # represent a reaction.
    save_array(pairedrxn,'rxn_pair.txt')
    #save_array(groupedsmiles,'groupedsmiles.dat')