# rxnpca

Script is an attempt to use principal component analysis (PCA) to identify key
characteristics allowing to differ 'efficient' (working) reactions from
'inefficient' (not-working) ones.


## Description

For a random sample of chemical compounds, script employs retrosynthetic
transforms developed at NU, to generate a space of all possible reactions. 
Then, it performs PCA for a predefined set of reaction descriptors.

Currently only two toy reaction descriptors are used: number of carbons in
reaction SMILES and molecular mass of reactants. 


## Remarks

Currently, script is in an early alpha stage. Expect crashes regularly.
