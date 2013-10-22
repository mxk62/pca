# rxnpca

Script is an attempt to use principal component analysis (PCA) to identify key
characteristics allowing to differ 'efficient' (working) reactions from
'inefficient' (not-working) ones.


## Description

For a random sample of chemical compounds, script employs retrosynthetic
transforms developed at NU, to generate a space of all possible reactions. 
Then, it performs PCA for a predefined set of reaction descriptors.

Currently about a 100 of molecular descriptors are available. For each
reactions a sum of values of a given molecular descriptor is calculated
separately for its reactants and products. Each such a pair is treated as
*reaction* descriptor corresponding to a given molecular descriptor.
