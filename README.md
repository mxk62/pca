# rxnpca

Script is an attempt to use principal component analysis (PCA) to identify key
characteristics allowing to differ 'efficient' (working) reactions from
'inefficient' (not-working) ones.


## Description

For a random sample of chemical compounds, script employs retrosynthetic
transforms developed at NU, to generate a space of all possible reactions. 
Then, it performs PCA for a predefined set of reaction descriptors.

Currently about a 100 of molecular descriptors are available. For each
reaction a sum of values of a given molecular descriptor is calculated
separately for its reactants and products. Compounds which H-depleted molecular
graphs contains a single atom only (e.g. water 'O') are *excluded*, as many
topological descriptors are not defined for them.


Each such a pair is treated as *reaction* descriptor corresponding to a given molecular descriptor.


## Usage

To run the script type (currently all arguments are optional):

    rxnpca.py --size <size> --seed <seed> --selection-type <type>

where

+	**size**

	Sample size, i.e. number of chemicals for which reactions will be
	generated. Defaults to 1000.

+	**seed**

	Pseudo-random number generator seed. Default to `None`, meaning that
	current system time will be used.

+	**type**

	Selection method of reactions sharing the same product. Available
	choices are:
	-	`all`: all availawill be selected,
	-	`random`: a random pair of published and unpublished reactions will
		be selected for each product,
	-	`popular`: a  pair of published and unpublished reactions having the
	  	highest value of popularity index will be picked.
	If absent, default to `all`.

