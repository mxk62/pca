#/usr/bin/env python2

"""
rxnpca.py

Scripts performs principal component analysis (PCA) for a given sample of
chemical reactions and their descriptors.
"""

import sys
import mdp
from numpy import array
from pymongo import MongoClient
from pymongo.errors import ConnectionFailure
from sample import Sample
from reaction import Reaction, Transform
from chemical import Chemical


def save_array(a, filename):
    """Saves array in a file."""
    with open(filename, 'w') as f:
        for row in a:
            for e in row:
                f.write(str(e) + ' ')
            f.write('\n')


def is_published(db, rxn):
    """Returns true if reaction is published.

    Functions scans a collection looking for a reaction having exactly the same
    SMILES as input reaction. Returns True, if found.
    """
    doc = db.chemical.find_one({'smiles': rxn.prod_smis[0]})
    if doc is not None:
        for uid in doc['reactions']['produced']:
            if db.reaction.find_one({'reaction_smarts': rxn.smi}) is not None:
                return 1
    return 0


# Initialize connection with the database.
try:
    client = MongoClient()
except ConnectionFailure:
    sys.err.write('Error: cannot connect to the database.')
    sys.exit()
db = client['data']

# Initialize available transforms.
transforms = []
for rec in db['retro'].find():
    t = Transform(rec['reaction_smarts'].encode('ascii'))
    transforms.append(t)

# Get a random sample of chemical compounds.
samplesize = 100
sample = Sample(db['chemical'])

# For each chemical perform a single retrosynthetic step.
reactions = set([])
for rec in sample.get(samplesize):
    chem = Chemical(rec['smiles'].encode('ascii'))
    for transform in transforms:
        if len(transform.retrons) != 1:
            continue
        rxn_smis = chem.make_retrostep(transform)
        if rxn_smis is not None:
            reactions.update(rxn_smis)

# Calculate reactions descriptors.
indata = []
ispublished = []
for smi in reactions:
    rxn = Reaction(smi)
    ispublished.append(is_published(db, rxn))
    indata.append(rxn.get_descriptors())
save_array(indata, 'descriptors.dat')

with open('ispublished.dat', 'w') as f:
    for e in ispublished:
        f.write(str(e) + '\n')

# Here goes PCA analysis.
outdata = mdp.pca(array(indata))
save_array(outdata, 'pca.dat')
