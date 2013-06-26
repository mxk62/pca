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
samplesize = 1000
sample = Sample(db['chemical'])

# For each chemical perform a single retrosynthetic step.
reactions = {}
for rec in sample.get(samplesize):
    chem = Chemical(rec['smiles'].encode('ascii'))

    # Add existing incoming reactions of the chemical to the pool of all
    # possible reactions allowing for sythnesizing the compound.
    for uid in rec['reactions']['produced']:
        smi = db['reaction'].find_one({'_id': uid}).get('smiles', None)
        if smi is not None:
            reactions[smi.encode('ascii')] = 1

    # Then, iterate over available transforms to generate any possible reaction
    # leading to it, adding only unkonwn new ones to the pool.
    #
    # Note:
    # Using setdefault() method enusres that existing reactions, already in the
    # pool, will not their status overwritten.
    for transform in transforms:
        if len(transform.retrons) != 1:
            continue
        for smi in chem.make_retrostep(transform):
            reactions.setdefault(smi, 0)

# Calculate reactions descriptors.
indata = []
status = []
for smi, stat in reactions.items():
    rxn = Reaction(smi)
    indata.append(rxn.get_descriptors())
    status.append([stat])
save_array(indata, 'descriptors.dat')
save_array(status, 'status.dat')

# Here goes PCA analysis.
outdata = mdp.pca(array(indata))
save_array(outdata, 'pca.dat')
