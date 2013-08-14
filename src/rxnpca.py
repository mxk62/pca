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
from rdkit import Chem
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


# Create a map between functional group SMARTS and their representation as
# query molecules (used in RDKit's HasSubstructMatch() function).
groups = {}
with open('../data/functional_groups.txt', 'r') as f:
    for line in f:
        gid, smarts = line.split()
        try:
            pattern = Chem.MolFromSmarts(smarts)
        except Exception:
            continue
        groups[smarts] = pattern

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
sample = Sample(db['chemical'], rng_seed=1)

# For each chemical in the sample perform a single retrosynthetic step.
reactions = {}
for rec in sample.get(samplesize):
    smi = rec['smiles'].encode('ascii')
    chem = Chemical(smi)

    # Add existing incoming reactions of the chemical to the pool of all
    # possible reactions allowing for synthesizing the compound.
    for uid in rec['reactions']['produced']:
        smi = db['reaction'].find_one({'_id': uid}).get('smiles', None)
        if smi is not None:
            smi = smi.encode('ascii')
            reactions[smi] = Reaction(smi)
            reactions[smi].is_published = True

    # Then, iterate over available transforms to generate any possible reaction
    # leading to it, adding only unknown new ones to the pool.
    for t in [t for t in transforms if len(t.retrons) == 1]:
        for smi in chem.make_retrostep(t):
            if smi not in reactions:
                reactions[smi] = Reaction(smi)
                reactions[smi].is_published = False

# Find functional groups which are present on chemicals involved in each
# reaction, either from the database or finding them from a scratch.
field = 'functional_groups'
bad_smis = []
for smi, rxn in reactions.items():
    chemicals = rxn.reactants + rxn.products

    if None in [chem.mol for chem in chemicals]:
        bad_smis.append(smi)
        continue

    for chem in chemicals:
        rec = db.chemical.find_one({'smiles': chem.smiles})
        if rec is not None and field in rec.keys():
            chem.functional_groups = rec[field]
        else:
            chem.find_groups(groups)

# For now, remove any reaction which has None among it reactants or products.
for smi in bad_smis:
    del reactions[smi]

# Calculate reactions descriptors.
indata = []
status = []
for smi, rxn in reactions.items():
    indata.append(rxn.get_descriptors())
    #indata.append(rxn.get_group_descriptor(groups))
    status.append([int(rxn.is_published)])
save_array(indata, 'descriptors.dat')
save_array(status, 'status.dat')


# Here goes PCA analysis with help of Modular toolkit for Data Processing
# (MDP):
#     http://mdp-toolkit.sourceforge.net
#
# A node is the basic building block of an MDP application. It represents a
# data processing element, here node performs PCA.
#
# Number of principal components which will be kept can be limited by
# 'output_dim' parameter (default: None), e.g.:
#  - output_dim=10 will keep first 10 principal components
#  - output_dim=0.8 will keep the number of principal components which can
#    account for 80% of the input variance.
pcanode = mdp.nodes.PCANode(output_dim=0.95)

# The node need to be trained to perform its task. PCA algorithm requires the
# computation of the mean and covariance matrix of a set of training data from
# which the principal eigenvectors of the data distribution are estimated.
pcanode.train(array(indata))

# To make PCA algorithm to compute the selected principal eigenvectors,
# training must be finalized.
pcanode.stop_training()

# Once the training is finished, it is possible to execute the node, i.e. to
# project the input data on the principal components learned in the training
# phase.
outdata = pcanode.execute(array(indata))
for i, x in enumerate(pcanode.d, start=1):
    print i, x

save_array(outdata, 'projection.dat')
