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


# Initialize connection with the database.
try:
    client = MongoClient()
except ConnectionFailure:
    sys.stderr.write('Error: cannot connect to the database.')
    sys.exit(1)
db = client['data']

# Initialize available transforms. Ignore ones having more than one retron.
transforms = []
for rec in db['retro'].find():
    try:
        t = Transform(rec['reaction_smarts'].encode('ascii'))
    except ValueError:
        continue
    if len(t.retrons) == 1:
        transforms.append(t)

# Get a random sample of chemical compounds.
sample = Sample(db['chemical'], size=1000, rng_seed=1)

# For each chemical in the sample perform a single retrosynthetic step.
reactions = {}
for chem_rec in sample.get():
    smi = chem_rec['smiles'].encode('ascii')
    try:
        chem = Chemical(smi)
    except ValueError:
        sys.stderr.write('Invalid chemical SMILES: {0}. Ignoring'.format(smi))
        continue

    # Ignore single element reactions, e.g. ionization. Many descriptors
    # cannot be calculated for them since their adjacency and distance
    # matrices are equal to zero in such cases.
    if len(chem.mol.GetAtoms()) == 1:
        continue

    # For chemical at hand, iterate over available transforms to generate any
    # possible incoming reactions leading to it.
    new_rxns = {}
    for t in transforms:
        possible = chem.make_retrostep(t)
        for smi in set(possible).difference(reactions):
            try:
                new_rxns[smi] = Reaction(smi)
            except ValueError:
                continue
    reactions.update(new_rxns)

    # Then use database to find out which of those reactions are already
    # published.
    for uid in chem_rec['reactions']['produced']:
        rxn_rec = db['reaction'].find_one({'_id': uid})
        if rxn_rec is None:
            continue

        # RX.ID are not unique thus the 'rxid' field in the database is a
        # list. Pick the first eleemnt, if not empty. Use -1 to mark
        # published reaction without a known RX.ID.
        rxid = rxn_rec.get('rxid', [-1])[0]

        old_smi = rxn_rec.get('smiles', None)
        if old_smi is not None:
            try:
                old_rxn = Reaction(old_smi.encode('ascii'), rxnid=rxid)
            except ValueError:
                continue
            for new_rxn in [r for r in new_rxns.values() if r == old_rxn]:
                del reactions[new_rxn.smiles]
                reactions[old_rxn.smiles] = old_rxn


# Find functional groups which are present on chemicals involved in each
# reaction, either from the database or finding them from a scratch.
#
# Create a map between known functional group SMARTS and their
# representation as query molecules (used in RDKit's HasSubstructMatch()
# function).
#groups = {}
#with open('../data/functional_groups.txt', 'r') as f:
#    for line in f:
#        gid, smarts = line.split()
#        try:
#            pattern = Chem.MolFromSmarts(smarts)
#        except Exception:
#            continue
#        groups[smarts] = pattern
#
#field = 'functional_groups'
#for smi, rxn in reactions.items():
#    for chem in rxn.reactants + rxn.products:
#        rec = db.chemical.find_one({'smiles': chem.smiles})
#        if rec is not None and field in rec.keys():
#            chem.functional_groups = rec[field]
#        else:
#            chem.find_groups(groups)

# Calculate reactions descriptors.
descriptors = []
rxids = []
smiles = []
status = []
for smi, rxn in reactions.items():
    try:
        descriptors.append(rxn.get_descriptors())
    except RuntimeWarning:
        sys.stderr.write('{}'.format(rxn.smiles))
    #indata.append(rxn.get_group_descriptor(groups))
    status.append([1 if rxn.rxnid is not None else 0])
    smiles.append([rxn.smiles])
    rxids.append([rxn.rxnid])
save_array(descriptors, 'descriptors.dat')
save_array(rxids, 'rxids.dat')
save_array(smiles, 'smiles.dat')
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
#pcanode = mdp.nodes.PCANode(output_dim=0.95)

# The node need to be trained to perform its task. PCA algorithm requires the
# computation of the mean and covariance matrix of a set of training data from
# which the principal eigenvectors of the data distribution are estimated.
#pcanode.train(array(indata))

# To make PCA algorithm to compute the selected principal eigenvectors,
# training must be finalized.
#pcanode.stop_training()

# Once the training is finished, it is possible to execute the node, i.e. to
# project the input data on the principal components learned in the training
# phase.
#outdata = pcanode.execute(array(indata))
#for i, x in enumerate(pcanode.d, start=1):
#    print i, x

#save_array(outdata, 'projection.dat')
