#/usr/bin/env python2

"""
rxnpca.py

Script performs principal component analysis (PCA) for random chemicals in the
database.
"""

import argparse
import sys
import mdp
import random
from numpy import array
from pymongo import MongoClient
from pymongo.errors import ConnectionFailure
from rdkit import Chem
from sample import Sample
from reaction import Reaction, Transform
from chemical import Chemical


def get_random(reactions):
    """Returns a random reaction."""
    return [reactions[random.randint(0, len(reactions) - 1)]]


def get_popular(reactions):
    """Returns a reaction with the highest popularity index."""
    p = [r.popularity for r in rxns]
    return [reactions[p.index(max(p))]]


def get_all(reactions):
    """Returns all reactions."""
    return reactions


def save_array(a, filename):
    """Saves array in a file."""
    with open(filename, 'w') as f:
        for row in a:
            for e in row:
                f.write(str(e) + ' ')
            f.write('\n')


# Parse command line options.
parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=int, default=None,
                    help='generator seed, defaults to None')
parser.add_argument('--size', type=int, default=1000,
                    help='sample size, defaults to 1000')
parser.add_argument('--selection-type', type=str, default='all',
                    choices=['all', 'popular', 'random'],
                    help='reaction selection method (all, popular, random);'
                         'default to \'all\'')
args = parser.parse_args()

# Initialize connection with the database.
try:
    client = MongoClient('130.203.235.117')
except ConnectionFailure:
    sys.stderr.write('Error: cannot connect to the database.')
    sys.exit(1)
db = client['data']

# Initialize available transforms. Ignore ones having more than one retron.
transforms = []
for rec in db['retro'].find():
    try:
        t = Transform(rec['reaction_smarts'].encode('ascii'),
                      dbid=rec['old_ID'])
    except ValueError:
        continue
    if len(t.retrons) == 1:
        transforms.append(t)

# Get a random sample of chemical compounds.
sample = Sample(db['chemical'], size=args.size, rng_seed=args.seed)

# Initialize data structures for stats gathering.
keys = ['accepted', 'valid', 'total']
new_stats = dict(zip(keys, len(keys) * [0]))
old_stats = dict(zip(keys, len(keys) * [0]))

# For each chemical in the sample perform a single retrosynthetic step.
reactions = {}
disconnected_count = 0
for chem_rec in sample.get():
    smi = chem_rec['smiles'].encode('ascii')

    # Ignore disconnected chemicals, e.g. such as c1cc([O-].[Na+])ccc1,
    # as we treat '.' (dot) as a compound separator.
    if '.' in smi:
        disconnected_count += 1
        continue

    try:
        chem = Chemical(smi)
    except ValueError:
        sys.stderr.write('Invalid chemical SMILES: {0}. Ignoring'.format(smi))
        continue

    # Ignore single element reactions, e.g. ionization. Many descriptors
    # cannot be calculated for them since their adjacency and distance
    # matrices are equal to zero in such cases.
    if chem.a.size == 1:
        continue

    # For chemical at hand, iterate over available transforms to generate any
    # possible incoming reactions leading to it.
    candidates = {}
    for t in transforms:
        candidates.update({s: [t.popularity,t.id] for s in chem.make_retrostep(t)})
    new_stats['total'] += len(candidates)

    # Then add valid ones which were not yet generated to the main pool.
    new_rxns = {}
    for smi in set(candidates.keys()).difference(reactions):
        try:
            rxn = Reaction(smi, popularity=candidates[smi][0], 
                           tid = candidates[smi][1])
        except ValueError:
            continue
        new_rxns[smi] = rxn
    new_stats['accepted'] += len(new_rxns)
    new_stats['valid'] += len(new_rxns)
    reactions.update(new_rxns)

    # Finally, use the database to find out which of those newly  added
    # reactions are already published.
    cursor = db['reaction'].find(
        {'_id': {'$in': chem_rec['reactions']['produced']}})
    for rxn_rec in cursor:
        old_stats['total'] += 1

        # RX.ID are not unique thus the 'rxid' field in the database is a
        # list. Pick the first element, if not empty. Use -1 to mark
        # published reaction without a known RX.ID.
        rxid = rxn_rec.get('rxid', [-1])[0]
        
        # Pick the first year from the list of publication years, otherwise
        # set it to zero if it is not available.
        year = rxn_rec.get('year', [0])[0]

        old_smi = rxn_rec.get('smiles', None)
        if old_smi is not None:
            try:
                old_rxn = Reaction(old_smi.encode('ascii'), rxnid=rxid,
                                   year=year)
            except ValueError:
                continue
            old_stats['valid'] += 1
            for new_rxn in [r for r in new_rxns.values() if r.equals(old_rxn)]:
                reactions[new_rxn.smiles].rxnid = rxid
                reactions[new_rxn.smiles].year = year
                old_stats['accepted'] += 1

# Each reaction is associated with a transform by which it was generated.
# For example T applied on P generates R1 --> P 
# Now extract all possible reactions generated by T on P.
# This is stored in dictionary 'possible'. key: transform ID and
# value is list of status of all possible reactions generated by T
possible = {}
for rxn in reactions.values():
    tid = rxn.tid
    status = 1 if rxn.rxnid is not None else 0
    possible.setdefault(tid,[status]).append(status)

# Get the statistics  
# key: Transform ID
# value: List[number of possible reactions, number of published reactions]  
tstat = {}
for tid in possible.keys():
    pos_rxns = len(possible[tid])
    p = possible[tid].count(1)
    tstat.setdefault(tid,[pos_rxns,p])
    
print "Reaction stats (accepted, valid, total):"
print " * published:", ", ".join(str(old_stats[k]) for k in keys)
print " * unpublished:", ", ".join(str(new_stats[k]) for k in keys)
print
print 'Disconnected chemicals: ', disconnected_count

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

# Sort reaction based on their products and category
#
# Sorted reactions are stored in a dictionary in which a key is a product
# SMILES and value are reactions sharing the same products. Those reactions
# are also stored in a dictionary which may have at most two keys: 'published'
# and 'unpublished'. If key exists, its corresponding value is a list of
# reactions of a given type. It other words it looks like this:
#
#     { product SMILES: { 'published': [], 'unpublished': [] } }
sorted_rxns = {}
for rxn in reactions.values():
    smi = rxn.prod_smis[0]
    sorted_rxns.setdefault(smi, {})

    category = 'unpublished' if rxn.rxnid is None else 'published'
    sorted_rxns[smi].setdefault(category, []).append(rxn)

# Discard entries which do not have either published or unpublished reactions.
sorted_rxns = {key: val for key, val in sorted_rxns.items()
               if len(sorted_rxns[key]) == 2}

# For each available product, from each category, select a reaction based on
# specified method, either all, picked at random or according to its
# popularity.
reactions = {}
method = {'all': get_all, 'popular': get_popular, 'random': get_random,}
for entry in sorted_rxns.values():
    for rxns in entry.values():
        reactions.update({r.smiles: r for r in method[args.selection_type](rxns)})

# Calculate reactions descriptors and save them to a file. Dump reactions'
# auxiliary data (id, year, popularity, etc.) too.
descriptors = []
auxiliaries = []
smiles = []
for rxn in reactions.values():
    descriptors.append(rxn.get_descriptors())

    # MUST stay commented out if functional groups are not initialized above.
    #descriptors.append(rxn.get_group_descriptor())

    status = 1 if rxn.rxnid is not None else 0
    possibility = tstat[rxn.tid][0]
    published = tstat[rxn.tid][1]
    auxiliaries.append([rxn.rxnid, rxn.year, rxn.popularity, 
                        possibility, published, status])
    smiles.append([rxn.smiles])
save_array(descriptors, 'descriptors.dat')
save_array(auxiliaries, 'auxiliaries.dat')
save_array(smiles, 'smiles.dat')

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
