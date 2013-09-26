import random
import sys
from pymongo import MongoClient
from pymongo.errors import ConnectionFailure


class Sample:
    """Represents a random sample from a given database collection."""

    def __init__(self, collection, size=10, rng_seed=None):
        """Sets a database collection which will be used."""

        # Find out database collection size.
        self.coll = collection
        self.collsize = collection.count()

        # Initialize pseudo-random number generator.
        self.rng = random.Random()
        self.rng.seed(rng_seed)

        self.sample = {}
        self.fetch(size)

    def clear(self):
        """Clears the conentent of the sample."""
        self.sample.clear()

    def fetch(self, size):
        """Fetches a random sample from the collection.

        Populates sample having 'size' of random, non-repeating elements.
        Each successive call, replaces previous sample content.
        """
        size = min(size, self.collsize)
        if self.sample:
            self.clear()
        while len(self.sample) < size:
            uid = self.rng.randint(0, self.collsize)
            if uid in self.sample:
                continue

            record = self.coll.find_one({'_id': uid})
            if record is not None:
                self.sample[uid] = record

    def get(self):
        """Returns the content of the sample."""
        return self.sample.values()

    def read(self, filename):
        """Read and return the sample saved in a file."""
        self.clear()
        with open(filename, 'r') as f:
            line = f.readline()
            for uid in [int(s) for s in line.strip().split('\t')]:
                record = self.coll.find_one({'_id': uid})
                if record is not None:
                    self.sample[uid] = record
        return self.sample.values()

    def save(self, filename):
        """Write the sample to a file."""
        with open(filename, 'w') as f:
            f.write('\t'.join(uid for uid in self.sample.keys()))


class ExampleData:
    """Represents test data."""

    def __init__(self, document):
        self.doc = document

    def get_descriptors(self):
        return [float(x) for x in self.doc['data']]


if __name__ == '__main__':

    # Initialize connection with the database.
    try:
        client = MongoClient()
    except ConnectionFailure:
        sys.stderr.write('Error: connection failure. Exiting.')
        sys.exit(1)

    # Create a test database if it does not exists.
    data = [(2.5, 2.4), (0.5, 0.7), (2.2, 2.9), (1.9, 2.2), (3.1, 3.0),
            (2.3, 2.7), (2.0, 1.6), (1.0, 1.1), (1.5, 1.6), (1.1, 0.9)]

    dbs = client.database_names()
    if 'pca' in dbs:
        client.drop_database('pca')
    db = client['pca']
    coll = db['data']
    recs = [{"_id": uid, "data": list(entry)}
            for uid, entry in enumerate(data)]
    coll.insert(recs)

    # Retrive a sample from the collection and print entry descriptors.
    sample = Sample(coll)
    for doc in sample.get():
        datum = ExampleData(doc)
        print datum.get_descriptors()
