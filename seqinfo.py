import pandas as pd
import numpy as np
from fasta_reader import IndexedFastaReader
from fna2faa_gmsc import translate

BASE_DIR = 'gmsc-db/'


class SeqInfo:
    def __init__(self, database):
        if database not in ('90AA', '100AA'):
            raise NotImplementedError(f'Database was {database}! Only "90AA" and "100AA" are supported')
        self.seqix = IndexedFastaReader(f'{BASE_DIR}/GMSC10.{database}.fna.xz')
        self.database = database
        self.habitat = pd.read_table(f'{BASE_DIR}/{database}_ref_multi_general_habitat_index.tsv',
                                    index_col=0,
                                    header=None,
                                    names=['seq_ix', 'habitat']
                                    ).squeeze().values
        self.habitat_ix = np.load(f'{BASE_DIR}/{database}_habitat.npy', mmap_mode='r')

        self.taxonomy = pd.read_table(f'{BASE_DIR}/{database}_ref_taxonomy_index.tsv',
                                    index_col=0,
                                    header=None,
                                    names=['seq_ix', 'taxonomy']
                                    ).squeeze().values
        self.taxonomy_ix = np.load(f'{BASE_DIR}/{database}_taxonomy.npy', mmap_mode='r')


    def get_seqinfo(self, seq_id):
        _,db,ix = seq_id.split('.')
        ix = int(ix)
        if db != self.database:
            raise IndexError(f'Only IDs for database "{self.database}" are accepted (got "{seq_id}"')

        nuc = self.seqix.get(seq_id).decode('ascii')
        return {
                "seq_id": seq_id,
                "nucleotide": nuc,
                "aminoacid": translate(nuc),
                'habitat': self.habitat[self.habitat_ix[ix]],
                'taxonomy': self.taxonomy[self.taxonomy_ix[ix]],
                }


