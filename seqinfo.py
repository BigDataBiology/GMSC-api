import pandas as pd
import gzip
import numpy as np
from fasta_reader import IndexedFastaReader
from fna2faa_gmsc import translate

BASE_DIR = 'gmsc-db/'
MAX_FILTER_RESULTS = 100

def with_digits(prefix, n):
    n = f'{n:09}'
    return f'{prefix}.{n[:3]}_{n[3:6]}_{n[6:9]}'


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
                                    ).squeeze()
        self.habitat_ix = np.load(f'{BASE_DIR}/{database}_habitat.npy', mmap_mode='r')

        self.taxonomy = pd.read_table(f'{BASE_DIR}/{database}_ref_taxonomy_index.tsv',
                                    index_col=0,
                                    header=None,
                                    names=['seq_ix', 'taxonomy']
                                    ).squeeze()
        self.taxonomy_ix = np.load(f'{BASE_DIR}/{database}_taxonomy.npy', mmap_mode='r')

        hqs = [line.strip() for line in gzip.open(f'{BASE_DIR}/90AA_highquality.txt.gz', 'rt')]
        hq_ixs = [int(hq.split('.')[2]) for hq in hqs]
        self.is_hq = np.zeros(len(self.habitat_ix), dtype=bool)
        self.is_hq[list(hq_ixs)] = True


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
                'habitat': self.habitat.values[self.habitat_ix[ix]],
                'taxonomy': self.taxonomy.values[self.taxonomy_ix[ix]],
                }

    def seq_filter(self, hq_only, habitat_q, taxonomy_q):
        matches = self.habitat.str.contains(habitat_q).values[self.habitat_ix]
        if hq_only:
            matches &= self.is_hq
        if taxonomy_q is not None:
            match_taxonomy = self.taxonomy.str.contains(taxonomy_q).values[self.taxonomy_ix]
            matches &= match_taxonomy
        [ixs] = np.where(matches)
        rs = []
        for i,ix in enumerate(ixs):
            seq_id = with_digits(f'GMSC10.{self.database}', ix)
            if i < MAX_FILTER_RESULTS:
                rs.append(self.get_seqinfo(seq_id))
            else:
                rs.append({'seq_id': seq_id})
        return rs

