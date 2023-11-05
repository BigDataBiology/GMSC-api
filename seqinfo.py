import xz
import pandas as pd
import gzip
import lzma
import numpy as np
from os import path
from fna2faa_gmsc import translate
from typing import List


BASE_DIR = 'gmsc-db/'
MAX_THICK_RESULTS = 20
MAX_TOTAL_RESULTS = 1000

def with_digits(prefix, n):
    n = f'{n:09}'
    return f'{prefix}.{n[:3]}_{n[3:6]}_{n[6:9]}'


class IndexedFastaReader:
    def __init__(self, ifile):
        if ifile.endswith('.xz'):
            self.seqfile = xz.open(ifile, 'rb')
            ifile = ifile[:-len('.xz')]
        else:
            self.seqfile = open(ifile, 'rb')
        self.sindex = np.load(ifile + '.starts.npy', mmap_mode='r')

    def get(self, ix):
        self.seqfile.seek(int(self.sindex[ix]))
        data = self.seqfile.read(int(self.sindex[ix+1] - self.sindex[ix]))
        _h, seq, _empty = data.split(b'\n')
        return seq

class SeqInfo:
    def __init__(self, database):
        if database not in ('90AA', '100AA'):
            raise NotImplementedError(f'Database was {database}! Only "90AA" and "100AA" are supported')
        self.seqix = IndexedFastaReader(
                f'{BASE_DIR}/GMSC10.{database}.fna'
                if path.exists(f'{BASE_DIR}/GMSC10.{database}.fna')
                else f'{BASE_DIR}/GMSC10.{database}.fna.xz'
                )
        self.database = database
        self.habitat = pd.read_table(f'{BASE_DIR}/GMSC10.{database}.habitat.index.tsv',
                                    index_col=0,
                                    header=None,
                                    names=['seq_ix', 'habitat']
                                    ).squeeze()
        self.habitat_ix = np.load(f'{BASE_DIR}/GMSC10.{database}.habitat.npy', mmap_mode='r')
        
        self.taxonomy = pd.read_table(f'{BASE_DIR}/GMSC10.{database}.taxonomy.index.tsv',
                                    index_col=0,
                                    header=None,
                                    names=['seq_ix', 'taxonomy']
                                    ).squeeze()
        self.taxonomy_ix = np.load(f'{BASE_DIR}/GMSC10.{database}.taxonomy.npy', mmap_mode='r')

        hqs = [line.strip() for line in lzma.open(f'{BASE_DIR}/GMSC10.90AA.high_quality.tsv.xz', 'rt')]
        hq_ixs = [int(hq.split('.')[2]) for hq in hqs]
        self.is_hq = np.zeros(len(self.habitat_ix), dtype=bool)
        self.is_hq[list(hq_ixs)] = True


    def get_seqinfo(self, seq_id):
        _,db,ix = seq_id.split('.')
        ix = int(ix)
        if db != self.database:
            raise IndexError(f'Only IDs for database "{self.database}" are accepted (got "{seq_id}"')

        nuc = self.seqix.get(ix).decode('ascii')
        return {
                "seq_id": seq_id,
                "nucleotide": nuc,
                "aminoacid": translate(nuc),
                'habitat': self.habitat.values[self.habitat_ix[ix]],
                'taxonomy': self.taxonomy.values[self.taxonomy_ix[ix]],
                }

    def seq_filter(self, hq_only : bool, habitat_q : List[str], taxonomy_q : str):
        if habitat_q:
            habitat_r = self.habitat.str.contains(habitat_q[0]).values
            for q in habitat_q[1:]:
                habitat_r &= self.habitat.str.contains(q).values
            matches = habitat_r[self.habitat_ix]
        else:
            matches = np.ones(len(self.habitat_ix), dtype=bool)
        if hq_only:
            matches &= self.is_hq
        if taxonomy_q is not None:
            match_taxonomy = self.taxonomy.str.contains(taxonomy_q).values[self.taxonomy_ix]
            matches &= match_taxonomy
        [ixs] = np.where(matches)
        # Highest numbers are best
        ixs = ixs[::-1]
        rs = []
        for i,ix in enumerate(ixs[:MAX_TOTAL_RESULTS]):
            seq_id = with_digits(f'GMSC10.{self.database}', ix)
            if i < MAX_THICK_RESULTS:
                rs.append(self.get_seqinfo(seq_id))
            else:
                rs.append({'seq_id': seq_id})
        return rs

class ClusterIx:
    def __init__(self):
        self.ix = np.load('gmsc-db/GMSC10.cluster.index.npy', mmap_mode='r')
        self.data = np.load('gmsc-db/GMSC10.cluster.data.npy', mmap_mode='r')

    def get_cluster_members(self, n : int):
        return self.data[self.ix[n]:self.ix[n+1]]
