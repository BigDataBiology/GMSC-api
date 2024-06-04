import xz
import polars as pl
import pandas as pd

import gzip
import numpy as np
from os import path
from fna2faa_gmsc import translate
from typing import List, Optional


BASE_DIR = 'gmsc-db/'
INDEX_DIR = 'gmsc-db-index/'
MAX_THICK_RESULTS = 20
MAX_TOTAL_RESULTS = 1000

def with_digits(prefix, n):
    n = f'{n:09}'
    return f'{prefix}.{n[:3]}_{n[3:6]}_{n[6:9]}'

def get_hits(matches, max_results):
    n_matches = len(matches)
    # Highest numbers are best
    matches = matches[::-1]
    [ixs] = np.where(matches[:max_results])
    if len(ixs) == max_results:
        ixs *= -1
        ixs += n_matches - 1
        return ixs
    max_results -= len(ixs)
    chunks = [ixs]
    ch_size = max_results
    while max_results > 0:
        matches = matches[ch_size:]
        if not len(matches):
            break
        [ixs] = np.where(matches[:ch_size])
        ixs += n_matches - len(matches)
        chunks.append(ixs)
        max_results -= len(ixs)
        ch_size *= 2
    ixs = np.concatenate(chunks)
    ixs *= -1
    ixs += n_matches - 1
    return ixs


class IndexedFastaReader:
    def __init__(self, ifile):
        if ifile.endswith('.xz'):
            self.seqfile = xz.open(ifile, 'rb')
            ifile = ifile[:-len('.xz')]
        else:
            self.seqfile = open(ifile, 'rb')
        self.sindex = np.load(ifile.replace(BASE_DIR, INDEX_DIR) + '.starts.npy', mmap_mode='r')

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
        self.habitat = pd.read_table(f'{INDEX_DIR}/GMSC10.{database}.general_habitat.index.tsv',
                                    index_col=0,
                                    header=None,
                                    names=['seq_ix', 'habitat']
                                    ).squeeze()
        self.habitat_ix = np.load(f'{INDEX_DIR}/GMSC10.{database}.general_habitat.npy', mmap_mode='r')

        self.taxonomy = pd.read_table(f'{INDEX_DIR}/GMSC10.{database}.taxonomy.index.tsv',
                                    index_col=0,
                                    header=None,
                                    names=['seq_ix', 'taxonomy']
                                    ).squeeze()
        self.taxonomy_ix = np.load(f'{INDEX_DIR}/GMSC10.{database}.taxonomy.npy', mmap_mode='r')

        self.quality_metrics = None
        if path.exists(f'{INDEX_DIR}/GMSC10.{database}.quality_test.parquet'):
            self.quality_metrics = pl.read_parquet(f'{INDEX_DIR}/GMSC10.{database}.quality_test.parquet')

        self.is_hq = None
        if path.exists(f'{INDEX_DIR}/GMSC10.{database}.high_quality_ix.npy'):
            hq_ixs = np.load(f'{INDEX_DIR}/GMSC10.{database}.high_quality_ix.npy', mmap_mode='r')
            self.is_hq = np.zeros(len(self.habitat_ix), dtype=bool)
            self.is_hq[hq_ixs] = True


    def get_seqinfo(self, seq_id):
        _,db,ix = seq_id.split('.')
        ix = int(ix)
        if db != self.database:
            raise IndexError(f'Only IDs for database "{self.database}" are accepted (got "{seq_id}"')

        nuc = self.seqix.get(ix).decode('ascii')
        quality = None
        if self.quality_metrics is not None:
            quality = dict(zip(
                self.quality_metrics.columns,
                self.quality_metrics.row(ix)
                ))
        return {
                "seq_id": seq_id,
                "nucleotide": nuc,
                "aminoacid": translate(nuc),
                'habitat': self.habitat.values[self.habitat_ix[ix]],
                'taxonomy': self.taxonomy.values[self.taxonomy_ix[ix]],
                'quality': quality,
                }

    def seq_filter(self,
                   hq_only : bool,
                   habitat_q : List[str],
                   taxonomy_q : str,
                   *,
                   quality_antifam : Optional[bool] = None,
                   quality_terminal : Optional[bool] = None,
                   quality_rnacode : Optional[float] = None,
                   quality_metap : Optional[int] = None,
                   quality_metat : Optional[int] = None,
                   quality_riboseq : Optional[float] = None,
                   ):
        matches = None
        if habitat_q:
            habitat_r = self.habitat.str.contains(habitat_q[0]).values
            for q in habitat_q[1:]:
                habitat_r &= self.habitat.str.contains(q).values
            matches = habitat_r[self.habitat_ix]
        if hq_only:
            if self.is_hq is None:
                raise ValueError('High quality information not loaded')
            if matches is None:
                matches = self.is_hq.copy()
            else:
                matches &= self.is_hq
        if taxonomy_q:
            match_taxonomy = self.taxonomy.str.contains(taxonomy_q).values[self.taxonomy_ix]
            if matches is None:
                matches = match_taxonomy
            else:
                matches &= match_taxonomy
        advanced_conditions = []
        if quality_antifam is not None:
            if quality_antifam:
                advanced_conditions.append(pl.col('antifam'))
            else:
                advanced_conditions.append(pl.col('antifam').not_())
        if quality_terminal is not None:
            if quality_terminal:
                advanced_conditions.append(pl.col('terminal'))
            else:
                advanced_conditions.append(pl.col('terminal').not_())
        if quality_rnacode is not None:
            advanced_conditions.append(pl.col('rnacode') <= quality_rnacode)
        if quality_metap is not None:
            advanced_conditions.append(pl.col('metap') >= quality_metap)
        if quality_metat is not None:
            advanced_conditions.append(pl.col('metat') >= quality_metat)
        if quality_riboseq is not None:
            advanced_conditions.append(pl.col('riboseq') >= quality_riboseq)
        if advanced_conditions:
            if self.quality_metrics is None:
                raise ValueError('Quality metrics not loaded')
            if len(advanced_conditions) == 1:
                [advanced_conditions] = advanced_conditions
            else:
                advanced_conditions = advanced_conditions[0].and_(*advanced_conditions[1:])
            sel = self.quality_metrics.select(advanced_conditions.alias('matched'))
            if matches is None:
                matches = sel['matched'].to_numpy()
            else:
                matches &= sel['matched'].to_numpy()

        if matches is not None:
            ixs = get_hits(matches, MAX_TOTAL_RESULTS)
        else:
            ixs = np.arange(len(self.habitat_ix)-1, len(self.habitat_ix)-MAX_TOTAL_RESULTS-1, -1)

        rs = []
        for i,ix in enumerate(ixs):
            seq_id = with_digits(f'GMSC10.{self.database}', ix)
            if i < MAX_THICK_RESULTS:
                rs.append(self.get_seqinfo(seq_id))
            else:
                rs.append({'seq_id': seq_id})
        return rs

class ClusterIx:
    def __init__(self):
        self.ix = np.load(f'{INDEX_DIR}/GMSC10.cluster.index.npy', mmap_mode='r')
        self.data = np.load(f'{INDEX_DIR}/GMSC10.cluster.data.npy', mmap_mode='r')

    def get_cluster_members(self, n : int):
        return self.data[self.ix[n]:self.ix[n+1]]
