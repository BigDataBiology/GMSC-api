from os import path
from time import sleep

DB_DIR = 'gmsc-db'

def parse_gmsc_mapper_results(basedir):
    import pandas as pd
    try:
        alignment = pd.read_table(f'{basedir}/alignment.out.smorfs.tsv', header=None)
        habitat = pd.read_table(f'{basedir}/habitat.out.smorfs.tsv', index_col=0)
        quality = pd.read_table(f'{basedir}/quality.out.smorfs.tsv', index_col=0)
        taxonomy = pd.read_table(f'{basedir}/taxonomy.out.smorfs.tsv', index_col=0)
    except pd.errors.EmptyDataError:
        return {}
    meta = pd.concat([habitat, quality, taxonomy], axis=1)
    alignment.columns = 'qseqid,sseqid,full_qseq,full_sseq,qlen,slen,length,qstart,qend,sstart,send,bitscore,pident,evalue,qcovhsp,scovhsp'.split(',')

    result = meta.to_dict('index')

    for k in result.keys():
        loc = alignment.query('qseqid == @k')
        result[k]['aminoacid'] = loc.head(1).full_qseq.values[0]
        result[k]['hits'] = \
                loc[['sseqid', 'evalue', 'pident']].rename(columns=
                    {'sseqid': 'id',
                    'pident': 'identity'},).to_dict('records')
    return result


def save_search_result(r, sid):
    import json
    with open(f'search-results/{sid}.json', 'wt') as out:
        json.dump(r, out)

def load_search_results():
    from glob import glob
    import json
    rs = glob('search-results/*.json')
    rs.sort(key = lambda k: int(k.split('/')[1].split('-')[0]))
    res = {}
    for r in rs:
        k = r.split('/')[1].split('.')[0]
        res[k] = json.load(open(r))
    return res


def do_search(seqdata, sid, is_contigs, nr_threads):
    import tempfile
    import subprocess
    from time import time
    with tempfile.TemporaryDirectory() as tdir:
        if not path.exists(DB_DIR):
            sleep(10)
            return parse_gmsc_mapper_results('./demo_gmsc_mapper_output')
        fname = path.join(tdir, ('seqs.fna' if is_contigs else 'seqs.faa'))
        with open(fname, "w") as f:
            f.write(seqdata)
        sleep(1)
        subprocess.check_call(
                ['gmsc-mapper',
                 ('--input' if is_contigs else '--aa-genes'), fname,
                 '-o', path.join(tdir, "output"),
                 '--threads', str(nr_threads),
                 '--db', f'{DB_DIR}/GMSC10.90AA.diamonddb.dmnd',
                 '--habitat', f'{DB_DIR}/GMSC10.90AA.habitat.npy',
                 '--habitat-index', f'{DB_DIR}/GMSC10.90AA.habitat.index.tsv',
                 '--quality', f'{DB_DIR}/GMSC10.90AA.high_quality.tsv.xz',
                 '--taxonomy', f'{DB_DIR}/GMSC10.90AA.taxonomy.npy',
                 '--taxonomy-index', f'{DB_DIR}/GMSC10.90AA.taxonomy.index.tsv',
                 ],
                )
        r = parse_gmsc_mapper_results(path.join(tdir, "output"))
        save_search_result(r, sid)
        return r

