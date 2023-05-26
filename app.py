from flask import Flask, request, jsonify
from datetime import datetime
from os import path, environ
import threading
from time import sleep
from concurrent.futures import ThreadPoolExecutor
from collections import namedtuple
from flask_cors import CORS

from seqinfo import SeqInfo

NR_THREADS_GMSC_MAPPER = 2

DB_DIR = 'gmsc-db'
if path.exists(DB_DIR):
    seqinfo90  = SeqInfo( '90AA')
    seqinfo100 = SeqInfo('100AA')
    IS_DEMO = False
else:
    import sys
    sys.stderr.write(f'WARNING: Database directory {DB_DIR} not found\n')
    sys.stderr.write(f'WARNING: Using demo database\n')
    IS_DEMO = True

app = Flask('GMSC')
CORS(app)


@app.get("/v1/seq-info/<seq_id>")
def get_seq_info(seq_id):
    tokens = seq_id.split(".")
    if len(tokens) != 3:
        return {"error": "Invalid sequence ID (only 'GMSC10' database supported)"}, 400
    (db, cluster_level, seq_ix) = tokens
    if db != 'GMSC10':
        return {"error": "Invalid sequence ID"}, 400
    if cluster_level not in ('100AA', '90AA'):
        return {"error": "Invalid sequence ID (seq type must be '90AA' or '100AA')"}, 400
    if IS_DEMO:
        from demo import get_demo_seqinfo
        return get_demo_seqinfo(int(seq_ix))
    return (seqinfo90 if cluster_level == '90AA' else seqinfo100).get_seqinfo(seq_id)


@app.post("/v1/seq-info-multi/")
def get_seq_info_multi():
    entries = request.json.get('seq_ids')
    if entries is None:
        return {"error": "Missing seq_ids parameter"}, 400
    if len(entries) > 100:
        return {"error": "Too many seq_ids"}, 400
    rs = []
    for seq_id in entries:
        tokens = seq_id.split(".")
        if len(tokens) != 3:
            return {"error": "Invalid sequence ID (only 'GMSC10' database supported)"}, 400
        (db, cluster_level, seq_ix) = tokens
        if db != 'GMSC10':
            return {"error": "Invalid sequence ID"}, 400
        if cluster_level != '90AA':
            return {"error": "Invalid sequence ID: only 90AA identifiers can be used for seq-info-multi"}, 400
        rs.append(seqinfo90.get_seqinfo(e))


def parse_bool(s : str):
    if s is None:
        return False
    if s in (True, False):
        return s
    s = s.lower()
    if s in ('true', '1', 'yes'):
        return True
    if s in ('false', '0', 'no'):
        return False
    return None

@app.post('/v1/seq-filter/')
def get_seq_filter():
    hq_only = request.form.get('hq_only', False)
    habitat = request.form.get('habitat')
    if habitat is None:
        return {"status": "error", "msg": "Invalid habitat value"}, 400
    taxonomy = request.form.get('taxonomy')
    if taxonomy is not None :
        if '__' not in taxonomy:
            return {"status": "error", "msg": "Invalid taxonomy value"}, 400
    results = seqinfo90.seq_filter(parse_bool(hq_only), habitat, taxonomy)
    return jsonify({
        "status": "Ok",
        "results": results,
        })

class SearchIDGenerator:
    def __init__(self):
        self.next_id = 0
    def get_next_id(self):
        import random
        from string import ascii_lowercase
        randstr = ''.join(random.choice(ascii_lowercase) for i in range(4))
        self.next_id += 1
        return f"{self.next_id}-{randstr}"

    @classmethod
    def get_index(cls, search_id : str) -> int:
        ix, _ = search_id.split("-")
        return int(ix)

    def get_cur_index(self) -> int:
        return self.next_id

searches = {}
next_search_id = SearchIDGenerator()
search_lock = threading.Lock()
searcher = ThreadPoolExecutor(2)
SearchObject = namedtuple("SearchObject", ["start_time", "future"])

def parse_gmsc_mapper_results(basedir):
    import pandas as pd
    alignment = pd.read_table(f'{basedir}/alignment.out.smorfs.tsv', header=None)
    habitat = pd.read_table(f'{basedir}/habitat.out.smorfs.tsv', index_col=0)
    quality = pd.read_table(f'{basedir}/quality.out.smorfs.tsv', index_col=0)
    taxonomy = pd.read_table(f'{basedir}/taxonomy.out.smorfs.tsv', index_col=0)
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

def do_search(seqdata):
    import tempfile
    import subprocess
    with tempfile.TemporaryDirectory() as tdir:
        if not path.exists(DB_DIR):
            sleep(10)
            return parse_gmsc_mapper_results('./demo_gmsc_mapper_output')
        with open(path.join(tdir, "seqs.faa"), "w") as f:
            f.write(seqdata)
        sleep(1)
        subprocess.check_call(
                ['gmsc-mapper',
                 '--aa-genes', path.join(tdir, "seqs.faa"),
                 '-o', path.join(tdir, "output"),
                 '--threads', str(NR_THREADS_GMSC_MAPPER),
                 '--db', f'{DB_DIR}/90AA_GMSC.dmnd',
                 '--habitat', f'{DB_DIR}/90AA_habitat.npy',
                 '--habitat-index', f'{DB_DIR}/90AA_ref_multi_general_habitat_index.tsv',
                 '--quality', f'{DB_DIR}/90AA_highquality.txt.gz',
                 '--taxonomy', f'{DB_DIR}/90AA_taxonomy.npy',
                 '--taxonomy-index', f'{DB_DIR}/90AA_ref_taxonomy_index.tsv',
                 ],
                )

        return parse_gmsc_mapper_results(path.join(tdir, "output"))


@app.post('/internal/seq-search/')
def seq_search():
    now = datetime.now()
    seqdata = request.form.get('sequence_faa')
    with search_lock:
        sid = next_search_id.get_next_id()
        searches[sid] = SearchObject(now, searcher.submit(do_search, seqdata))
    return jsonify({
        "search_id": sid,
        "status": "Ok",
        })

@app.get('/internal/seq-search/<search_id>')
def seq_search_results(search_id):
    if search_id not in searches:
        return {"error": "Invalid search ID"}, 400
    sdata = searches[search_id].future
    if not sdata.done():
        sleep(1)
        return {
                "search_id": search_id,
                "status": 'Running' if sdata.running() else 'Queued'
                }
    r = sdata.result()
    return jsonify({
        "search_id": search_id,
        "status": "Done",
        "results": r
        })

@app.post('/internal/seq-search-list/')
def seq_search_list():
    secret = environ.get('GMSC_API_INTERNAL_PWD', None)
    pwd = request.form.get('pwd')
    if secret is None:
        return {"error": "Invalid search ID"}, 500
    if pwd != secret:
        return {"error": "Wrong password"}, 403
    def status_for(f):
        if f.done():
            return "Done"
        if f.running():
            return "Running"
        return "Queued"
    return {
            "status": "Ok",
            "searches": [
                {"search_id": k,
                "status": status_for(v.future),
                } for k,v in searches.items()],
            }

