from flask import Flask, request, jsonify
from datetime import datetime
from os import path, environ
import threading
from time import sleep
from concurrent.futures import ThreadPoolExecutor
from collections import namedtuple

from seqinfo import SeqInfo, ClusterIx, with_digits
from search import do_search, load_search_results

NR_THREADS_GMSC_MAPPER = 2

DB_DIR = 'gmsc-db'
if path.exists(DB_DIR):
    seqinfo90  = SeqInfo( '90AA')
    seqinfo100 = SeqInfo('100AA')
    clusterinfo = ClusterIx()
    IS_DEMO = False
else:
    import sys
    sys.stderr.write(f'WARNING: Database directory {DB_DIR} not found\n')
    sys.stderr.write(f'WARNING: Using demo database\n')
    IS_DEMO = True

app = Flask('GMSC')


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
        if cluster_level == '90AA':
            rs.append(seqinfo90.get_seqinfo(seq_id))
        elif cluster_level == '100AA':
            rs.append(seqinfo100.get_seqinfo(seq_id))
        else:
            return {"error": "Invalid sequence ID: does not match GMSC10 ID format"}, 400
    return rs


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
        habitat = []
    else:
        habitat = habitat.split(",")
    taxonomy = request.form.get('taxonomy')
    if taxonomy is None:
        taxonomy = ""
    results = seqinfo90.seq_filter(parse_bool(hq_only), habitat, taxonomy)
    return jsonify({
        "status": "Ok",
        "results": results,
        })


@app.get('/v1/cluster-info/<cluster_id>')
def get_cluster_info(cluster_id):
    tokens = cluster_id.split(".")
    if len(tokens) != 3:
        return {"error": "Invalid sequence ID (only 'GMSC10' database supported)"}, 400
    (db, cluster_level, seq_ix) = tokens
    if db != 'GMSC10':
        return {"error": "Invalid sequence ID"}, 400
    if cluster_level != '90AA':
        return {"error": "Invalid sequence ID: only 90AA identifiers can be used for cluster-info"}, 400
    members = clusterinfo.get_cluster_members(int(seq_ix))
    rs = []
    for m in members:
        m = with_digits('GMSC10.100AA', m)
        if len(rs) < 20:
            rs.append(seqinfo100.get_seqinfo(m))
        else:
            rs.append({'seq_id': m})
    return {
            'status': 'Ok',
            'cluster': rs,
            }



class SearchIDGenerator:
    def __init__(self, start):
        self.next_id = start
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

def identity(x):
    return x

search_lock = threading.Lock()
searcher = ThreadPoolExecutor(2)
SearchObject = namedtuple("SearchObject", ["start_time", "future"])
searches = {
        k:SearchObject(datetime.now(), searcher.submit(identity, v)) for k,v in load_search_results().items()
        }

next_search_id = SearchIDGenerator(len(searches))

@app.post('/internal/seq-search/')
def seq_search():
    now = datetime.now()
    seqdata = request.form.get('sequence_faa')
    is_contigs = parse_bool(request.form.get('is_contigs'))
    if seqdata is None:
        return {"error": "Missing sequence_faa parameter"}, 400
    with search_lock:
        sid = next_search_id.get_next_id()
        searches[sid] = SearchObject(now, searcher.submit(do_search, seqdata, sid, is_contigs, NR_THREADS_GMSC_MAPPER))
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
        status = 'Running' if sdata.running() else 'Queued'
        sleep(1)
        return {
                "search_id": search_id,
                "status": status,
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
        return {"error": "No secret set"}, 500
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

