from flask import Flask, request, jsonify
from fna2faa_gmsc import translate
from datetime import datetime
from os import path
import sqlite3
DB_PATH = 'gmsc-db/gmsc10hq.sqlite3'
if path.exists(DB_PATH):
    con = sqlite3.connect(DB_PATH, check_same_thread=False)
else:
    import sys
    sys.stderr.write(f'WARNING: Database file {DB_PATH} not found\n')
    sys.stderr.write(f'WARNING: Using demo database\n')
    con = sqlite3.connect('gmsc10-demo.sqlite3', check_same_thread=False)

app = Flask('GMSC')

sequences = [
    ("GMSC10.100AA.000_000_000", "GTGGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCGGCAGGTGGGGCAAGAGGAGAGAATTTTGACGAGAATAAAATAGACGCAGAAAGAGAAGCAGGAGTAGACGTACAAGTCGATCGCGGGGTGCTGCTGCTGTTGCTGTTGATCCTACTGCTGCTATTGCTGCTGCTGCTGCTGCTGTTGGTGCTGGTGACGCTGGCCGCCGTGCTCCCTTGTCGGGATAAGGGCGGGGATTGA", "water associated", ""),

    ("GMSC10.100AA.000_000_001", "ATGGCTGCTGCAGCAGCAGCTGCTGCTGCTGCTGCCGCCGCCTTTGTAGTTTTCGTAGGTTTTCCATCGTCATCTTTCTTGTCTGATGATTTTCTAACATTCTTGTCAGAATCCGGAAAGGATGGCAGGGCGGGTTTCCGGCTCCGAGGAAATCGAGGTGGCGGAGCTCTAGTTTCTTGGTATCCTGCCAATATCTGTTGGGCATCAATAGCAGCATCTACTAGATTGGAAATGTACCGGGCTGTCAATTCTGAGAGGAGTTGA", "marine", ""),
    ("GMSC10.100AA.000_000_002", "GTGGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCCAAAGGGAGTGCTGTGCGCAATATATTTTCGACGGGCATGCTTGAGACACTGCGAGCACTGGCTGACGGTGTTGCAGCATCTGATATACTTCCAAGAGTACTGGCAGTACGTCATCGTTCAATACTGCCCATACCGCGAGTGAGTATACTCGCAGTATTTCGGATTATGTATGTCTGTACTGCCTGA", "marine", ""),
    ("GMSC10.100AA.000_000_003", "ATGGCCGCCGCCGCAGCCGCCGCCGCCGCCGCCGCCGCCGTCGGAACCTCCTTGGACACCGGGCCGGACCCGCTGATGGCGATGCTCGCCGGCGATTTCCAGCCCGCGGCGGGGTTCAAGTCGAGCGACTCGGCGCCGGAGCTCGTCGAGGTTGTCAAGATTCACTCTCCGGTTGTCGCGACGGGGGAGAACGGCGCCGCGAGCGCGATGGATGTCCTGACGGAGCCGGCTTGA", "marine", ""),
    ("GMSC10.100AA.000_000_004", "ATGGCCGCCGCCGCCGCCGCCGCCGCCGCAGCCGCCGCCGTCGGAACCTCCGTGGACACCGGGCCGGACCCGCTGATGGCGATGCTCGCCGGCGATTTCCAGCCCGCTGCGGGGTTCAAGTCGAGCGACTCGGCGCCGGAGCTCGTCGAGGTTGTCAAGATTCACTCGCCGGTTGTCGCGACGGGGGAGAACGGCGCCGCGAGCGCGATGGATGTCCTGACGGAGCCGGCTTGA", "soil", "d__Bacteria;p__Actinobacteriota;c__Actinomycetia;o__Streptosporangiales;f__Streptosporangiaceae"),
    ("GMSC10.100AA.000_000_005", "GTGGCGGCTGCGGCCGCAGCTGCCGCAGCGGCGGCCTGGACCTGGGCTTTTTCCTTCAGGTCGGCCATCTTGTTGTTCATGTCCTGCAGTTTTTCATCGGCGGCCTGGCGAATGCTGTCCAGATCGACCTTTTCAGAAAACTCGCACCAGGTGAGCACACCCACCATCAGCGGCAGCAAGAAGACGAAAGCCGAGGCCAGCGCAAACACCAGGCTGAAGCCCACGCCTGCCCCCATCATGTGGCCCGAGGCGCCGCCCATCATGATGGCCATGGGGTTAAAGCCACCCATTCCATAA", "soil", "d__Bacteria;p__Actinobacteriota;c__Actinomycetia"),
]


@app.get("/v1/seq-info/<seq_id>")
def get_seq_info(seq_id):
    tokens = seq_id.split(".")
    if len(tokens) != 3:
        return {"error": "Invalid sequence ID (only 'GMSC10' database supported)"}, 400
    (db, sample_type, seq_id) = tokens
    if db != 'GMSC10':
        return {"error": "Invalid sequence ID"}, 400
    if sample_type not in ('100AA', '90AA'):
        return {"error": "Invalid sequence ID (seq type must be '90AA' or '100AA')"}, 400
    seq_id = int(seq_id, 10)
    if seq_id >= len(sequences):
        return {"error": "Invalid sequence ID (too large)"}, 400
    (seq_id, nuc, habitat, taxonomy) = sequences[seq_id]
    aa = translate(nuc)
    return jsonify({
        "seq_id": seq_id,
        "nucleotide": nuc,
        "aminoacid": aa,
        "habitat": habitat,
        "taxonomy": taxonomy,
        })


def with_digits(prefix, n):
    n = f'{n:09}'
    return f'{prefix}.{n[:3]}_{n[3:6]}_{n[6:9]}'

def parse_bool(s : str):
    s = s.lower()
    if s in ('true', '1', 'yes'):
        return True
    if s in ('false', '0', 'no'):
        return False
    return None

@app.route('/v1/seq-filter/', methods=['POST'])
def get_seq_filter():
    cur = con.cursor()
    hq_only = request.form.get('hq_only', False)
    if not parse_bool(hq_only):
        return {"status": "error", "msg": "Only HQ searching supported (at the moment)"}, 400
    habitat = request.form.get('habitat')
    if habitat is None:
        return {"status": "error", "msg": "Invalid habitat value"}, 400
    q = "SELECT * FROM smorf90aa " \
                    "WHERE " \
                        "instr(habitat, ?) > 0 "
    qargs = [habitat]
    taxonomy = request.form.get('taxonomy')
    if taxonomy is not None :
        if '__' not in taxonomy:
            return {"status": "error", "msg": "Invalid taxonomy value"}, 400
        q += "AND instr(taxonomy, ?) > 0 "
        qargs.append(taxonomy)
    q += "LIMIT 1001"
    results = [
        { "seq_id": with_digits("GMSC10.90AA", i),
           "habitat": habitat,
           "taxonomy": taxonomy}
           for (i, habitat, taxonomy) in
                cur.execute(q, qargs).fetchall()]
    return jsonify({
        "status": "Ok",
        "results": results,
        })

searches = {}

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

next_search_id = SearchIDGenerator()

@app.route('/internal/seq-search/', methods=['POST'])
def seq_search():
    sid = next_search_id.get_next_id()
    searches[sid] = datetime.now()
    return jsonify({
        "search_id": sid,
        "status": "Ok",
        })

@app.get('/internal/seq-search/<search_id>')
def seq_search_results(search_id):
    if search_id not in searches:
        return {"error": "Invalid search ID"}, 400
    sdata = searches[search_id]
    if (datetime.now() - sdata).seconds < 10:
        return {"search_id": search_id, "status": "Running"}
    if (datetime.now() - sdata).seconds > 120:
        return {"search_id": search_id, "status": "Expired"}
    return jsonify({
        "search_id": search_id,
        "status": "Done",
        "results": [
            {
                "query_id": "query_1",
                "aminoacid": "MHEDVIQFARNEVWSLV....",
                "taxonomy": "s__Bacteroides_vulgatus",
                "habitat": "human gut",
                "hits": [
                    { "id": "GMSC10.xxAA.xxx_xxx_xxxx",
                      "e_value": "2.1e-23",
                      "aminoacid": "MHEELIQFARNEV...",
                      "identity": "98.4"
                    },
                ]
            },
        ]})
