from jug import Task, TaskGenerator, barrier
from jug.utils import jug_execute
import os

INDEX_DIR = 'gmsc-db'
HASHES = {
    'GMSC10.cluster.sorted2.tsv.xz': '6de8fe07f523b48d59d133fb4ffabc14',
    'GMSC10.100AA.fna.xz': '18a9e27e976f27082e4d17450b15961d',
    'GMSC10.100AA.annotation.tsv.xz': '50577cb91086e38d59f8053d2fbc65c5',
    'GMSC10.90AA.annotation.tsv.xz': 'b2cb758cb1892909b68b8463350e8be8',
    'GMSC10.90AA.txt.xz': '83470ca958d84b4df374b1f9d0338638',
    'fna2faa_gmsc.py': '30cf3d9f82ecf9e4d89cfb476d88b60d',
    }


def md5_file(f):
    import hashlib
    h = hashlib.md5()
    with open(f, 'rb') as ifile:
        while ch := ifile.read(1024*1024):
            h.update(ch)
    return h.hexdigest()

@Task
def make_index_dir():
    os.makedirs(INDEX_DIR, exist_ok=True)

barrier()

@TaskGenerator
def download_file_if_needed(f):
    from os import path
    import requests
    target = f'{INDEX_DIR}/{f}'
    if path.exists(target):
        print(f'File exists. Checking hash...')
        if md5_file(target) == HASHES[f]:
            print('Hash ok')
            return target
        raise IOError(f'Unexpected hash for {target}')

    url = f'https://zenodo.org/records/7944371/files/{f}?download=1'
    r = requests.get(url, allow_redirects=True, stream=True)
    with open(target, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    return target

@TaskGenerator
def create_90aa(fna_file, sel90):
    import lzma
    import pandas as pd
    import numpy as np
    import subprocess

    to90 = np.empty(1_000_000_000, np.int64)
    to90.fill(-1)


    def with_digits_90aa(n):
        n = f'{n:09}'
        return f'GMSC10.90AA.{n[:3]}_{n[3:6]}_{n[6:9]}'

    ix = 0
    for ch in pd.read_table(sel90, chunksize=100_000, header=None):
        ch = ch[0].str.split('.').str[2].map(int).reset_index()
        ch = ch.values
        to90[ch.T[1]] = ch.T[0]
        ix += 1

    ifile = lzma.open(fna_file, 'rt')
    tmp_file = 'tmp_quasi_fa'
    with open(tmp_file, 'wt') as ofile:
        while True:
            header = ifile.readline()
            if not header:
                break
            seq = ifile.readline()

            ix = int(header[14:-1],10)
            if (p := to90[ix]) != -1:
                ofile.write(f'>{with_digits_90aa(p)}\t{seq}')

    del to90
    subprocess.check_call(['sort', '-S', '12G', tmp_file, '-o', tmp_file+'.sorted'])
    os.unlink(tmp_file)

    oname = f'{INDEX_DIR}/GMSC10.90AA.fna'
    with open(oname, 'wt') as ofile:
        for line in open(tmp_file+'.sorted'):
            line = line.replace('\t', '\n')
            ofile.write(line)
    os.unlink(tmp_file+'.sorted')
    return oname

@TaskGenerator
def fna2faa(fna):
    import subprocess
    assert fna.endswith('.fna')
    faa = fna[:-len('.fna')]+'.faa'
    with open(fna, 'rb') as ifile:
        with open(faa, 'wb') as ofile:
            subprocess.check_call(
                    ['python', f'{INDEX_DIR}/fna2faa_gmsc.py'],
                    stdin=ifile,
                    stdout=ofile,
                    )
    return faa

download_file_if_needed('GMSC10.cluster.sorted2.tsv.xz')
download_file_if_needed('fna2faa_gmsc.py')
download_file_if_needed('GMSC10.100AA.annotation.tsv.xz')
download_file_if_needed( 'GMSC10.90AA.annotation.tsv.xz')
fna90 = create_90aa(
        download_file_if_needed('GMSC10.100AA.fna.xz'),
        download_file_if_needed('GMSC10.90AA.txt.xz'),
        )
jug_execute(['xz', '--threads=0', '--keep', fna90])
faa90 = fna2faa(fna90)
jug_execute(['diamond', 'makedb', '--db', f'{INDEX_DIR}/GMSC10.90AA.diamonddb.dmnd', '--in', faa90])
