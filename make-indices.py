try:
    from jug import TaskGenerator
except ImportError:
    import sys
    sys.stderr.write('Jug is necessary to run this script\n')
    sys.stderr.write('Please install it.\n')
    raise

from os import path

@TaskGenerator
def make_start_index(ifname, index_dir):
    import numpy as np
    import lzma
    CHUNK_SIZE = 1_000_000
    LT = ord(b'>')

    if ifname.endswith('.xz'):
        ifname = ifname[:-len('.xz')]
        op = lambda : lzma.open(ifname + '.xz', 'rb')
    else:
        op = lambda : open(ifname, 'rb')
    with op() as ifile:
        start = 0
        starts = []
        while ch := ifile.read(CHUNK_SIZE):
            ch = np.frombuffer(ch, np.uint8)
            [cur] = np.where(ch==LT)
            cur += start
            starts.append(cur)
            start += len(ch)

    starts.append([start])
    starts = np.concatenate(starts)

    ofname = f'{index_dir}/{path.basename(ifname)}.starts.npy'
    np.save(ofname, starts)
    return ofname


def get_ix(n):
    _,_,n = n.split('.')
    return int(n)


@TaskGenerator
def get_cluster_sizes():
    import lzma
    total_n = 0
    max_90aa = -1
    prev = -1
    with lzma.open('gmsc-db/GMSC10.cluster.sorted2.tsv.xz', 'rt') as f:
        for line in f:
            total_n += 1
            _,n = line.split()
            n = get_ix(n)
            if n > max_90aa:
                max_90aa = n
            if n < prev:
                raise ValueError(f'Line {total_n} not sorted')
            if n > (prev +1):
                raise ValueError(f'Line {total_n} skipped')
            prev = n
    return total_n, max_90aa

@TaskGenerator
def make_cluster_index(sizes, index_dir):
    import lzma
    import numpy as np
    total_n, max90aa = sizes
    ix = np.zeros(max90aa + 2, dtype=np.uint64)
    data = np.zeros(total_n, dtype=np.uint64)
    prev = -1
    with lzma.open('gmsc-db/GMSC10.cluster.sorted2.tsv.xz', 'rt') as f:
        for cur100,line in enumerate(f):
            ix100,ix90 = line.split()
            ix100 = get_ix(ix100)
            ix90 = get_ix(ix90)
            if ix90 != prev:
                ix[ix90] = cur100
                prev = ix90
            data[cur100] = ix100
    ix[-1] = len(data)
    np.save(f'{index_dir}/GMSC10.cluster.index.npy', ix)
    np.save(f'{index_dir}/GMSC10.cluster.data.npy', data)

def convert_to_npy(infile):
    import numpy as np
    outfile = infile[:-len('.idx.tsv')] + '.npy'
    df = np.loadtxt(infile, dtype=int, usecols=1)
    np.save(outfile, df)

@TaskGenerator
def create_tax_index(infile, index_dir):
    import lzma
    import numpy as np

    tax_set = set()

    nr_lines = 0
    with lzma.open(infile,'rt') as f:
        for line in f:
            nr_lines += 1
            linelist = line.strip().split('\t')
            if len(linelist) == 2:
                tax_set.add(linelist[1])
            else:
                tax_set.add('Unknown')

    tax_order = sorted(list(tax_set))
    tax_dict = {}

    assert infile.endswith('.tsv.xz')
    basename = path.basename(infile)[:-len('.tsv.xz')]
    outfile1 = f'{index_dir}/{basename}.index.tsv'
    outfile2 = f'{index_dir}/{basename}.npy'

    with open(outfile1,'wt') as out:
        for n,item in enumerate(tax_order):
            tax_dict[item] = n
            out.write(f'{n}\t{item}\n')

    odata = np.zeros(nr_lines, int)
    with lzma.open(infile,'rt') as f:
        for ix,line in enumerate(f):
            tokens = line.strip().split('\t')
            if len(tokens) == 2:
                odata[ix] = tax_dict[tokens[1]]
    np.save(outfile2, odata)
    return outfile2

@TaskGenerator
def create_habitat_index(infile, index_dir):
    import lzma
    habitat_set = set()

    with lzma.open(infile,'rt') as f:
        for line in f:
            smorf,habitat = line.strip().split('\t')
            habitat_set.add(habitat)

    habitat_order = sorted(list(habitat_set))
    habitat_dict = {}

    assert infile.endswith('.tsv.xz')
    basename = path.basename(infile)[:-len('.tsv.xz')]
    outfile1 = f'{index_dir}/{basename}.index.tsv'
    outfile2 = f'{index_dir}/{basename}.idx.tsv'

    with open(outfile1,'wt') as out:
        for n,item in enumerate(habitat_order):
            habitat_dict[item] = n
            out.write(f'{n}\t{item}\n')

    with open(outfile2,'wt') as out:
        with lzma.open(infile,'rt') as f:
            for line in f:
                smorf,habitat = line.strip().split('\t')
                out.write(f'{smorf}\t{habitat_dict[habitat]}\n')
    convert_to_npy(outfile2)
    return outfile2


INDEX_DIRECTORY = 'gmsc-db-index'

make_start_index('gmsc-db/GMSC10.100AA.fna.xz', INDEX_DIRECTORY)
make_start_index('gmsc-db/GMSC10.90AA.fna.xz', INDEX_DIRECTORY)


sizes = get_cluster_sizes()
make_cluster_index(sizes, INDEX_DIRECTORY)


create_tax_index('gmsc-db/GMSC10.100AA.taxonomy.tsv.xz', INDEX_DIRECTORY)
create_tax_index('gmsc-db/GMSC10.90AA.taxonomy.tsv.xz', INDEX_DIRECTORY)

create_habitat_index('gmsc-db/GMSC10.100AA.general_habitat.tsv.xz', INDEX_DIRECTORY)
create_habitat_index('gmsc-db/GMSC10.90AA.general_habitat.tsv.xz', INDEX_DIRECTORY)
