try:
    from jug import TaskGenerator
except ImportError:
    import sys
    sys.stderr.write('Jug is necessary to run this script\n')
    sys.stderr.write('Please install it.\n')
    raise


@TaskGenerator
def make_start_index(ifname):
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

    ofname = f'{ifname}.starts.npy'
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
def make_cluster_index(sizes):
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
    np.save('gmsc-db/GMSC10.cluster.index.npy', ix)
    np.save('gmsc-db/GMSC10.cluster.data.npy', data)

make_start_index('gmsc-db/GMSC10.100AA.fna.xz')
make_start_index('gmsc-db/GMSC10.90AA.fna')


sizes = get_cluster_sizes()
make_cluster_index(sizes)

