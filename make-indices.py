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

    starts = np.concatenate(starts)

    ofname = f'{ifname}.starts.npy'
    np.save(ofname, starts)
    return ofname

make_start_index('gmsc-db/GMSC10.100AA.fna.xz')
make_start_index('gmsc-db/GMSC10.90AA.fna')
