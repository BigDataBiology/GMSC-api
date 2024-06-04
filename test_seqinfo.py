import seqinfo
def test_get_hits():
    import numpy as np
    arr = np.random.random(100_000)

    for t in [.01, .2, .4, .7, .9, 1.2]:
        matches = arr < t
        for n in [10, 100, 1000, 10_000, 100_000]:
            assert matches[seqinfo.get_hits(matches, n)].all()
