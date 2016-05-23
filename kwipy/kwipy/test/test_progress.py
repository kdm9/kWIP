from kwipy.progress import (
    human_readable,
    ProgressLogger,
)

def test_human_readable():
    units = ['', 'K', 'M']

    for i, u in enumerate(units):
        i = 1000 ** i
        assert human_readable(i) == "1{}".format(u)

    assert human_readable(1e9, maxsuffix='M') == '1000M'
    assert human_readable(1e9, maxsuffix='G') == '1G'
